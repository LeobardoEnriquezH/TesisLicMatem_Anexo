// =====================================================================
//  Proceso Gaussiano (GP) SEMIPARAMETRICO para la media de respuesta
//      mu(A, X) = m(A, X) + f(A, X),   f ~ GP(0, C)
//  donde m(A, X) = gamma0 + tau*A + X'beta es una MEDIA LINEAL y f es el GP.
//
//  Kernel (identico a la ecuacion gp_kernel de la tesis):
//      C_ij = eta * exp( -rho * || z_i - z_j ||^2 )  +  0.01 * delta_ij
//  con la convencion EXACTA del documento:
//      * rho  : tasa de DECAIMIENTO.  rho mayor  =>  la correlacion cae MAS
//               rapido con la distancia (NO es la length-scale de cov_exp_quad).
//      * eta  : amplitud = varianza marginal del proceso   (eta = alpha^2).
//      * phi  : varianza del ruido de observacion           (phi = sigma^2).
//      * 0.01 : jitter de la teoria, incluido EXPLICITAMENTE en la covarianza.
//
//  POR QUE LA MEDIA LINEAL (anti-encogimiento hacia cero):
//  Con media 0, alli donde los datos escasean o el kernel acopla los brazos de
//  tratamiento, el GP revierte a 0 y el ATE se atenua artificialmente. Con media
//  lineal, el GP modela DESVIACIONES respecto a la regresion lineal: la reversion
//  ocurre hacia un efecto NO NULO (la parte lineal), no hacia cero. El kernel sigue
//  siendo EXACTAMENTE el de la tesis. El ATE = parte lineal (via tau) + desviacion
//  flexible del GP. Se mantiene A tambien en el kernel para captar heterogeneidad /
//  no linealidad del tratamiento sin colapsar a cero.
//
//  POR QUE EL KERNEL A MANO (squared_distance) Y NO cov_exp_quad():
//  cov_exp_quad(z, alpha, rho) = alpha^2 * exp( -||.||^2 / (2 rho^2) ): su "rho"
//  es una length-scale (rho mayor => decaimiento MAS LENTO, lo contrario a la
//  teoria) y su varianza marginal es alpha^2, no eta. Escribirlo a mano hace que
//  (eta, rho) coincidan LITERALMENTE con la teoria.
//
//  NOTA DE VERSION DE STAN
//  -----------------------
//  Usa la sintaxis de arreglos "clasica"  vector[d] x1[N1]  y multi_normal_cholesky,
//  compatible con la configuracion con la que ya trabajabas. Si tu Stan es >= 2.33 y
//  aparece un error de sintaxis en las declaraciones de arreglos, cambia SOLO:
//      vector[d] x1[N1]   ->   array[N1] vector[d] x1
//      vector[d] x2[N2]   ->   array[N2] vector[d] x2
//      vector[] x   (en las firmas de las funciones)  ->   array[] vector x
//  El resto no cambia (ya no se usa cov_exp_quad, tambien retirado en 2.33).
// =====================================================================

functions {

  // Kernel cuadratico-exponencial de la tesis, SIN jitter ni ruido
  // (el 0.01 y phi=sigma^2 los anade a la diagonal quien llama).
  // rho mayor  =>  decaimiento mas rapido.
  matrix gp_kernel(vector[] x, real eta, real rho) {
    int N = size(x);
    matrix[N, N] K;
    for (i in 1:N) {
      K[i, i] = eta;                       // exp(0) = 1
      for (j in (i + 1):N) {
        real kij = eta * exp(-rho * squared_distance(x[i], x[j]));
        K[i, j] = kij;
        K[j, i] = kij;
      }
    }
    return K;
  }

  // Covarianza cruzada entre dos conjuntos de puntos (solo kernel).
  matrix gp_kernel_cross(vector[] xa, vector[] xb, real eta, real rho) {
    int Na = size(xa);
    int Nb = size(xb);
    matrix[Na, Nb] K;
    for (i in 1:Na)
      for (j in 1:Nb)
        K[i, j] = eta * exp(-rho * squared_distance(xa[i], xb[j]));
    return K;
  }

  // Prediccion posterior del GP (Rasmussen & Williams, Alg. 2.1) con MEDIA LINEAL.
  // Devuelve  m2 + (prediccion GP de los residuos r1 = y1 - media_lineal_train).
  //   jitter = 0.01 (parte de la covarianza del GP) ;  delta = estabilizacion extra.
  vector gp_pred_rng(vector[] x2, vector r1, vector[] x1, vector m2,
                     real eta, real rho, real sigma,
                     real jitter, real delta) {
    int N1 = rows(r1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K;
      matrix[N1, N1] L_K;
      vector[N1] K_div_r1;
      matrix[N1, N2] k_x1_x2;
      matrix[N1, N2] v_pred;
      vector[N2] f2_mu;
      matrix[N2, N2] cov_f2;

      // Cov(residuos) = K(x1,x1) + 0.01 I + sigma^2 I
      K = gp_kernel(x1, eta, rho);
      for (n in 1:N1)
        K[n, n] = K[n, n] + jitter + square(sigma);
      L_K = cholesky_decompose(K);

      K_div_r1 = mdivide_left_tri_low(L_K, r1);
      K_div_r1 = mdivide_right_tri_low(K_div_r1', L_K)';

      // Covarianza cruzada test-train: kernel puro (el jitter es diagonal).
      k_x1_x2 = gp_kernel_cross(x1, x2, eta, rho);

      // Media predictiva = media lineal en test + correccion GP de los residuos
      f2_mu = m2 + (k_x1_x2' * K_div_r1);

      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);

      // Prior latente en test: K(x2,x2) + 0.01 I
      cov_f2 = gp_kernel(x2, eta, rho) - (v_pred' * v_pred);
      for (n in 1:N2)
        cov_f2[n, n] = cov_f2[n, n] + jitter;

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_matrix(rep_vector(delta, N2)));
    }
    return f2;
  }
}

data {
  int<lower=1> N1;            // n de entrenamiento
  int<lower=1> d;             // dimension del vector de regresores del kernel
  vector[d] x1[N1];           // diseno del kernel (train); continuas estandarizadas
  vector[N1] y1;              // outcome ESTANDARIZADO (media 0, sd 1)

  int<lower=1> N2;            // n de prueba (= 2*N1: bloques A=1 y A=0)
  vector[d] x2[N2];           // diseno del kernel (contrafactual)

  real<lower=0> gpscale;      // mediana de las distancias^2 del diseno (calibra rho)

  int<lower=1> q;             // # columnas de la media lineal [1, A, X]
  matrix[N1, q] W1;           // diseno de la media lineal (train)
  matrix[N2, q] W2;           // diseno de la media lineal (test)
  int<lower=1, upper=q> a_col;// indice de la columna de A (para el prior N(0,1) de tau)
}

transformed data {
  real jitter = 0.01;         // jitter de la teoria (parte de la covarianza del GP)
  real delta  = 1e-6;         // estabilizacion numerica extra para multi_normal_rng

  vector[q] prior_sd;         // desviaciones del prior de los coeficientes de la media
  for (j in 1:q)
    prior_sd[j] = 5;          // difuso para intercepto y betas (Y estandarizado)
  prior_sd[a_col] = 1;        // tau ~ N(0,1) : efecto de tratamiento, como en la teoria
}

parameters {
  real<lower=0> rho;          // decaimiento  (convencion de la tesis)
  real<lower=0> alpha;        // amplitud en escala SD ;  eta = alpha^2  (varianza marginal)
  real<lower=0> sigma;        // ruido en escala SD    ;  phi = sigma^2  (varianza del ruido)
  vector[q] gamma;            // coeficientes de la media lineal ; gamma[a_col] = tau
}

model {
  matrix[N1, N1] K;
  matrix[N1, N1] L_K;
  real eta = square(alpha);

  // ---- Priors ----
  // rho ~ Gamma(2, 2*gpscale): media = 1/gpscale, de modo que a la distancia
  // tipica (mediana) la correlacion es exp(-1) ~ 0.37. Solo depende de la geometria
  // de los datos => se reescala solo para el conjunto recortado (PS en [0.1, 0.9]).
  rho   ~ gamma(2, 2 * gpscale);
  alpha ~ normal(0, 1);       // half-normal por <lower=0> ; Y estandarizado
  sigma ~ normal(0, 1);       // half-normal por <lower=0> ; Y estandarizado
  gamma ~ normal(0, prior_sd);// tau ~ N(0,1) ; intercepto y betas difusos

  // Verosimilitud marginal del GP con MEDIA LINEAL (la funcion latente f se integra):
  //   y1 ~ N( W1*gamma ,  K(x1,x1) + 0.01 I + sigma^2 I )
  K = gp_kernel(x1, eta, rho);
  for (n in 1:N1)
    K[n, n] = K[n, n] + jitter + square(sigma);
  L_K = cholesky_decompose(K);

  y1 ~ multi_normal_cholesky(W1 * gamma, L_K);
}

generated quantities {
  vector[N2] f_test;          // mu(A,X) = media lineal + desviacion GP, en test (para el ATE)
  vector[N2] y_test;          // (opcional) realizacion del outcome con ruido

  {
    vector[N1] r1 = y1 - W1 * gamma;   // residuos respecto a la media lineal
    vector[N2] m2 = W2 * gamma;        // media lineal en los puntos de prueba
    f_test = gp_pred_rng(x2, r1, x1, m2, square(alpha), rho, sigma, jitter, delta);
  }

  for (n in 1:N2)
    y_test[n] = normal_rng(f_test[n], sigma);
}
