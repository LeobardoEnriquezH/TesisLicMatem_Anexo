############################################################################
###----------------------------------------------------------------------###
### Autor: Leobardo Enriquez, basado en Oganisian (2020)
### Inferencia causal bayesiana no parametrica del ATE
###   (DPM + BART + GP(Stan) + Regresion lineal bayesiana de referencia)
###
### Identificacion por estandarizacion + Bayesian bootstrap (mismo esquema
### de inferencia para los cuatro modelos), fiel a la teoria de la tesis.
###
### RESUMEN DE CORRECCIONES RESPECTO A LA VERSION PREVIA
### ---------------------------------------------------------------------
###  (1) Tipos de variables. state (codigos INEGI nominales) y pib_class
###      (quintil) entraban como columnas NUMERICAS: el DPM les ponia un
###      local gaussiano y BART partia sobre umbrales de codigos sin sentido.
###      Ahora se codifican como dummies 0/1 (k-1, categoria de referencia =
###      la mas frecuente), de forma UNIFORME en DPM, BART, GP y LM. Asi el
###      DPM usa locales Bernoulli (respeta el soporte, como pide la teoria)
###      y BART/GP/LM ven indicadores y no codigos.
###  (2) Escalas. Se estandarizan las covariables CONTINUAS y el outcome Y con
###      estadisticas de los datos cargados; binarias y dummies se dejan 0/1.
###      Esto vuelve coherentes los priors escritos en la teoria (theta~N(0,1),
###      amplitud del GP ~N(0,1)) y la geometria del kernel del GP. Los ATE se
###      reexpresan en puntos porcentuales multiplicando por sd(Y).
###  (3) GP/Stan. El kernel ahora se evalua tal cual la teoria
###      C_ij = eta*exp(-rho*||.||^2) + 0.01, con priors calibrados a la escala
###      de distancias del diseno (se pasa gpscale a Stan). Se corrige la
###      contradiccion de parametrizacion de cov_exp_quad (ver el .stan).
###  (4) Regresion lineal de referencia. La teoria pide una regresion lineal
###      BAYESIANA (theta~N(0,1), beta~N(mu0,Sigma0), phi~IG(a0,b0)) cuyo theta
###      es el ATE. Se implementa con un muestreador de Gibbs conjugado
###      (Normal-Inverse-Gamma) y se alimenta al MISMO Bayesian bootstrap, en
###      lugar del bootstrap frecuentista de OLS previo.
###  (5) Robustez. Nada se fija a n=281: todo se adapta a nrow(datos), se
###      eliminan columnas de varianza cero y las dummies se construyen con los
###      niveles PRESENTES. El mismo script corre tal cual sobre el conjunto
###      recortado (propensity score en [0.1, 0.9]).
###  (6) Resultados. Las extracciones psi se guardan en una LISTA con nombre y
###      todos los resumenes se calculan por vector (se elimina el cbind de
###      vectores de longitudes distintas que reciclaba en silencio).
###----------------------------------------------------------------------###
############################################################################

## ---------------------------------------------------------------------
## Paquetes
## ---------------------------------------------------------------------
## ChiRP: la rama "fast-reg" ya no existe (devuelve 404). Se instala desde la
## rama por defecto del repositorio.
if (!requireNamespace("ChiRP", quietly = TRUE)) {
  devtools::install_github("stablemarkets/ChiRP")
}
library(ChiRP)
library(BayesTree)
library(rstan)
library(gtools)      # rdirichlet (Bayesian bootstrap)
library(coda)        # HPDinterval
library(latex2exp)   # TeX en la grafica

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## ---------------------------------------------------------------------
## Ruta de datos
## ---------------------------------------------------------------------
## Para el conjunto RECORTADO (propensity score en [0.1, 0.9]) basta cambiar
## esta ruta por la del archivo recortado; todo el pipeline se adapta a la n.
DATA_PATH <- "dataset_trimmed.xlsx"

## ---------------------------------------------------------------------
## Bayesian bootstrap  (Rubin, 1981): estandarizacion del contraste causal
## ---------------------------------------------------------------------
## mu_a1, mu_a0: matrices  n x M  (filas = unidades, columnas = draws MCMC).
## Devuelve M realizaciones posteriores de Psi.
bayes_boot <- function(mu_a1, mu_a0) {
  n <- nrow(mu_a1)
  M <- ncol(mu_a1)
  psi_post <- numeric(M)
  for (m in 1:M) {
    bb_weights  <- as.numeric(rdirichlet(1, rep(1, n)))
    psi_post[m] <- sum(bb_weights * (mu_a1[, m] - mu_a0[, m]))
  }
  psi_post
}

## ---------------------------------------------------------------------
## Dummies de referencia (k-1) a partir de los niveles PRESENTES
## ---------------------------------------------------------------------
## Referencia = categoria mas frecuente (robusta al recorte). Devuelve NULL si
## la variable no tiene variacion en los datos cargados.
make_dummies <- function(x, prefix) {
  x  <- factor(x)
  lv <- levels(x)
  if (length(lv) <= 1) return(NULL)
  ref    <- names(which.max(table(x)))   # categoria de referencia
  lv_use <- setdiff(lv, ref)
  M <- vapply(lv_use, function(l) as.numeric(x == l), numeric(length(x)))
  if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  colnames(M) <- paste0(prefix, "_", lv_use)
  M
}

#### -------------------------- Datos ----------------------------------- ####
set.seed(1)
dat <- readxl::read_excel(DATA_PATH)

## Mapeo a una sola tabla con nombres cortos
d <- data.frame(
  Y = dat$outcome,            # cambio en informalidad 2023:2024 - 2017:2018 (puntos %)
  A = dat$trt,                # tratamiento (1 = tratado)            [binaria]
  V = dat$state,              # estado, codigos INEGI                [NOMINAL -> dummies]
  I = dat$informal,           # informalidad base                   [continua]
  M = dat$maquila,            # maquila importante (1)               [binaria]
  S = dat$imssissste,         # acceso IMSS/ISSSTE (%)               [continua]
  G = dat$pib_class,          # quintil de PIB 2018                  [ordinal -> dummies]
  U = dat$urban,              # municipio urbano (1)                 [binaria]
  D = dat$lpibpc,             # log(PIB per capita) 2018             [continua]
  P = dat$poverty,            # pobreza (%)                          [continua]
  E = dat$share_manufacture,  # participacion manufactura (%)        [continua]
  O = dat$lpopulation         # log(poblacion)                       [continua]
)
N <- nrow(d)

## (Opcional) diagnosticos rapidos como en la version original
plot(d$I, d$Y, col = ifelse(d$A == 1, "red", "black"),
     xlab = "Informalidad base", ylab = "Outcome")

## ---------------------------------------------------------------------
## Diseno UNIFICADO (mismas columnas para DPM, BART, GP y LM)
## ---------------------------------------------------------------------
cont_vars <- c("I", "S", "D", "P", "E", "O")   # continuas (se estandarizan)
bin_vars  <- c("M", "U")                        # binarias (se dejan 0/1)

## Estandarizacion de continuas con estadisticas de los datos cargados
Xcont   <- as.matrix(d[, cont_vars, drop = FALSE])
x_mean  <- colMeans(Xcont)
x_sd    <- apply(Xcont, 2, sd)
Xcont_s <- scale(Xcont, center = x_mean, scale = x_sd)

## Binarias 0/1
Xbin <- as.matrix(d[, bin_vars, drop = FALSE])

## Dummies de state (V) y pib_class (G)
Vd <- make_dummies(d$V, "V")
Gd <- make_dummies(d$G, "G")
Xdum <- do.call(cbind, Filter(Negate(is.null), list(Vd, Gd)))

## Matriz de covariables (sin A y sin intercepto)
Xcov_mat <- cbind(Xcont_s, Xbin, Xdum)

## Eliminar columnas de varianza cero (robustez ante el recorte)
sds  <- apply(Xcov_mat, 2, sd)
keep <- is.finite(sds) & sds > 0
Xcov_mat <- Xcov_mat[, keep, drop = FALSE]
Xcov_df  <- as.data.frame(Xcov_mat, check.names = FALSE)

## Estandarizacion del outcome (los cuatro modelos se ajustan en esta escala;
## el ATE se reexpresa en puntos porcentuales multiplicando por y_sd)
y_mean <- mean(d$Y)
y_sd   <- sd(d$Y)
y_s    <- as.vector((d$Y - y_mean) / y_sd)

## Disenos de entrenamiento y contrafactuales (A = 1 y A = 0 sobre las mismas X)
A  <- d$A
design_train <- cbind(A = A, Xcov_mat)                 # GP / BART (sin intercepto)
A1 <- cbind(A = 1, Xcov_mat)
A0 <- cbind(A = 0, Xcov_mat)
design_test  <- rbind(A1, A0)                          # 2N x p (bloque A=1, luego A=0)

#### ----------------------- Dirichlet Process -------------------------- ####
## fDPMix detecta el tipo por columna: 2 valores unicos -> local Bernoulli;
## >2 -> local gaussiana. Con el diseno unificado, A y las dummies son 0/1
## (Bernoulli) y las continuas estandarizadas son gaussianas, respetando el
## soporte como exige la teoria. NO se agrega intercepto (fDPMix lo incluye
## internamente; las dummies k-1 evitan colinealidad con ese intercepto).
d_train_dp <- data.frame(Y = y_s, A = A, Xcov_df, check.names = FALSE)
d_test_dp  <- rbind(
  data.frame(A = 1, Xcov_df, check.names = FALSE),
  data.frame(A = 0, Xcov_df, check.names = FALSE)
)
dp_formula <- reformulate(c("A", colnames(Xcov_mat)), response = "Y")

## NOTA SOBRE EL ENCOGIMIENTO DEL DPM HACIA CERO
## El DPM modela la DENSIDAD CONJUNTA (Y, A, X). Por defecto, los coeficientes de
## la regresion local de cada cluster tienen prior centrado en el OLS global y
## dispersion beta_var_scale (default 1000 ~ casi sin agrupamiento parcial). Bajo
## fuerte desbalance y tratamiento escaso (43/281), los clusters segregan tratados
## y controles: dentro de cada cluster casi no hay variacion en A, el efecto local
## queda mal identificado y la diferencia contrafactual se "lava" hacia cero.
##   Palanca: beta_var_scale mas pequeno => agrupamiento parcial (los efectos por
##   cluster se acercan al efecto OLS global), lo que reduce el lavado hacia cero a
##   costa de algo de flexibilidad. Aqui se usa un valor moderado; subelo para mas
##   heterogeneidad o bajalo para mas pooling. El remedio PRINCIPAL del encogimiento
##   del DPM sigue siendo el solape: corre este mismo script sobre el conjunto
##   recortado (PS en [0.1, 0.9]), donde cada cluster contiene ambos brazos.
set.seed(2)
res <- fDPMix(
  d_train        = d_train_dp,
  formula        = dp_formula,
  d_test         = d_test_dp,
  iter           = 2000,
  burnin         = 1000,
  init_k         = 10,
  mu_scale       = 2,      # dispersion de las medias de cluster (geometria del agrupamiento)
  beta_var_scale = 25      # agrupamiento parcial de los efectos por cluster (default 1000)
)

## res$test: filas 1..N -> bloque A=1 ; filas (N+1)..2N -> bloque A=0
psi_dp_std <- bayes_boot(res$test[1:N, , drop = FALSE],
                         res$test[(N + 1):(2 * N), , drop = FALSE])

#### ----------------------- BART ---------------------------------------- ####
## Diseno numerico con dummies alineadas (train y test comparten construccion).
## BART reescala y internamente; pasamos y estandarizado por uniformidad.
set.seed(3)
bart_res <- bart(
  x.train = design_train,
  y.train = y_s,
  x.test  = design_test,
  ndpost  = 500,
  nskip   = 500,
  verbose = FALSE
)

bart_pred  <- t(bart_res$yhat.test)                    # 2N x draws
psi_bart_std <- bayes_boot(bart_pred[1:N, , drop = FALSE],
                           bart_pred[(N + 1):(2 * N), , drop = FALSE])

#### ----------------------- Gaussian Process (Stan) -------------------- ####
## Diseno del kernel = A + continuas estandarizadas + binarias/dummies 0/1 (sin
## intercepto, vacuo en un kernel RBF). Las dummies NO se estandarizan: estandarizar
## indicadores poco frecuentes infla sus distancias y distorsiona el kernel isotropico.
##
## MEDIA LINEAL (anti-encogimiento). En la version previa el GP tenia media 0: donde
## los datos escasean o el kernel acopla los brazos, revierte a 0 y el ATE se atenua.
## Aqui el GP modela DESVIACIONES respecto a una media lineal m(A,X) = gamma0 + tau*A +
## X'beta (con tau ~ N(0,1), como el efecto de la teoria). Asi la reversion ocurre
## hacia la regresion lineal (efecto no nulo), no hacia cero, y el kernel sigue siendo
## EXACTAMENTE el de la tesis eta*exp(-rho*||.||^2)+0.01. El ATE = parte lineal (tau) +
## desviacion flexible del GP; mantiene A tambien en el kernel para captar
## heterogeneidad/no linealidad del tratamiento sin colapsar a cero.
## (Alternativa aun mas anclada: quitar A del kernel para que TODO el efecto pase por
##  tau; quedaria practicamente igual al modelo lineal.)
##
## gpscale = mediana de las distancias^2 del diseno (misma escala que x1): calibra el
## prior de rho a la geometria de ESTOS datos (full o recortados).
d2      <- as.matrix(dist(design_train))^2
gpscale <- median(d2[upper.tri(d2)])

## Diseno de la media lineal: [intercepto, A, X]  (train y test)
WG_train <- cbind(Intercept = 1, design_train)
WG_test  <- cbind(Intercept = 1, design_test)
a_col_gp <- which(colnames(WG_train) == "A")     # indice de A (para el prior N(0,1) de tau)

## Convertir matriz a arreglo vector[d] de Stan
mat_to_vec_array <- function(X) lapply(seq_len(nrow(X)), function(i) as.numeric(X[i, ]))

stan_data <- list(
  N1      = nrow(design_train),
  d       = ncol(design_train),
  x1      = mat_to_vec_array(design_train),
  y1      = y_s,
  N2      = nrow(design_test),
  x2      = mat_to_vec_array(design_test),
  gpscale = gpscale,
  q       = ncol(WG_train),
  W1      = WG_train,
  W2      = WG_test,
  a_col   = a_col_gp
)

stan_mod <- stan_model("gaussian_process.stan")

set.seed(1000)
stan_res <- sampling(
  stan_mod, data = stan_data,
  chains = 4, warmup = 1000, iter = 2000, seed = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

## f_test: media de respuesta mu(X) = media lineal + desviacion GP, en los contrafactuales
f_test  <- rstan::extract(stan_res, pars = "f_test")$f_test   # draws x 2N
gp_pred <- t(f_test)                                          # 2N x draws

psi_gp_std <- bayes_boot(gp_pred[1:N, , drop = FALSE],
                         gp_pred[(N + 1):(2 * N), , drop = FALSE])

#### ----------------------- Regresion lineal bayesiana ----------------- ####
## Teoria (subseccion de referencia):
##   Y_i ~ N(theta*A_i + X_i' beta, phi),  theta ~ N(0,1),
##   beta ~ N(mu0, Sigma0),  phi ~ IG(a0, b0).
## Bajo linealidad/aditividad, theta es el ATE; aqui ademas se post-procesa con
## el MISMO Bayesian bootstrap para mantener identico el esquema de inferencia.
## Muestreador de Gibbs conjugado (Normal-Inverse-Gamma), sin dependencias extra.
set.seed(123)

W  <- cbind(Intercept = 1, design_train)               # [1, A, X]   (N x p)
W1 <- cbind(Intercept = 1, A1)                          # contrafactual A = 1
W0 <- cbind(Intercept = 1, A0)                          # contrafactual A = 0
p  <- ncol(W)
a_idx <- which(colnames(W) == "A")                      # columna de theta

## Prior:  theta ~ N(0,1) ; intercepto y betas ~ N(0, 100) (difuso en escala estand.)
mu0   <- rep(0, p)
s0    <- rep(100, p); s0[a_idx] <- 1
S0inv <- diag(1 / s0, p)
a0 <- 0.01; b0 <- 0.01                                   # IG poco informativo

niter <- 2000; burn <- 1000; keepn <- niter - burn
XtX <- crossprod(W)                                      # W'W
Xty <- crossprod(W, y_s)                                 # W'y

mu_a1_lm <- matrix(0, N, keepn)
mu_a0_lm <- matrix(0, N, keepn)
phi <- var(y_s); s <- 0
for (it in 1:niter) {
  Vpost <- solve(XtX / phi + S0inv)
  mpost <- Vpost %*% (Xty / phi + S0inv %*% mu0)
  gamma <- as.vector(mpost + t(chol(Vpost)) %*% rnorm(p))   # gamma | phi, Y
  resid <- as.vector(y_s - W %*% gamma)
  phi   <- 1 / rgamma(1, shape = a0 + N / 2, rate = b0 + 0.5 * sum(resid^2))  # phi | gamma, Y
  if (it > burn) {
    s <- s + 1
    mu_a1_lm[, s] <- as.vector(W1 %*% gamma)
    mu_a0_lm[, s] <- as.vector(W0 %*% gamma)
  }
}
psi_lm_std <- bayes_boot(mu_a1_lm, mu_a0_lm)

#### ----------------------- Resultados --------------------------------- ####
## Reexpresion a puntos porcentuales (la media de Y se cancela en la diferencia,
## por lo que solo se multiplica por y_sd). Se guardan en una LISTA con nombre;
## todos los resumenes se calculan POR VECTOR (sin cbind de longitudes distintas).
psi_list <- list(
  "BART"              = psi_bart_std * y_sd,
  "Dirichlet Process" = psi_dp_std   * y_sd,
  "Gaussian Process"  = psi_gp_std   * y_sd
  #"Linear Model"      = psi_lm_std   * y_sd
)

post_mean <- sapply(psi_list, mean)
p_pos     <- sapply(psi_list, function(v) mean(v > 0))   # P(ATE > 0)
p_neg     <- sapply(psi_list, function(v) mean(v < 0))   # P(ATE < 0)  (teoria)
hpdi90    <- sapply(psi_list, function(v) {
  h <- HPDinterval(as.mcmc(v), prob = 0.90)
  c(lo = h[1], hi = h[2])
})

tab <- data.frame(
  Model     = names(psi_list),
  n         = N,
  E_ATE     = round(post_mean, 3),
  P_ATE_pos = round(p_pos, 3),
  P_ATE_neg = round(p_neg, 3),
  HPDI90_lo = round(hpdi90["lo", ], 3),
  HPDI90_hi = round(hpdi90["hi", ], 3),
  row.names = NULL
)
print(tab)

#### ----------------------- Grafica ------------------------------------ ####
mods <- names(psi_list)
K    <- length(mods)
rng  <- range(unlist(psi_list))
## Espacio extra ARRIBA para alojar la tabla sin encimarse con las bandas
ylim <- c(rng[1] - 0.10 * diff(rng), rng[2] + 0.45 * diff(rng))

op <- par(no.readonly = TRUE)
par(mar = c(5, 4, 2, 1))

plot(seq_len(K), post_mean, pch = 20, col = "blue",
     xlim = c(0.5, K + 0.5), ylim = ylim,
     ylab = TeX("$ATE$ (puntos porcentuales)"),
     xlab = "Modelo", axes = FALSE)
axis(1, at = seq_len(K), labels = mods, padj = 0.5)

## Eje vertical con mas valores: etiquetas en cada entero y marcas menores cada 0.5
y_lo <- floor(ylim[1]); y_hi <- ceiling(ylim[2])
axis(2, at = seq(y_lo, y_hi, by = 1), las = 1)                      # etiquetas (enteros)
axis(2, at = seq(y_lo, y_hi, by = 0.5), labels = FALSE, tcl = -0.25) # marcas menores

## Bandas creibles por modelo (calculadas de cada vector)
#colfunc <- colorRampPalette(c("white", "lightblue"))
ci_perc <- c(.90, .85, .80)
#offsets <- c(-0.04, 0, 0.04)
offsets <- c(0, 0, 0)

colvec <- c("lightblue1", "lightblue3", "lightgray")
lwd_vec <- c(8, 8, 8)

for (i in seq_along(ci_perc)) {
  for (j in seq_len(K)) {
    q <- quantile(
      psi_list[[j]],
      probs = c((1 - ci_perc[i]) / 2, (1 + ci_perc[i]) / 2),
      na.rm = TRUE
    )
    
    segments(
      x0 = j + offsets[i], y0 = q[1],
      x1 = j + offsets[i], y1 = q[2],
      col = colvec[i],
      lwd = lwd_vec[i],
      lend = "round"
    )
  }
}

points(seq_len(K), post_mean, col = "steelblue", pch = 17, cex = 1.2)
abline(h = 0, col = "red", lty = 2)
legend("topleft",
       legend = c("Media posterior + Intervalos creibles centrales", "Sin efecto"),
       col = c("steelblue", "red"), pch = c(17, NA), lty = c(NA, 2), bty = "n")

## ---- Tabla embebida en la esquina superior derecha ----
tab_lines <- c(
  sprintf("%-17s %6s %5s %5s %6s %6s",
          "Modelo", "ATE", "P>0", "P<0", "L90", "U90"),
  apply(tab, 1, function(r)
    sprintf("%-17s %6.2f %5.2f %5.2f %6.2f %6.2f",
            r[["Model"]], as.numeric(r[["E_ATE"]]), as.numeric(r[["P_ATE_pos"]]),
            as.numeric(r[["P_ATE_neg"]]), as.numeric(r[["HPDI90_lo"]]),
            as.numeric(r[["HPDI90_hi"]])))
)
old_family <- par("family"); par(family = "mono")     # fuente monoespaciada -> columnas alineadas
legend("topright",
       legend = tab_lines,
       title  = "Resumen del ATE (puntos porcentuales; HPDI 90%)",
       bty = "o", bg = rgb(1, 1, 1, 0.9), box.col = "grey70",
       cex = 0.72, adj = 0)
par(family = old_family)

par(op)
