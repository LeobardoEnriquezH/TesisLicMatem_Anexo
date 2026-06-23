data {
  int<lower=1> N;                 // number of observations
  int<lower=1> Pv;                // number of strata (dimension of V)
  int<lower=1> Pw;                // number of confounders (dimension of W)
  matrix[N, Pw] W;                // confounders (include intercept column if you want)
  vector[N] A;                    // treatment indicator (0/1)
  vector[N] Y;                    // continuous outcome
  array[N] int<lower=1, upper=Pv> state; 
  array[Pv] int<lower=0> n_v;         
  array[Pv+1] int<lower=0> ind;    // start indices (0-based), so v_start = ind[v]+1
  int<lower=1> Kt;
  matrix[N, Kt] Wt;
}

parameters {
  vector[Pw] beta_w;          // incluye intercepto global en W
  vector[Pv] z_theta;         // efectos estandarizados
  real mu;                    // media global del efecto tratamiento
  real alpha0;
  real<lower=0> tau_alpha;
  vector[Pv] z_alpha;
  real<lower=0> tau;          // dispersión entre estados
  real<lower=0> sigma_y;      // escala outcome
  real<lower=0> nu_minus_two;
  vector[Kt] gamma;
}

transformed parameters {
  vector[Pv] theta;
  vector[Pv] alpha_state;
  real nu;

  theta = mu + tau * z_theta;
  alpha_state = alpha0 + tau_alpha * z_alpha;
  nu = nu_minus_two + 2;
}


model {
  beta_w ~ normal(0, 2.5);
  
  z_theta ~ normal(0, 1);
  mu ~ normal(0, 2.5);
  tau ~ normal(0, 2);      // half-normal por lower=0
  
  alpha0 ~ normal(0, 2.5);
  tau_alpha ~ normal(0, 2);
  z_alpha ~ normal(0, 1);
  
  sigma_y ~ normal(0, 2);  // half-normal por lower=0
  nu_minus_two ~ exponential(0.1);
  gamma ~ normal(0, 1);
  
  for (i in 1:N) {
    real eta;
    eta = W[i] * beta_w
    + alpha_state[state[i]]
    + (theta[state[i]] + Wt[i] * gamma) * A[i];

    Y[i] ~ student_t(nu, eta, sigma_y);
  }
}


generated quantities {
  vector[Pv] ate;
  real overall_ate;
  vector[N] y_rep;

  // CATE por estado con Bayesian Bootstrap
  for (v in 1:Pv) {
    int nv = n_v[v];
    int v_start = ind[v] + 1;
    int v_end = ind[v + 1];
    vector[nv] diff_v;
    vector[nv] w_bb;

    for (j in 1:nv) {
      int i = v_start + j - 1;

      // Diferencia individual bajo tratamiento vs control
      diff_v[j] = theta[v] + Wt[i] * gamma;
    }

    w_bb = dirichlet_rng(rep_vector(1.0, nv));
    ate[v] = dot_product(w_bb, diff_v);
  }

  // Efecto promedio global con Bayesian Bootstrap
  {
    vector[N] diff_all;
    vector[N] w_all;

    for (i in 1:N) {
      diff_all[i] = theta[state[i]] + Wt[i] * gamma;
    }

    w_all = dirichlet_rng(rep_vector(1.0, N));
    overall_ate = dot_product(w_all, diff_all);
  }

  // Réplicas posteriores
  for (i in 1:N) {
    real eta;
    eta = W[i] * beta_w
        + alpha_state[state[i]]
        + (theta[state[i]] + Wt[i] * gamma) * A[i];

    y_rep[i] = student_t_rng(nu, eta, sigma_y);
  }
}

