###--------------------------------------------------------------------------###
### Author: Leobardo Enriquez, based on Oganisian (2020)
### Conduct partial pooling analysis
###--------------------------------------------------------------------------###

## Load Packages
library(rstan)
library(LaplacesDemon)
library(latex2exp)
library(coda)
set.seed(66)


####------------------------ ---------Data ---------------------------------####
warmup = 2000
iter = 4000
n_draws = iter - warmup

#import confounder (I,M,S,U, etc), treatment (A), and outcome (Y)
dat<-readxl::read_xlsx("dataset.xlsx")
N = nrow(dat) # sample size

W_df <- data.frame(
  I = dat$informal,
  M = dat$maquila,
  S = dat$imssissste,
  G = dat$pib_class,
  U = dat$urban,
  D = dat$lpibpc,
  P = dat$poverty,
  E = dat$share_manufacture,
  O = dat$lpopulation
)

# estandarizar covariables continuas
W_df$I <- as.numeric(scale(W_df$I))
W_df$S <- as.numeric(scale(W_df$S))
W_df$D <- as.numeric(scale(W_df$D))
W_df$P <- as.numeric(scale(W_df$P))
W_df$E <- as.numeric(scale(W_df$E))
W_df$O <- as.numeric(scale(W_df$O))

Wmat <- model.matrix(~ 0 + I + M + S + G + U + D + P + E + O, data = W_df)

# Covariables que modificarán el efecto del tratamiento
Wtrt <- as.matrix(W_df[, c("I", "S", "U")])
Kt <- ncol(Wtrt)


# --------------------------------------------------------------------------
# Variable categórica de estratificación
# --------------------------------------------------------------------------
# Cambiar la clasificación que se quiere graficar:
# "state", "pib_class", "informal_quintile", etc.
stratum_var  <- "state"
stratum_xlab <- "Estado"

#stratum_var  <- "pib_class"
#stratum_xlab <- "Quintil de PIB"

#stratum_var  <- "informal_class"
#stratum_xlab <- "Quintil de informalidad inicial"

# Función para conservar orden lógico:
# - si es factor, respeta el orden de sus niveles;
# - si es numérica, ordena de menor a mayor;
# - si es texto, ordena alfabéticamente.
get_stratum_levels <- function(x) {
  if (is.factor(x)) {
    levels(droplevels(x))
  } else if (is.numeric(x) || is.integer(x)) {
    sort(unique(x[!is.na(x)]))
  } else {
    sort(unique(as.character(x[!is.na(x)])))
  }
}

v_raw    <- dat[[stratum_var]]
v_levels <- get_stratum_levels(v_raw)
v_fac    <- factor(v_raw, levels = v_levels)
v_labels <- as.character(v_levels)

# Identificador entero que Stan necesita.
# Se conserva el nombre 'state' porque el archivo .stan
# usa una variable llamada state.
state <- as.integer(v_fac)
Pv <- nlevels(v_fac)

if (anyNA(state)) {
  stop("Hay observaciones con NA en la variable de estratificación: ", stratum_var)
}

A=dat$trt #Treatment indicator

Y=as.numeric(dat$outcome)#Outcome:change in informal employment from 2017–2018 to 2025



ord <- order(state)

state_ord <- state[ord]
A_ord <- dat$trt[ord]
Y_ord <- as.numeric(dat$outcome[ord])
W_ord <- Wmat[ord, , drop = FALSE]
Wtrt_ord <- Wtrt[ord, , drop = FALSE]

# Tamaño de muestra por estrato, respetando el orden de v_labels
n_v <- as.numeric(table(factor(state_ord, levels = seq_len(Pv))))

stan_data <- list(
  N = N,
  Pv = Pv,
  Pw = ncol(W_ord),
  Kt = Kt,
  state = state_ord,
  W = W_ord,
  Wt = Wtrt_ord,
  A = A_ord,
  Y = Y_ord,
  n_v = n_v,
  ind = c(0, cumsum(n_v))
)

####------------------------ Sample Posterior    ---------------------------####
partial_pool_model = stan_model(file = "partial_pool_BB.stan")


stan_res <- sampling(
  partial_pool_model,
  data = stan_data,
  pars = c("ate", "overall_ate", "theta", "mu", "tau", "beta_w", "gamma",
           "alpha_state", "alpha0", "tau_alpha", "sigma_y", "nu", "y_rep"),
  warmup = warmup,
  iter = iter,
  chains = 4,
  cores = 4,
  seed = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

ATE_draws = extract(stan_res, pars='ate')[[1]]
overall_mean = extract(stan_res, pars='overall_ate')[[1]]

print(stan_res, pars = c("mu", "tau", "alpha0", "tau_alpha", "sigma_y", "nu", "gamma"))
rstan::traceplot(stan_res, pars = c("mu", "tau", "alpha0", "tau_alpha", "sigma_y", "nu", "gamma"))

####------------------- Compute Frequentist Estimates   --------------------####


freq_df <- data.frame(
  Y = as.numeric(dat$outcome),
  A = dat$trt,
  V = v_fac,
  W_df
)

freq_reg <- lm(
  Y ~ A * V + A:(I + S + U) + I + M + S + G + U + D + P + E + O,
  data = freq_df
)

Psi_freq <- numeric(Pv)

for (vf in seq_len(Pv)) {
  lev <- levels(freq_df$V)[vf]
  
  new1 <- freq_df[freq_df$V == lev, , drop = FALSE]
  new0 <- new1
  
  new1$A <- 1
  new0$A <- 0
  
  p1 <- predict(freq_reg, newdata = new1)
  p0 <- predict(freq_reg, newdata = new0)
  
  Psi_freq[vf] <- mean(p1 - p0)
}


############# convergence check

print(stan_res, pars = c("mu", "tau", "alpha0", "tau_alpha", "sigma_y", "nu", "theta"))
summary(stan_res)$summary[, c("Rhat", "n_eff")]




stan_rhat <- summary(stan_res)$summary[, "Rhat"]
stan_neff <- summary(stan_res)$summary[, "n_eff"]

print(range(stan_rhat, na.rm = TRUE))
print(range(stan_neff, na.rm = TRUE))


#######################posterior predictive checking
####################### posterior predictive checking #######################

y_rep_draws <- rstan::extract(stan_res, pars = "y_rep")[[1]]

# media y desviación observadas
y_obs_mean <- mean(Y_ord)
y_obs_sd   <- sd(Y_ord)

# medias y desviaciones replicadas
y_rep_means <- apply(y_rep_draws, 1, mean)
y_rep_sds   <- apply(y_rep_draws, 1, sd)

cat("Posterior predictive p-value for mean: ",
    mean(y_rep_means >= y_obs_mean), "\n")
cat("Posterior predictive p-value for sd: ",
    mean(y_rep_sds >= y_obs_sd), "\n")

# Histograma de medias replicadas
hist(y_rep_means, breaks = 30,
     main = "Posterior predictive check: media de Y",
     xlab = "Media replicada")
abline(v = y_obs_mean, lwd = 2, lty = 2)

# Histograma de desviaciones replicadas
hist(y_rep_sds, breaks = 30,
     main = "Posterior predictive check: sd de Y",
     xlab = "Desviación estándar replicada")
abline(v = y_obs_sd, lwd = 2, lty = 2)


######################### Plot  #########################

######################### Plot  #########################

# Number of strata actually present in the posterior
Pv_post <- ncol(ATE_draws)

# Etiquetas de los estratos efectivamente graficados
v_labs_post <- v_labels[seq_len(Pv_post)]

# Posterior mean by stratum
post_mean <- colMeans(ATE_draws)

# Posterior probabilities by stratum
p_ate_gt0 <- colMeans(ATE_draws > 0)
p_ate_lt0 <- colMeans(ATE_draws < 0)

# Posterior quantiles
q85 <- apply(ATE_draws, 2, quantile, probs = 0.85)
q90 <- apply(ATE_draws, 2, quantile, probs = 0.90)

# ---------------- HPDI (coda) ----------------

hpdi_levels <- c(.90, .85, .80)

get_hpdi <- function(draws_mat, level) {
  out <- sapply(seq_len(ncol(draws_mat)), function(j) {
    h <- coda::HPDinterval(coda::as.mcmc(draws_mat[, j]), prob = level)
    c(lower = h[1], upper = h[2])
  })
  out
}

hpdi_list <- lapply(hpdi_levels, function(lv) get_hpdi(ATE_draws, lv))
names(hpdi_list) <- as.character(hpdi_levels)

# ---------------- Summary table ----------------

# Sample sizes by stratum, respetando el orden de v_labels
n_byV <- as.numeric(table(factor(state, levels = seq_len(Pv))))
n_byV <- n_byV[seq_len(Pv_post)]

tab <- data.frame(
  v = v_labs_post,
  n = n_byV,
  post_mean = post_mean,
  p_ate_gt0 = p_ate_gt0,
  p_ate_lt0 = p_ate_lt0,
  q85 = q85,
  q90 = q90,
  freq_MLE = Psi_freq[seq_len(Pv_post)],
  check.names = FALSE
)

# HPDI principal para reportar en el cuadro
lv0 <- as.character(hpdi_levels[1])

tab$HPDI_lo <- hpdi_list[[lv0]]["lower", ]
tab$HPDI_hi <- hpdi_list[[lv0]]["upper", ]

names(tab)[names(tab) == "HPDI_lo"] <- paste0("HPDI", lv0, "_lo")
names(tab)[names(tab) == "HPDI_hi"] <- paste0("HPDI", lv0, "_hi")

# Redondeo para el cuadro
tab_fmt <- tab
num_cols <- setdiff(names(tab_fmt), c("v", "n"))
tab_fmt[num_cols] <- lapply(tab_fmt[num_cols], function(x) round(x, 3))

hpdi_lo_name <- paste0("HPDI", lv0, "_lo")
hpdi_hi_name <- paste0("HPDI", lv0, "_hi")

tab_small <- tab_fmt[, c(
  "v",
  "n",
  "post_mean",
  "p_ate_gt0",
  "p_ate_lt0",
  hpdi_lo_name,
  hpdi_hi_name
)]

colnames(tab_small) <- c(
  "v",
  "n",
  "E[ATE]",
  "P(ATE>0)",
  "P(ATE<0)",
  paste0("HPDI", lv0, "_lo"),
  paste0("HPDI", lv0, "_hi")
)

# ---------------- Plot + table in the same graph ----------------
# ---------------- Plot + table in the same graph ----------------

op <- par(no.readonly = TRUE)

# Más margen inferior para nombres largos de entidades o categorías
par(mar = c(8, 4, 2, 1))

# Texto de la tabla
txt <- capture.output(
  print(tab_small, row.names = FALSE, right = TRUE)
)

# ---------------- Rango vertical de los datos ----------------

data_ylim <- range(
  unlist(hpdi_list),
  post_mean,
  Psi_freq[seq_len(Pv_post)],
  0,
  mean(overall_mean),
  na.rm = TRUE
)

data_pad <- 0.06 * diff(data_ylim)
data_ylim <- c(data_ylim[1] - data_pad, data_ylim[2] + data_pad)

# Reservar una franja inferior dentro de la gráfica para la tabla
n_table_lines <- length(txt)

table_gap  <- 0.08 * diff(data_ylim)
table_band <- max(
  0.30 * diff(data_ylim),
  0.055 * diff(data_ylim) * n_table_lines
)

ylim_range <- c(
  data_ylim[1] - table_gap - table_band,
  data_ylim[2]
)

# Rango horizontal: solo la región de categorías
# ya no se reserva espacio a la derecha para la tabla
xlim_range <- c(0.75, Pv_post + 0.25)

plot(
  seq_len(Pv_post),
  post_mean,
  type = "n",
  xlim = xlim_range,
  ylim = ylim_range,
  ylab = TeX("$\\psi(v)$"),
  xlab = "",
  axes = FALSE
)

# ---------------- Eje X con nombres de categorías ----------------

axis_labs <- paste0(v_labs_post, "\n(n=", n_byV, ")")

axis(
  side = 1,
  at = seq_len(Pv_post),
  labels = FALSE
)

usr <- par("usr")

text(
  x = seq_len(Pv_post),
  y = usr[3] - 0.045 * diff(usr[3:4]),
  labels = axis_labs,
  srt = 45,
  adj = 1,
  xpd = TRUE,
  cex = 0.8
)

mtext(stratum_xlab, side = 1, line = 6.5)

yticks <- pretty(ylim_range)

axis(
  side = 2,
  at = yticks,
  labels = yticks
)

# ---------------- HPDI bands ----------------

hpdi_plot_levels <- sort(hpdi_levels, decreasing = TRUE)

colfunc <- colorRampPalette(c("white", "lightblue"))
colvec <- colfunc(length(hpdi_plot_levels))
names(colvec) <- as.character(hpdi_plot_levels)

# Dibujar bandas HPDI de la más amplia a la más estrecha
for (lv in hpdi_plot_levels) {
  h <- hpdi_list[[as.character(lv)]]
  
  segments(
    x0 = seq_len(Pv_post),
    y0 = h["lower", ],
    x1 = seq_len(Pv_post),
    y1 = h["upper", ],
    col = colvec[as.character(lv)],
    lwd = 10,
    lend = "round"
  )
}

# Medias posteriores
points(
  seq_len(Pv_post),
  post_mean,
  col = "steelblue",
  pch = 17,
  cex = 1.2
)

# Estimaciones frecuentistas / MLE
points(
  seq_len(Pv_post),
  Psi_freq[seq_len(Pv_post)],
  col = "black",
  pch = 19,
  cex = 1.0
)

# ---------------- Líneas de referencia ----------------

segments(
  x0 = 0.75,
  y0 = mean(overall_mean),
  x1 = Pv_post + 0.2,
  y1 = mean(overall_mean),
  col = "steelblue",
  lty = 2
)

segments(
  x0 = 0.75,
  y0 = 0,
  x1 = Pv_post + 0.2,
  y1 = 0,
  col = "red",
  lty = 3
)

legend(
  "topleft",
  legend = c(
    paste0(
      "Posterior mean + HPDI bands (",
      paste(hpdi_plot_levels * 100, "%", collapse = ", "),
      ")"
    ),
    "MLE",
    "Overall effect",
    "No effect"
  ),
  col = c("steelblue", "black", "steelblue", "red"),
  pch = c(17, 19, NA, NA),
  lty = c(NA, NA, 2, 3),
  bty = "n",
  cex = 0.9
)

# ---------------- Tabla dentro de la gráfica, abajo ----------------

# Tamaño inicial de la tabla
table_cex <- 0.62

# Ajuste automático si la tabla es demasiado ancha o alta
available_width <- 0.92 * diff(xlim_range)

table_y_top <- data_ylim[1] - 0.45 * table_gap
table_y_bottom_limit <- ylim_range[1] + 0.03 * diff(ylim_range)
available_height <- table_y_top - table_y_bottom_limit

repeat {
  max_width <- max(strwidth(txt, units = "user", cex = table_cex))
  line_h <- strheight("M", units = "user", cex = table_cex) * 1.45
  table_height <- length(txt) * line_h
  
  if (
    (max_width <= available_width && table_height <= available_height) ||
    table_cex <= 0.35
  ) {
    break
  }
  
  table_cex <- table_cex * 0.95
}

# Recalcular altura final de la tabla
line_h <- strheight("M", units = "user", cex = table_cex) * 1.45
table_height <- length(txt) * line_h

table_x <- mean(xlim_range)
table_y <- table_y_top

# Fondo blanco para que la tabla se lea bien
rect(
  xleft   = xlim_range[1] + 0.03 * diff(xlim_range),
  xright  = xlim_range[2] - 0.03 * diff(xlim_range),
  ybottom = table_y - table_height - 0.02 * diff(ylim_range),
  ytop    = table_y + 0.02 * diff(ylim_range),
  col     = adjustcolor("white", alpha.f = 0.92),
  border  = "gray80",
  lwd     = 0.7
)

text(
  x = table_x,
  y = table_y,
  labels = paste(txt, collapse = "\n"),
  adj = c(0.5, 1),
  family = "mono",
  cex = table_cex
)

par(op)
