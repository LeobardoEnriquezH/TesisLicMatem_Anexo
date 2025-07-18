###Source: Arman Oganisian https://github.com/stablemarkets
## Load Packages
library(rstan)
library(LaplacesDemon)
library(latex2exp)
set.seed(340)


####------------------------ Data ---------------------------------####
N = 276 # #Num of municipalities (sample size)
warmup = 1000
iter = 2000
n_draws = iter - warmup

#import confounder (L), dose (A), and outcome (Y)
dat<-read.csv("data_partial_pool_CATE.csv")

W=dat$W #Informal occupied in baseline year
M=dat$M #State 19 Nuevo Leon (state 4 in this case)
S=dat$S #IMSS or ISSSTE affiliation increment in the period of study
R=dat$R #PEA increment in the period 

Wmat = cbind(1, W, M, S, R)

V=dat$V #Th states: 1:Baja Calif, 2:Chih, 3:Coah, 4:NL, 5:Sonora, 6:Tamaulipas

A=dat$A #Treatment indicator

Y=dat$Y #Outcome: informal employment increment

## sample size in each stratum
n_v = as.numeric(table(V))

## get indices of each stratum
ind_list = lapply(sort(unique(V)), function(v) c(1:N)[V==v] )

stan_data = list(Y=Y[order(V)], A=A[order(V)], V=V[order(V)], 
                 W = Wmat[order(V), ],
                 Pw = ncol(Wmat), 
                 Pv = length(unique(V)), 
                 N=N, n_v=n_v,
                 ind = c(0, cumsum(n_v)))

####------------------------ Sample Posterior    ---------------------------####
partial_pool_model = stan_model(file = "partial_pool.stan")

stan_res = sampling(partial_pool_model, data = stan_data, 
                    pars = c("odds_ratio", 'mu', "overall_odds_ratio" ),
                    warmup = warmup, iter = iter, chains=1, seed=1)

Psi_draws = extract(stan_res, pars='odds_ratio')[[1]]
overall_mean = extract(stan_res, pars='overall_odds_ratio')[[1]]

####------------------- Compute Frequentist Estimates   --------------------####

vfac = factor(V)

## estimate outcome model, which we'll then integrate over p(W)
freq_reg = glm(Y ~  W + vfac + A + vfac*A, family=binomial('logit') )

## loop through strata of interest
Psi_freq = numeric(6) ## shell to contain causal odds ratios
for(vf in 1:6){
  vval= as.factor(vf)
  
  ## standardize over empirical distribution p(W)
  p1 = predict(freq_reg, data.frame(W, vfac=vval, A=1 ), type='response')
  p0 = predict(freq_reg, data.frame(W, vfac=vval, A=0 ), type='response')
  
  ## take mean over empirical distribution among V=v
  marg_mean_y1 = mean( p1[V==vf] )
  marg_mean_y0 = mean( p0[V==vf] )
  
  ## compute causal odds ratio
  Psi_freq[vf] = (marg_mean_y1/(1-marg_mean_y1)) / (marg_mean_y0/(1-marg_mean_y0))
}


####-------------------         Plot Results            --------------------####
v_strata = 1:length(unique(V))
post_mean = colMeans(Psi_draws)

#png("ppooling_plot.png", width = 600, height = 500)

plot(post_mean, pch=20, col='blue', ylim=c(0,3),
     ylab=TeX("$\\Psi(v)$"), xlab=TeX("$V$"), axes=F )

axis_labs = paste0(v_strata, "\n(n=",table(V),")")

axis(1, at = 1:6, labels = axis_labs, padj = .5 )
axis(2, at = 0:6, labels = 0:6 )


### Plot posterior credible Intervals 
colfunc <- colorRampPalette(c("white", "lightblue"))
#ci_perc = seq(.95,.05,-.05)
ci_perc = c(.95,.90,.85)
colvec = colfunc(length(ci_perc))
names(colvec) = ci_perc

for(i in ci_perc){
  pci = apply(Psi_draws, 2, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
  segments(1:6,pci[1,], 1:6, pci[2,], col=colvec[as.character(i)], lwd=10 )
}
###

points(post_mean, col='steelblue', pch=20, lwd=8)
points(Psi_freq, col='black', pch=20, lwd=5)
abline(h= mean(overall_mean), col='steelblue', lty=2)
abline(h= 1, col='red', lty=2)

legend('topleft', 
       legend = c('Posterior Mean/Credible Band', 'MLE', "Overall Effect" , "No Effect"),
       col = c('steelblue', 'black', 'steelblue','red' ), pch=c(20,20,NA, NA),
       lty = c(NA,NA,2,2), bty='n')

dev.off()

