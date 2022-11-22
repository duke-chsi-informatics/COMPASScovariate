#Source for now but will want to call using COMPASS.R eventually
source("~/project_repos/COMPASScovariate/R/COMPASS-covariate.R")
source("~/project_repos/COMPASScovariate/R/updatebeta.R")

library(COMPASS)

I <- 200 ## sample size i.e., #subjects; I = 30
K <- 10 ## number of cell categories: K = 16
K1 = K-1
p <- 5 # number of covariates # larger p

set.seed(123)
gamma <- matrix(0, nrow = I, ncol =  K) ## 0/1 response

X <- matrix(rnorm(I*p), ncol = p, nrow = I) # covariates
mu_true = rnorm(K, 0, 1) # the true intercept
alpha_true = matrix(rbinom(p*K1, 1, 0.7), nrow = K1) # variable selection indicators

B_true = matrix(rnorm(K1*p, 0, 5), nrow = K1)
Balpha = B_true*alpha_true # coefficients with zeros


for(k in 1:K1){
  z = mu_true[k] + X%*%(Balpha[k,]) # linear combination with a bias
  pr = 1/(1+exp(-z))   # pass through an inv-logit function
  gamma[,k] = unlist(lapply(pr, function(x) rbinom(1,1, x)))
}

gamma[,K] = rbinom(I,1,0.5)
gamma_true = gamma

pu_true = MCMCpack::rdirichlet(I, rep(10, K1))
n_u = array(0, dim = c(I, K))
Nu = ceiling(rnorm(I,20000, 10))
for(i in 1:I){
  n_u[i,1:K1] = rmultinom(1, 5000, pu_true[i,])
}
n_u[,K] = Nu - rowSums(n_u)
ps_true = pu_true

for(i in 1:I){
  l1 = which(gamma[i,1:K1]==1)
  l0 = which(gamma[i,1:K1]==0)
  ps_true[i,l1] = MCMCpack::rdirichlet(1, rep(20, length(l1)))*sum(pu_true[i,l1])
}

n_s = array(0, dim = c(I, K))
Ns = ceiling(rnorm(I,30000, 10))
for(i in 1:I){
  n_s[i,1:K1] = rmultinom(1, 6000, ps_true[i,])
}
n_s[,K] = Ns - rowSums(n_s)

N = 10000
K1 <- K - 1L
p1 = p + 1L #
beta =  array(0, dim=c(K1,p1, N)) #K-1 subsets, p features + 1 intercept
nu =  invgamma::rinvgamma(p, shape = 0.5, rate = 1)
u = invgamma::rinvgamma(K1, shape = 0.5, rate = 1)
lambda2 = invgamma::rinvgamma(p, shape = 0.5, rate = 1/nu)
tau2 = invgamma::rinvgamma(K1, shape = 0.5, rate = 1/u)
W =  matrix(sapply(1:(I*K1), function(x) pgdraw::pgdraw(1,0)), nrow = I)
sig2k = rep(10, K1) # var for normal prior of intercept
X1 = cbind(1,X)



fit_cov = .COMPASS.covariate(n_s, n_u, 100000, 10, X)
for(tt in 2:N){
  res4 = updatebeta(gamma, X1, nu, u, lambda2, tau2, W, K1,p1, sig2k)
  W = res4$W
  u = res4$u
  nu = res4$nu
  tau2 = res4$tau2
  lambda2 = res4$lambda2
  beta[,,tt] = res4$beta
}
