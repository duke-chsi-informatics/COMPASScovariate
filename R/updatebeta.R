updatebeta = function(gammak, X1, nu, u, lambda2, tau2, W, K1,p1, sig2k){
  ### generate beta from MN
  beta_new = array(0, dim=c(K1,p1))
  p = p1-1
  K = K1+1
  for (k in 1:K1){
    Omegak = diag(W[,k])
    Bk_inv = diag(1/c(sig2k[k], lambda2*tau2[k]))
    tmp = t(X1)%*%Omegak%*%X1 + Bk_inv
    cm = chol(tmp)
    Vk = Matrix::chol2inv(cm) # posterior variance matrix 
    kappak = gammak[,k] - 0.5
    mk = Vk%*%(t(X1)%*%kappak)
    beta_new[k,] = mvtnorm::rmvnorm(1, mk, Vk)
  }
  
  ### updating the hyperparameters
  beta2 = 0.5*(beta_new[,-1])^2
  ## lambda2
  tau2M = (matrix(replicate(p,tau2),nrow=K1))
  lambda2 = invgamma::rinvgamma(p, shape = 0.5*K, rate = 1/nu + colSums(beta2/tau2M))
  ## tau2
  lambda2M = t(matrix(replicate(K1,lambda2),nrow=p))
  tau2 = invgamma::rinvgamma(K1, shape = 0.5*p1, rate = 1/u + rowSums(beta2/lambda2M))
  ## nu
  #set.seed(123)
  nu = invgamma::rinvgamma(p, shape = 1, rate = 1+1/lambda2)
  #set.seed(123)
  #sapply(1+1/lambda2, function(x) invgamma::rinvgamma(1, shape = 1, rate = x))
  ## u
  u = invgamma::rinvgamma(K1, shape = 1, rate = 1+1/tau2)
  ## W
  W = matrix(sapply(as.vector(X1%*%t(beta_new)), function(x) pgdraw::pgdraw(1,x)), nrow = I)
  
  return(list("W" = W, "u" = u, "nu" = nu, "tau2" = tau2, "lambda2" = lambda2, "beta" = beta_new ))
}
