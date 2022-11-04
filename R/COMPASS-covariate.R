.COMPASS.covariate <- function(n_s, n_u, iterations, replications, X,
                              verbose=TRUE, ...) {
  
  vmessage <- function(...)
    if (verbose) message(...) else invisible(NULL)
  
  ## Initial parameters
  vmessage("Initializing parameters...")
  N <- iterations ## The number of iterations to run the model
  ttt <- 2000 ## The step size in mode search (fixed)
  SS <- 1L
  
  N_s <- rowSums(n_s)
  N_u <- rowSums(n_u)
  
  I <- nrow(n_s)
  K <- ncol(n_u)
  K1 <- K - 1L
  p <- ncol(X) # number of features
  p1 = p + 1L # with intercept
  X1 = cbind(1, X) # feature matrix with intercept X = Ixp
  init_with_fisher = FALSE
  if (!init_with_fisher) {
    indi = array(1, dim = c(I, K)) # 0 indicate that gamma_ik=0
    for (k in 1:K1) {
      #non-responders when p_u >= p_s
      l2 = which((n_s[, k] / N_s) - (n_u[, k] / N_u) <= 0)
      indi[l2, k] <- 0;
    }
  }
  #alternately initialize indicators from Fisher's test.
  if (init_with_fisher) {
    indi = array(0, dim = c(I, K)) # 0 indicate that gamma_ik=0
    for (k in 1:K1) {
      for (i in 1:I) {
        indi[i, k] <- as.integer(fisher.test(matrix(c(n_s[i, k], N_s[i], n_u[i, k], N_u[i]), ncol = 2, nrow = 2))$conf.int[1] > 1)
      }
    }
  }
  indi[, K] <- rowSums(indi[, 1:K1, drop = FALSE])
  indi <- matrix(as.integer(indi), nrow = I)
  
  
  #############################################
  mk = array(as.integer(0), dim = c(1,K1));
  Istar = 0;
  mKstar = 0;
  
  gamma =array(as.integer(0), dim=c(I,K,N));
  
  alpha_u = array(0, dim=c(N,K));
  alpha_s = array(0, dim=c(N,K));
  
  varp_s1 = array(sqrt(5),dim=c(K,1)); # sqrt(var)
  varp_s2 = array(sqrt(15),dim=c(K,1)); # sqrt(var)
  varp_s2[K] = sqrt(25);
  varp_s1[K] = sqrt(15);
  
  pvar_s = array(0.8,dim=c(K,1));
  pvar_s[K] = 0.6;
  varp_u = array(sqrt(10),dim=c(K,1));
  
  pp = array(0.65, dim = c(I, 1))
  pb1 <- clamp(1.5 / median(indi[, K]), 0, 0.9)
  pb2 <- clamp(5.0 / median(indi[, K]), 0, 0.9)
  lambda_s = rep(0, K);  
  lambda_s[1:K1] = (10 ^ -2) * max(N_s, N_u)
  lambda_s[K] = max(N_s, N_u) - sum(lambda_s[1:K1])
  lambda_u = lambda_s
  
  alpha_u[1, 1:(K - 1)] = 10 #initializaion 
  alpha_u[1, K] = 150
  
  alpha_s[1, 1:(K - 1)] = 10 #initialization 
  alpha_s[1, K] = 100
  
  #### related to X
  beta =  array(0, dim=c(K1,p1, N)) #K-1 subsets, p features + 1 intercept
  nu =  invgamma::rinvgamma(p, shape = 0.5, rate = 1)
  u = invgamma::rinvgamma(K1, shape = 0.5, rate = 1)
  lambda2 = invgamma::rinvgamma(p, shape = 0.5, rate = 1/nu)
  tau2 = invgamma::rinvgamma(K1, shape = 0.5, rate = 1/u)
  W =  matrix(sapply(1:(I*K1), function(x) pgdraw::pgdraw(1,0)), nrow = I)
  sig2k = rep(10, K1) # var for normal prior of intercept
  #################### acceptance rate ###########################
  A_gm = array(as.integer(0), dim=c(I,N));
  A_alphau = array(as.integer(0), dim=c(K,N));
  A_alphas = array(as.integer(0), dim=c(K,N));
  
  vmessage("Computing initial parameter estimates...")
  for (tt in 2:N) {
    
    if (tt %% 1000 == 0) vmessage("Iteration ", tt, " of ", iterations, ".")
    
     # update alphau
     res2 <- updatealphau_noPu_Exp(alphaut = alpha_u[tt-1,],n_s = n_s,n_u=n_u, I=I, K=K, lambda_u = lambda_u, var_p = varp_u, ttt = ttt,gammat =gamma[,,tt-1])
     alpha_u[tt,] = res2$alphau_tt;
     A_alphau[,tt] = res2$Aalphau; 

    #update gamma_covariates
    WW = -X1%*%t(beta[,,tt-1])#IxK1 # Xbeta
    logWK1 = -log(1+exp(WW)) #log(pr(gamma_ik) = 1)
    logWK0 = WW+logWK1
    res1 <- updategammak_cov(n_s = n_s,n_u=n_u,gammat = gamma[,,tt-1],I=I,K=K,SS = SS,alphau = alpha_u[tt,],alphas = alpha_s[tt-1,],alpha=1,mk=mk,Istar = Istar,
                              mKstar = mKstar,pp=pp, pb1 = pb1, pb2 = pb2, indi=indi, WK1 = logWK1,
                             WK0 = logWK0)
    gamma[,,tt] = res1$gamma_tt;
    A_gm[,tt] = res1$Ag;
    Istar = res1$xIstar;
    mk = res1$xmk;
    mKstar = res1$mKstar;
    
    #update alphas
    res3 <- updatealphas_Exp(alphast = alpha_s[tt-1,], n_s = n_s,  K=K, I=I, lambda_s = lambda_s, gammat =gamma[,,tt], var_1 = varp_s1,var_2 = varp_s2,p_var = pvar_s,ttt = ttt)
    alpha_s[tt,] = res3$alphas_tt;
    A_alphas[,tt] = res3$Aalphas;

    # update beta
    res4 = updatebeta(gamma[,-K,tt], X1, nu, u, lambda2, tau2, W, K1,p1, sig2k)
    W = res4$W
    u = res4$u
    nu = res4$nu
    tau2 = res4$tau2
    lambda2 = res4$lambda2
    beta[,,tt] = res4$beta
    ####
    if (tt>4000&&(tt %% 4000)==0) {
      intr = (tt-1000+1):tt
      for (kk in 1:K) {
        Au = mean(A_alphau[kk,intr])
        if (Au<0.001) {varp_u[kk] = varp_u[kk]*sqrt(0.1);}
        else if (Au <0.05) {varp_u[kk] = varp_u[kk]*sqrt(0.5);}
        else if (Au < 0.2) {varp_u[kk] = varp_u[kk]*sqrt(0.9);}
        else if (Au>0.5) {varp_u[kk] = varp_u[kk]*sqrt(1.1);}
        else if (Au>0.75) {varp_u[kk] = varp_u[kk]*sqrt(2);}
        else if (Au>0.95) {varp_u[kk] = varp_u[kk]*sqrt(10);}
        As = mean(A_alphas[kk,intr])
        if (As<0.001) {varp_s1[kk] = varp_s1[kk]*sqrt(0.1);}
        else if (As <0.05) {varp_s1[kk] = varp_s1[kk]*sqrt(0.5);}
        else if (As < 0.2) {varp_s1[kk] = varp_s1[kk]*sqrt(0.9);}
        else if (As>0.5) {varp_s2[kk] = varp_s2[kk]*sqrt(1.1);}
        else if (As>0.75) {varp_s2[kk] = varp_s2[kk]*sqrt(2);}
        else if (As>0.95) {varp_s2[kk] = varp_s2[kk]*sqrt(10);}
      }
      for ( i in 1:I) {
        Agm = mean(A_gm[i,intr]);
        if (Agm<0.001) {pp[i] = min(0.9,pp[i]*1.4);}
        else if (Agm<0.05) {pp[i] = min(0.9,pp[i]*1.2);}
        else if (Agm<0.2) {pp[i] = min(0.9,pp[i]*1.1);}
        else if (Agm>0.6) {pp[i] = max(0.1,pp[i]*0.8);}
        else if (Agm>0.75) {pp[i] = max(0.1,pp[i]*0.5);}
        else if (Agm>0.95) {pp[i] = max(0.1,pp[i]*0.2);}
      }
    }
  }
  alpha_u[1,] = alpha_u[N,]
  gamma[,,1] = gamma[,,N]
  alpha_s[1,] = alpha_s[N,]
  A_gm = array(as.integer(0), dim=c(I,N));
  A_alphau = array(as.integer(0), dim=c(K,N));
  A_alphas = array(as.integer(0), dim=c(K,N));
  
  sNN=replications ## number of 'replications' -- should be user defined
  vmessage("Fitting model with ", sNN, " replications.")
  for (stt in 1:sNN) {
    vmessage("Running replication ", stt, " of ", sNN, "...")
    for (tt in 2:N) {
      # update alphau
      res2 <- updatealphau_noPu_Exp(alphaut = alpha_u[tt-1,],n_s = n_s,n_u=n_u, I=I, K=K, lambda_u = lambda_u, var_p = varp_u, ttt = ttt,gammat =gamma[,,tt-1])
      alpha_u[tt,] = res2$alphau_tt;
      A_alphau[,tt] = res2$Aalphau;

      #update gamma
      WW = -X1%*%t(beta[,,tt-1])
      logWK1 = -log(1+exp(WW)) #IxK1 # Xbeta
      logWK0 = WW+logWK1
      res1 <- updategammak_cov(n_s = n_s,n_u=n_u,gammat = gamma[,,tt-1],I=I,K=K,SS = SS,alphau = alpha_u[tt,],alphas = alpha_s[tt-1,],alpha=1,mk=mk,Istar = Istar,
                               mKstar = mKstar,pp=pp, pb1 = pb1, pb2 = pb2, indi=indi, WK1 = logWK1,
                               WK0 = logWK0)
      gamma[,,tt] = res1$gamma_tt;
      A_gm[,tt] = res1$Ag;
      Istar = res1$xIstar;
      mk = res1$xmk;
      mKstar = res1$mKstar;
      
      # update alphas
      res3 <- updatealphas_Exp(alphast = alpha_s[tt-1,], n_s = n_s,  K=K, I=I, lambda_s = lambda_s, gammat =gamma[,,tt], var_1 = varp_s1,var_2 = varp_s2,p_var = pvar_s,ttt = ttt)
       
      alpha_s[tt,] = res3$alphas_tt;
      A_alphas[,tt] = res3$Aalphas;
      
      # update beta
      res4 = updatebeta(gamma[,-K,tt], X1, nu, u, lambda2, tau2, W, K1,p1, sig2k)
      W = res4$W
      u = res4$u
      nu = res4$nu
      tau2 = res4$tau2
      lambda2 = res4$lambda2
      beta[,,tt] = res4$beta
      if (tt %% 1000 == 0) vmessage("Iteration ", tt, " of ", iterations, ".")
    }
    if (stt == sNN) {break}
    alpha_u[1,] = alpha_u[N,];
    gamma[,,1] = gamma[,,N];
    alpha_s[1,] = alpha_s[N,];
    A_gm = array(as.integer(0), dim=c(I,N));
    A_alphau = array(as.integer(0), dim=c(K,N));
    A_alphas = array(as.integer(0), dim=c(K,N));
  }
  
  ######################################
  Nburn=0;
  Mgamma = mat.or.vec(I,K);
  Nseq = seq(Nburn+1,N,by=1)
  for (ttt in Nseq) {
    Mgamma = Mgamma + gamma[,,ttt]; #thining
  }
  Mgamma = Mgamma/(N-Nburn);
  
  vmessage("Done!")
  
  
  ## set names on the output
  output <- list(
    alpha_s=alpha_s,
    A_alphas=rowMeans(A_alphas),
    alpha_u=alpha_u,
    A_alphau=rowMeans(A_alphau),
    beta = beta,
    gamma=gamma,
    mean_gamma=Mgamma,
    A_gamma=rowMeans(A_gm),
    model="covariate"
  )
  
  return(output)
  
}
