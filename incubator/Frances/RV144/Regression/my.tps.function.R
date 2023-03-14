my.tps.function<-function(y,category=FALSE,scale=TRUE,dat,interaction=FALSE)
{
  if(scale)
  {
    y<-scale(y)
  }
  if(category)
  {
    y<-apply(y,2,function(x){x<-factor(cut(x,breaks=quantile(x,c(1,2/3,1/3,0)),include.lowest=T));levels(x)<-c("0","1","2");x})
  }
  
  res<-matrix(0,ncol(y),5)
  pval.int<-rep(NA,ncol(y))
  colnames(res)<-c("risk","CI.low","CI.up","p-value","t-stat")
  for(i in 1:ncol(y))
  {
    outcome<-y[,i]
    if(category==TRUE)
    {
      outcome<-as.factor(outcome)
    }
    
    kp <- dat$trt=="VACCINE" & !is.na(IgAprim)
    if(!interaction)
    {
      fit.tps <- try(tps(flrstatus ~ outcome + IgAprim + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE)  
    }else{
      fit.tps <- try(tps(flrstatus ~ outcome + IgAprim + outcome*IgAprim + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE)
      #print(fit.tps)
      temp <- 2*(1-pnorm(abs(fit.tps$coef[7]/sqrt(fit.tps$covm[7,7])))) 
      pval.int[i] <- ifelse(temp > 1.0, 1.0, temp) #  two-sided p-value for interaction        
    }
    if(class(fit.tps)!="try-error")
    {
      res[i,] <- as.matrix(cbind(round(exp(fit.tps$coef[2]),2),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),2),
                                 round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),2),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0)
                                                                                                    ,4),fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))
      if(category)
      {
        res[i,] <- as.matrix(cbind(round(1/exp(fit.tps$coef[2]),2),1/round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),2),
                                   1/round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),2),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0)
                                                                                                        ,4),fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))
      }
    }else
    {
      res[i,]<-rep(1,5)
    }
  }
  rownames(res)<-colnames(y)
  cbind(res,pval.int)
}


tps.multivariate<-function(y)
{
fit.tps <- tps(flrstatus ~ y + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE)
# x <- as.matrix(cbind(round(exp(fit.tps$coef[2]),2),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),2),
#                      round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),2),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0)
#                                                                                         ,4)))
# colnames(x)<-c("OR","CI.low","CI.up","p-value")
# rownames(x)<-colnames(y)

# x <- cbind(round(exp(fit.tps$coef),2),round(exp(fit.tps$coef - sqrt(diag(fit.tps$covm))*1.96),2),
#            round(exp(fit.tps$coef + sqrt(diag(fit.tps$covm))*1.96),2),round(pval2,4))

# Calculate multivariate p-value for all included coefficients
x <- x[2:(ncol(y)+1),]
temp <- t(fit.tps$coef[2:(ncol(y)+1)])%*%solve(fit.tps$covm[c(2:(ncol(y)+1)),c(2:(ncol(y)+1))])%*%fit.tps$coef[2:(ncol(y)+1)]
pval <- 1 - pchisq(temp, ncol(y))
return(pval)
}


tps.univariate<-function(y, IgA.include=FALSE)
{

  if(IgA.include){
    fit.tps <- tps(flrstatus ~ y + sex + risk.medium + risk.high + IgA, nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE)
  }
  else
  {
    fit.tps <- tps(flrstatus ~ y + sex + risk.medium + risk.high, nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE)
  }

x <- as.matrix(cbind(round(exp(fit.tps$coef[2]),2),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),2),
                     round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),2),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0)
                                                                                        ,4)))
colnames(x)<-c("OR","CI.low","CI.up","p-value")
rownames(x)<-colnames(y)
return(x)
}


tps.univariate.interaction<-function(y,z, IgA.include=FALSE)
{
  if(IgA.include==FALSE){
    fit.tps <- tps(flrstatus ~ y*z + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE)
    x <- as.matrix(cbind(round(exp(fit.tps$coef[2]),2),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),2),
                         round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),2),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0)
                                                                                            ,4)))
    temp <- 2*(1-pnorm(abs(fit.tps$coef[7]/sqrt(fit.tps$covm[7,7])))) 
    pval.int <- ifelse(temp > 1.0, 1.0, temp) #  two-sided p-value for interaction
    return(pval.int)
  }
  else{
    fit.tps <- tps(flrstatus ~ y*z + IgA + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE)
    temp <- 2*(1-pnorm(abs(fit.tps$coef[8]/sqrt(fit.tps$covm[8,8])))) 
    pval.int <- ifelse(temp > 1.0, 1.0, temp) #  two-sided p-value for interaction
    return(pval.int)
  }
}

tps.univariate.tert<-function(y)
{
  y<-apply(y,2,function(x){x<-factor(cut(x,breaks=quantile(x,c(1,2/3,1/3,0)),include.lowest=T));levels(x)<-c("0","1","2");x})
  
  fit.tps <- tps(flrstatus ~ as.factor(y) + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE)

  temp <- t(fit.tps$coef[2:3])%*%solve(fit.tps$covm[c(2:3),c(2:3)])%*%fit.tps$coef[2:3]
  pval <- 1 - pchisq(temp, 2) # overall p-value for tertiles
  
  #temp <- 2*(1-pnorm(abs(fit.tps$coef/sqrt(diag(fit.tps$covm))))) 
  #pval2 <- ifelse(temp > 1.0, 1.0, temp)[-1] #  two-sided p-values for all covariates excluding intercept
  
  #x <- as.matrix(cbind(round(exp(fit.tps$coef[-1]),2),round(exp(fit.tps$coef[-1] - sqrt(diag(fit.tps$covm)[-1])*1.96),2), round(exp(fit.tps$coef[-1] + sqrt(diag(fit.tps$covm)[-1])*1.96),2),round(pval2,4)))
  x <- as.matrix(cbind(0,0,0,round(pval,4)))
                       #,round(c(rep(pval,2),pval2[(length(pval2)-2):length(pval2)]),4)))
  colnames(x)<-c("OR","CI.low","CI.up","p-value")
  rownames(x)<-colnames(y)  
  return(x)
}


tps.univariate.matrix<-function(y, IgA.include=FALSE)
{
  n<-ncol(y)
  res<-matrix(0,4,n)
  for(i in 1:n)
  {
    res[,i]<-tps.univariate(y[,i,drop=FALSE], IgA.include=IgA.include)
  }
  colnames(res)<-colnames(y)
  rownames(res)<-c("OR","CI.low","CI.up","p-value")
  res
}

tps.univariate.matrix.tertile<-function(y)
{
  n<-ncol(y)
  res<-matrix(0,4,n)
  for(i in 1:n)
  {
    res[,i]<-tps.univariate.tert(y[,i,drop=FALSE])
  }
  colnames(res)<-colnames(y)
  rownames(res)<-c("OR","CI.low","CI.up","p-value")
  res
}

tps.interaction.table<-function(y,varI)
{
  fit.tps <- tps(flrstatus ~ y*varI + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE)
  varI.levels <- quantile(varI,probs=c(.2,.5,.8))
  logOR.1 <- fit.tps$coef[2] + fit.tps$coef[7]*varI.levels[1]
  var.1 <- fit.tps$covm[2,2] + fit.tps$covm[7,7]*varI.levels[1]^2 + fit.tps$covm[2,7]*2*varI.levels[1]
  logOR.2 <- fit.tps$coef[2] + fit.tps$coef[7]*varI.levels[2]
  var.2 <- fit.tps$covm[2,2] + fit.tps$covm[7,7]*varI.levels[2]^2 + fit.tps$covm[2,7]*2*varI.levels[2]
  logOR.3 <- fit.tps$coef[2] + fit.tps$coef[7]*varI.levels[3]
  var.3 <- fit.tps$covm[2,2] + fit.tps$covm[7,7]*varI.levels[3]^2 + fit.tps$covm[2,7]*2*varI.levels[3]
  logOR <- c(logOR.1,logOR.2,logOR.3)
  #print(logOR)
  varlogOR <- c(var.1,var.2,var.3)
  ll <- logOR - 1.96*sqrt(varlogOR)
  ul <- logOR + 1.96*sqrt(varlogOR)
  temp <- 2*(1-pnorm(abs(logOR/sqrt(varlogOR)))) 
  pval <- ifelse(temp > 1.0, 1.0, temp)  
  x <- data.frame(cbind(round(exp(logOR),2),round(exp(ll),2),round(exp(ul),2),round(pval,4)))
  names(x)<-c("OR","CI.low","CI.up","p-value")
  rownames(x)<-c("20%","50%","80%")
  return(x)
}

cch.univariate<-function(y,pd)
{
  fit <- cch(Surv(flrtime,flrstatus) ~ y + sex + risk.medium + risk.high,
             subcoh=in.subcohort,stratum=stratuminds,
             id=pd$ptid,cohort.size=cohortstratasizes,method="II.Borgan")

  x <- summary(fit)$coeff
  x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
             round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
  x<-x[1,,drop=FALSE]
  colnames(x)<-c("OR","CI.low","CI.up","p-value")
  rownames(x)<-colnames(y)
  return(x)
}


cch.univariate.matrix<-function(y,pd)
{
  n<-ncol(y)
  res<-matrix(0,4,n)
  for(i in 1:n)
  {
    res[,i]<-cch.univariate(y[,i,drop=FALSE],pd)
  }
  colnames(res)<-colnames(y)
  rownames(res)<-c("HR","CI.low","CI.up","p-value")
  res
}

cch.univariate.interaction<-function(y,z,pd)
{
  
  fit <- cch(Surv(flrtime,flrstatus) ~ y*z + sex + risk.medium + risk.high, subcoh=in.subcohort,stratum=stratuminds,id=pd$ptid,cohort.size=cohortstratasizes,method="II.Borgan")
    
  return(summary(fit)$coefficients["y:z","p"])
}
