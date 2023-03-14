##################################
# 1. univariate for each 6 primary and 8 toggle (quantitative) -- exclude toggle part from final report
# 2. univariate for each 6 primary and 8 toggle (tertitle) -- exclude toggle part from final report
# 3. multivariate for the 6 primary (quantitative)
# 4. multivariate for the 8 toggle (quantitative) -- exclude from final report
# 5. multivariate for the 5 primary + 1 toggle (quantitative)
# 6. multivariate for the 6 primary (tertile)
# 7. multivariate for 5 primary + 1 toggle (tertile)
##################################


#########
# 9/2/11 cum  inc code modified to use sampling wts as defined using controls only
#########
 

#############################################
# Univariate quantitative Cox model analyses:
#############################################
comment<-FALSE
if(comment)
{
kp <- dat$trt=="VACCINE" & !is.na(IgAprim)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
xIgAprim <- as.matrix(t(x[1,]))
pvalIgAprim <- x[1,4]

kp <- dat$trt=="VACCINE" & !is.na(Avidprim)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ Avidprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
xAvidprim <- as.matrix(t(x[1,]))
pvalAvidprim <- x[1,4]

kp <- dat$trt=="VACCINE" & !is.na(ADCCprim)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ ADCCprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
xADCCprim <- as.matrix(t(x[1,]))
pvalADCCprim <- x[1,4]

kp <- dat$trt=="VACCINE" & !is.na(NAbprim)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ NAbprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
xNAbprim <- as.matrix(t(x[1,]))
pvalNAbprim <- x[1,4]

kp <- dat$trt=="VACCINE" & !is.na(V2prim)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ V2prim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
xV2prim <- as.matrix(t(x[1,]))
pvalV2prim <- x[1,4]

kp <- dat$trt=="VACCINE" & CD4ICSprim_imp==0 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
xCD4prim <- as.matrix(t(x[1,]))
pvalCD4prim <- x[1,4]

# Obtain FWER-adjusted p-values and q-values
pvals <- c(pvalIgAprim,pvalAvidprim,pvalADCCprim,pvalNAbprim,pvalV2prim,pvalCD4prim)
holmpvalues <- round(p.adjust(pvals,"holm"), 4)
qvalues <- round(p.adjust(pvals,"BH"), 4)

xIgAprim <- cbind(xIgAprim,holmpvalues[1],qvalues[1])
xAvidprim <- cbind(xAvidprim,holmpvalues[2],qvalues[2])
xADCCprim <- cbind(xADCCprim,holmpvalues[3],qvalues[3])
xNAbprim <- cbind(xNAbprim,holmpvalues[4],qvalues[4])
xV2prim <- cbind(xV2prim,holmpvalues[5],qvalues[5])
xCD4prim <- cbind(xCD4prim,holmpvalues[6],qvalues[6])

write(t(xIgAprim),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAprim.tex",ncolumns=6,sep=" & ")
write(t(xAvidprim),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidprim.tex",ncolumns=6,sep=" & ")
write(t(xADCCprim),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCprim.tex",ncolumns=6,sep=" & ")
write(t(xNAbprim),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbprim.tex",ncolumns=6,sep=" & ")
write(t(xV2prim),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2prim.tex",ncolumns=6,sep=" & ")
write(t(xCD4prim),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSprim.tex",ncolumns=6,sep=" & ")


kp <- dat$trt=="VACCINE" & !is.na(IgAtog)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAtog[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- as.matrix(t(x[1,]))
write(t(x),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAtog.tex",ncolumns=4,sep=" & ")

kp <- dat$trt=="VACCINE" & !is.na(Avidtog)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ Avidtog[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- as.matrix(t(x[1,]))
write(t(x),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidtog.tex",ncolumns=4,sep=" & ")

kp <- dat$trt=="VACCINE" & !is.na(ADCCtog)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ ADCCtog[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- as.matrix(t(x[1,]))
write(t(x),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCtog.tex",ncolumns=4,sep=" & ")

kp <- dat$trt=="VACCINE" & !is.na(NAbAtog)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ NAbAtog[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- as.matrix(t(x[1,]))
write(t(x),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbAtog.tex",ncolumns=4,sep=" & ")

kp <- dat$trt=="VACCINE" & !is.na(NAbStog)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ NAbStog[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- as.matrix(t(x[1,]))
write(t(x),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbStog.tex",ncolumns=4,sep=" & ")

kp <- dat$trt=="VACCINE" & !is.na(NAbMtog)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ NAbMtog[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- as.matrix(t(x[1,]))
write(t(x),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbMtog.tex",ncolumns=4,sep=" & ")

kp <- dat$trt=="VACCINE" & !is.na(V2tog)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ V2tog[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- as.matrix(t(x[1,]))
write(t(x),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2tog.tex",ncolumns=4,sep=" & ")

kp <- dat$trt=="VACCINE" & CD4lumtog_imp==0 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ CD4lumtog[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- as.matrix(t(x[1,]))
write(t(x),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4lumtog.tex",ncolumns=4,sep=" & ")

##############################################
# Univariate low/med/high (primary variables)
##############################################

kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pvalIgAprim <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
xIgAprim <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pvalIgAprim,nrow(x)),4))

kp <- dat$trt=="VACCINE" & !is.na(Avidprimtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(Avidprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pvalAvidprim <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
xAvidprim <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pvalAvidprim,nrow(x)),4))


kp <- dat$trt=="VACCINE" & !is.na(ADCCprimtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(ADCCprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pvalADCCprim <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
xADCCprim <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pvalADCCprim,nrow(x)),4))

kp <- dat$trt=="VACCINE" & !is.na(NAbprimtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(NAbprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pvalNAbprim <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
xNAbprim <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pvalNAbprim,nrow(x)),4))

kp <- dat$trt=="VACCINE" & !is.na(V2primtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(V2primtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pvalV2prim <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
xV2prim <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pvalV2prim,nrow(x)),4))

kp <- dat$trt=="VACCINE" & CD4ICSprim_imp==0 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pvalCD4prim <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
xCD4prim <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pvalCD4prim,nrow(x)),4))


# Obtain FWER-adjusted p-values and q-values
pvals <- c(pvalIgAprim,pvalAvidprim,pvalADCCprim,pvalNAbprim,pvalV2prim,pvalCD4prim)
holmpvalues <- round(p.adjust(pvals,"holm"), 4)
qvalues <- round(p.adjust(pvals,"BH"), 4)
M <- nrow(xIgAprim)
holmpvalues <-  matrix(rep(holmpvalues,M),nrow=M,byrow=T)
qvalues <-  matrix(rep(qvalues,M),nrow=M,byrow=T)

xIgAprimadj <- cbind(xIgAprim,holmpvalues[,1],qvalues[,1])
xAvidprimadj <- cbind(xAvidprim,holmpvalues[,2],qvalues[,2])
xADCCprimadj <- cbind(xADCCprim,holmpvalues[,3],qvalues[,3])
xNAbprimadj <- cbind(xNAbprim,holmpvalues[,4],qvalues[,4])
xV2primadj <- cbind(xV2prim,holmpvalues[,5],qvalues[,5])
xCD4primadj <- cbind(xCD4prim,holmpvalues[,6],qvalues[,6])


write(xIgAprim[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAprimtert.tex",ncolumns=5,sep=" & ")
write(xIgAprim[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAprimtert.tex",ncolumns=5,sep=" & ")
write(xAvidprim[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidprimtert.tex",ncolumns=5,sep=" & ")
write(xAvidprim[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidprimtert.tex",ncolumns=5,sep=" & ")
write(xADCCprim[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCprimtert.tex",ncolumns=5,sep=" & ")
write(xADCCprim[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCprimtert.tex",ncolumns=5,sep=" & ")
write(xNAbprim[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbprimtert.tex",ncolumns=5,sep=" & ")
write(xNAbprim[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbprimtert.tex",ncolumns=5,sep=" & ")
write(xV2prim[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2primtert.tex",ncolumns=5,sep=" & ")
write(xV2prim[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2primtert.tex",ncolumns=5,sep=" & ")
write(xCD4prim[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSprimtert.tex",ncolumns=5,sep=" & ")
write(xCD4prim[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSprimtert.tex",ncolumns=5,sep=" & ")

write(xIgAprimadj[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAprimtertadj.tex",ncolumns=7,sep=" & ")
write(xIgAprimadj[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAprimtertadj.tex",ncolumns=7,sep=" & ")
write(xAvidprimadj[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidprimtertadj.tex",ncolumns=7,sep=" & ")
write(xAvidprimadj[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidprimtertadj.tex",ncolumns=7,sep=" & ")
write(xADCCprimadj[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCprimtertadj.tex",ncolumns=7,sep=" & ")
write(xADCCprimadj[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCprimtertadj.tex",ncolumns=7,sep=" & ")
write(xNAbprimadj[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbprimtertadj.tex",ncolumns=7,sep=" & ")
write(xNAbprimadj[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbprimtertadj.tex",ncolumns=7,sep=" & ")
write(xV2primadj[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2primtertadj.tex",ncolumns=7,sep=" & ")
write(xV2primadj[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2primtertadj.tex",ncolumns=7,sep=" & ")
write(xCD4primadj[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSprimtertadj.tex",ncolumns=7,sep=" & ")
write(xCD4primadj[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSprimtertadj.tex",ncolumns=7,sep=" & ")


####################################
# Univariate low/med/high (toggles) 
####################################

kp <- dat$trt=="VACCINE" & !is.na(IgAtogtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAtogtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pval,nrow(x)),4))
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAtogtert.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAtogtert.tex",ncolumns=5,sep=" & ")

kp <- dat$trt=="VACCINE" & !is.na(Avidtogtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(Avidtogtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pval,nrow(x)),4))
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidtogtert.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidtogtert.tex",ncolumns=5,sep=" & ")

kp <- dat$trt=="VACCINE" & !is.na(ADCCtogtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(ADCCtogtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pval,nrow(x)),4))
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCtogtert.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCtogtert.tex",ncolumns=5,sep=" & ")


kp <- dat$trt=="VACCINE" & !is.na(NAbAtogtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(NAbAtogtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pval,nrow(x)),4))
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbAtogtert.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbAtogtert.tex",ncolumns=5,sep=" & ")


kp <- dat$trt=="VACCINE" & !is.na(NAbStogtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(NAbStogtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pval,nrow(x)),4))
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbStogtert.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbStogtert.tex",ncolumns=5,sep=" & ")


kp <- dat$trt=="VACCINE" & !is.na(NAbMtogtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(NAbMtogtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pval,nrow(x)),4))
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbMtogtert.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbMtogtert.tex",ncolumns=5,sep=" & ")


kp <- dat$trt=="VACCINE" & !is.na(V2togtert)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(V2togtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pval,nrow(x)),4))
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2togtert.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2togtert.tex",ncolumns=5,sep=" & ")

kp <- dat$trt=="VACCINE" & CD4lumtog_imp==0 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(CD4lumtogtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(rep(pval,nrow(x)),4))
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4lumtogtert.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4lumtogtert.tex",ncolumns=5,sep=" & ")


##############################
# Multivariate model
# Quantitative variables 
##############################


kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(Avidprim) & !is.na(ADCCprim) & !is.na(NAbprim) &
                           !is.na(V2prim) & !is.na(CD4ICSprim) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + Avidprim[kp] + ADCCprim[kp] + NAbprim[kp] 
                + V2prim[kp] + CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))

x <- x[1:6,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvprim.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvprim.tex",ncolumns=5,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvprim.tex",ncolumns=5,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvprim.tex",ncolumns=5,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvprim.tex",ncolumns=5,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvprim.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvprim.tex",ncolumns=1)


############################################
# Sensitivity analysis: multivariate model 
# Quantitative and categorical variables 
############################################

# Original multivariate quantitative model
kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(Avidprim) & !is.na(ADCCprim) & !is.na(NAbprim) &
                           !is.na(V2prim) & !is.na(CD4ICSprim)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + Avidprim[kp] + ADCCprim[kp] + NAbprim[kp]
                + V2prim[kp] + CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:6,]
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)

# Remove the variables nowhere near significant on the univariate analyses:
kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(V2prim) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + 
                + V2prim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:2,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvprimsens1.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvprimsens1.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvprimsens1.tex",ncolumns=1)


# Remove gender and risk category:
kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(V2prim)
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] +
                + V2prim[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:2,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval <- 1 - pchisq(temp, 2)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvprimsens2.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvprimsens2.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvprimsens2.tex",ncolumns=1)

# Original multivariate categorical model
kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbprimtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) + factor(Avidprimtert[kp]) + factor(ADCCprimtert[kp]) + factor(NAbprimtert[kp]) 
                + factor(V2primtert[kp]) + factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")

# Remove variables nowhere near significant (keep only IgA and V2)
kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbprimtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) +  
                + factor(V2primtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")

x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
	round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:4,]
qval <- rep(round(p.adjust(c(pval1,pval2),method="BH"),4),each=2)
x <- cbind(x,round(c(pval1,pval1,pval2,pval2),4),qval)
 
temp <- t(fit$coef[1:4])%*%solve(fit$var[c(1:4),c(1:4)])%*%fit$coef[1:4]
pval.global <- 1 - pchisq(temp, 4)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvprimtertsens1.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvprimtertsens1.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvprimtertsens1.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvprimtertsens1.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtertsens1.tex",ncolumns=1)

# Remove gender and risk category
kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbprimtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) +  
                + factor(V2primtert[kp]),
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")

x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
	round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:4,]
qval <- rep(round(p.adjust(c(pval1,pval2),method="BH"),4),each=2)
x <- cbind(x,round(c(pval1,pval1,pval2,pval2),4),qval)
 
temp <- t(fit$coef[1:4])%*%solve(fit$var[c(1:4),c(1:4)])%*%fit$coef[1:4]
pval.global <- 1 - pchisq(temp, 4)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvprimtertsens2.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvprimtertsens2.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvprimtertsens2.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvprimtertsens2.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtertsens2.tex",ncolumns=1)





################################
# Toggled multivariate model
# Quantitative variables 
################################

kp <- dat$trt=="VACCINE" & !is.na(IgAtog) & !is.na(Avidprim) & !is.na(ADCCprim) & !is.na(NAbprim) &
                           !is.na(V2prim) & !is.na(CD4ICSprim) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAtog[kp] + Avidprim[kp] + ADCCprim[kp] + NAbprim[kp] 
                + V2prim[kp] + CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:6,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog1.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog1.tex",ncolumns=5,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog1.tex",ncolumns=5,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog1.tex",ncolumns=5,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog1.tex",ncolumns=5,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog1.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog1.tex",ncolumns=1)



kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(Avidtog) & !is.na(ADCCprim) & !is.na(NAbprim) &
                           !is.na(V2prim) & !is.na(CD4ICSprim) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + Avidtog[kp] + ADCCprim[kp] + NAbprim[kp] 
                + V2prim[kp] + CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:6,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog2.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog2.tex",ncolumns=5,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog2.tex",ncolumns=5,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog2.tex",ncolumns=5,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog2.tex",ncolumns=5,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog2.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog2.tex",ncolumns=1)


kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(Avidprim) & !is.na(ADCCtog) & !is.na(NAbprim) &
                           !is.na(V2prim) & !is.na(CD4ICSprim) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + Avidprim[kp] + ADCCtog[kp] + NAbprim[kp] 
                + V2prim[kp] + CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:6,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog3.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog3.tex",ncolumns=5,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog3.tex",ncolumns=5,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog3.tex",ncolumns=5,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog3.tex",ncolumns=5,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog3.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog3.tex",ncolumns=1)


kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(Avidprim) & !is.na(ADCCprim) & !is.na(NAbAtog) &
                           !is.na(V2prim) & !is.na(CD4ICSprim) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + Avidprim[kp] + ADCCprim[kp] + NAbAtog[kp] 
                + V2prim[kp] + CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:6,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog4.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog4.tex",ncolumns=5,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog4.tex",ncolumns=5,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog4.tex",ncolumns=5,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog4.tex",ncolumns=5,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog4.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog4.tex",ncolumns=1)


kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(Avidprim) & !is.na(ADCCprim) & !is.na(NAbStog) &
                           !is.na(V2prim) & !is.na(CD4ICSprim) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + Avidprim[kp] + ADCCprim[kp] + NAbStog[kp] 
                + V2prim[kp] + CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:6,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog5.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog5.tex",ncolumns=5,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog5.tex",ncolumns=5,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog5.tex",ncolumns=5,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog5.tex",ncolumns=5,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog5.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog5.tex",ncolumns=1)


kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(Avidprim) & !is.na(ADCCprim) & !is.na(NAbMtog) &
                           !is.na(V2prim) & !is.na(CD4ICSprim) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + Avidprim[kp] + ADCCprim[kp] + NAbMtog[kp] 
                + V2prim[kp] + CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:6,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog6.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog6.tex",ncolumns=5,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog6.tex",ncolumns=5,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog6.tex",ncolumns=5,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog6.tex",ncolumns=5,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog6.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog6.tex",ncolumns=1)


kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(Avidprim) & !is.na(ADCCprim) & !is.na(NAbprim) &
                           !is.na(V2tog) & !is.na(CD4ICSprim) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + Avidprim[kp] + ADCCprim[kp] + NAbprim[kp] 
                + V2tog[kp] + CD4ICSprim[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:6,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog7.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog7.tex",ncolumns=5,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog7.tex",ncolumns=5,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog7.tex",ncolumns=5,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog7.tex",ncolumns=5,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog7.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog7.tex",ncolumns=1)


kp <- dat$trt=="VACCINE" & !is.na(IgAprim) & !is.na(Avidprim) & !is.na(ADCCprim) & !is.na(NAbprim) &
                           !is.na(V2prim) & !is.na(CD4lumtog) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ IgAprim[kp] + Avidprim[kp] + ADCCprim[kp] + NAbprim[kp] 
                + V2prim[kp] + CD4lumtog[kp] + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4))
x <- x[1:6,]
x <- cbind(x,round(p.adjust(x[,4],"BH"), 4)) # q-values
temp <- t(fit$coef[1:6])%*%solve(fit$var[c(1:6),c(1:6)])%*%fit$coef[1:6]
pval <- 1 - pchisq(temp, 6)
write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog8.tex",ncolumns=5,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog8.tex",ncolumns=5,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog8.tex",ncolumns=5,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog8.tex",ncolumns=5,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog8.tex",ncolumns=5,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog8.tex",ncolumns=5,sep=" & ")
write(round(pval,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog8.tex",ncolumns=1)

###########################################
### Multivariate low/med/high model
########################################### 

kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbprimtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) + factor(Avidprimtert[kp]) + factor(ADCCprimtert[kp]) + factor(NAbprimtert[kp]) 
                + factor(V2primtert[kp]) + factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[5:6])%*%solve(fit$var[c(5:6),c(5:6)])%*%fit$coef[5:6]
pval3 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[7:8])%*%solve(fit$var[c(7:8),c(7:8)])%*%fit$coef[7:8]
pval4 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[9:10])%*%solve(fit$var[c(9:10),c(9:10)])%*%fit$coef[9:10]
pval5 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[11:12])%*%solve(fit$var[c(11:12),c(11:12)])%*%fit$coef[11:12]
pval6 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),
round(c(pval1,pval1,pval2,pval2,pval3,pval3,pval4,pval4,pval5,pval5,pval6,pval6,round(x[,4],4)[13:15]),4))
qval <- rep(round(p.adjust(c(pval1,pval2,pval3,pval4,pval5,pval6),method="BH"),4),each=2)
x <- x[1:12,]
x <- cbind(x,qval)

temp <- t(fit$coef[1:12])%*%solve(fit$var[c(1:12),c(1:12)])%*%fit$coef[1:12]
pval.global <- 1 - pchisq(temp, 12)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvprimtert.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvprimtert.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvprimtert.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidmvprimtert.tex",ncolumns=6,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvprimtert.tex",ncolumns=6,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCmvprimtert.tex",ncolumns=6,sep=" & ")
write(x[7,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvprimtert.tex",ncolumns=6,sep=" & ")
write(x[8,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbmvprimtert.tex",ncolumns=6,sep=" & ")
write(x[9,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvprimtert.tex",ncolumns=6,sep=" & ")
write(x[10,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvprimtert.tex",ncolumns=6,sep=" & ")
write(x[11,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvprimtert.tex",ncolumns=6,sep=" & ")
write(x[12,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSmvprimtert.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtert.tex",ncolumns=1)

################################
# Toggled low/med/high model 
################################

kp <- dat$trt=="VACCINE" & !is.na(IgAtogtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbprimtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAtogtert[kp]) + factor(Avidprimtert[kp]) + factor(ADCCprimtert[kp]) + factor(NAbprimtert[kp]) 
                + factor(V2primtert[kp]) + factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[5:6])%*%solve(fit$var[c(5:6),c(5:6)])%*%fit$coef[5:6]
pval3 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[7:8])%*%solve(fit$var[c(7:8),c(7:8)])%*%fit$coef[7:8]
pval4 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[9:10])%*%solve(fit$var[c(9:10),c(9:10)])%*%fit$coef[9:10]
pval5 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[11:12])%*%solve(fit$var[c(11:12),c(11:12)])%*%fit$coef[11:12]
pval6 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(c(pval1,pval1,pval2,pval2,pval3,pval3,pval4,pval4,pval5,pval5,pval6,pval6,round(x[,4],4)[13:15]),4))
qval <- rep(round(p.adjust(c(pval1,pval2,pval3,pval4,pval5,pval6),method="BH"),4),each=2)
x <- x[1:12,]
x <- cbind(x,qval)

temp <- t(fit$coef[1:12])%*%solve(fit$var[c(1:12),c(1:12)])%*%fit$coef[1:12]
pval.global <- 1 - pchisq(temp, 12)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidmvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCmvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[7,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[8,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbmvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[9,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[10,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[11,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog1tert.tex",ncolumns=6,sep=" & ")
write(x[12,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSmvtog1tert.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog1tert.tex",ncolumns=1)


kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidtogtert) & !is.na(ADCCprimtert) & !is.na(NAbprimtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) + factor(Avidtogtert[kp]) + factor(ADCCprimtert[kp]) + factor(NAbprimtert[kp]) 
                + factor(V2primtert[kp]) + factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[5:6])%*%solve(fit$var[c(5:6),c(5:6)])%*%fit$coef[5:6]
pval3 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[7:8])%*%solve(fit$var[c(7:8),c(7:8)])%*%fit$coef[7:8]
pval4 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[9:10])%*%solve(fit$var[c(9:10),c(9:10)])%*%fit$coef[9:10]
pval5 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[11:12])%*%solve(fit$var[c(11:12),c(11:12)])%*%fit$coef[11:12]
pval6 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(c(pval1,pval1,pval2,pval2,pval3,pval3,pval4,pval4,pval5,pval5,pval6,pval6,round(x[,4],4)[13:15]),4))
qval <- rep(round(p.adjust(c(pval1,pval2,pval3,pval4,pval5,pval6),method="BH"),4),each=2)
x <- x[1:12,]
x <- cbind(x,qval)

temp <- t(fit$coef[1:12])%*%solve(fit$var[c(1:12),c(1:12)])%*%fit$coef[1:12]
pval.global <- 1 - pchisq(temp, 12)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidmvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCmvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[7,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[8,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbmvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[9,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[10,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[11,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog2tert.tex",ncolumns=6,sep=" & ")
write(x[12,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSmvtog2tert.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog2tert.tex",ncolumns=1)

kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCtogtert) & !is.na(NAbprimtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) + factor(Avidprimtert[kp]) + factor(ADCCtogtert[kp]) + factor(NAbprimtert[kp]) 
                + factor(V2primtert[kp]) + factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[5:6])%*%solve(fit$var[c(5:6),c(5:6)])%*%fit$coef[5:6]
pval3 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[7:8])%*%solve(fit$var[c(7:8),c(7:8)])%*%fit$coef[7:8]
pval4 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[9:10])%*%solve(fit$var[c(9:10),c(9:10)])%*%fit$coef[9:10]
pval5 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[11:12])%*%solve(fit$var[c(11:12),c(11:12)])%*%fit$coef[11:12]
pval6 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(c(pval1,pval1,pval2,pval2,pval3,pval3,pval4,pval4,pval5,pval5,pval6,pval6,round(x[,4],4)[13:15]),4))
qval <- rep(round(p.adjust(c(pval1,pval2,pval3,pval4,pval5,pval6),method="BH"),4),each=2)
x <- x[1:12,]
x <- cbind(x,qval)

temp <- t(fit$coef[1:12])%*%solve(fit$var[c(1:12),c(1:12)])%*%fit$coef[1:12]
pval.global <- 1 - pchisq(temp, 12)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidmvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCmvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[7,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[8,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbmvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[9,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[10,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[11,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog3tert.tex",ncolumns=6,sep=" & ")
write(x[12,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSmvtog3tert.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog3tert.tex",ncolumns=1)


kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbAtogtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAtogtert[kp]) + factor(Avidprimtert[kp]) + factor(ADCCprimtert[kp]) + factor(NAbAtogtert[kp]) 
                + factor(V2primtert[kp]) + factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[5:6])%*%solve(fit$var[c(5:6),c(5:6)])%*%fit$coef[5:6]
pval3 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[7:8])%*%solve(fit$var[c(7:8),c(7:8)])%*%fit$coef[7:8]
pval4 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[9:10])%*%solve(fit$var[c(9:10),c(9:10)])%*%fit$coef[9:10]
pval5 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[11:12])%*%solve(fit$var[c(11:12),c(11:12)])%*%fit$coef[11:12]
pval6 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(c(pval1,pval1,pval2,pval2,pval3,pval3,pval4,pval4,pval5,pval5,pval6,pval6,round(x[,4],4)[13:15]),4))
qval <- rep(round(p.adjust(c(pval1,pval2,pval3,pval4,pval5,pval6),method="BH"),4),each=2)
x <- x[1:12,]
x <- cbind(x,qval)

temp <- t(fit$coef[1:12])%*%solve(fit$var[c(1:12),c(1:12)])%*%fit$coef[1:12]
pval.global <- 1 - pchisq(temp, 12)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidmvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCmvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[7,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[8,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbmvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[9,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[10,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[11,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog4tert.tex",ncolumns=6,sep=" & ")
write(x[12,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSmvtog4tert.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog4tert.tex",ncolumns=1)

kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbStogtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) + factor(Avidprimtert[kp]) + factor(ADCCprimtert[kp]) + factor(NAbStogtert[kp]) 
                + factor(V2primtert[kp]) + factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[5:6])%*%solve(fit$var[c(5:6),c(5:6)])%*%fit$coef[5:6]
pval3 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[7:8])%*%solve(fit$var[c(7:8),c(7:8)])%*%fit$coef[7:8]
pval4 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[9:10])%*%solve(fit$var[c(9:10),c(9:10)])%*%fit$coef[9:10]
pval5 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[11:12])%*%solve(fit$var[c(11:12),c(11:12)])%*%fit$coef[11:12]
pval6 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(c(pval1,pval1,pval2,pval2,pval3,pval3,pval4,pval4,pval5,pval5,pval6,pval6,round(x[,4],4)[13:15]),4))
qval <- rep(round(p.adjust(c(pval1,pval2,pval3,pval4,pval5,pval6),method="BH"),4),each=2)
x <- x[1:12,]
x <- cbind(x,qval)

temp <- t(fit$coef[1:12])%*%solve(fit$var[c(1:12),c(1:12)])%*%fit$coef[1:12]
pval.global <- 1 - pchisq(temp, 12)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidmvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCmvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[7,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[8,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbmvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[9,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[10,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[11,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog5tert.tex",ncolumns=6,sep=" & ")
write(x[12,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSmvtog5tert.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog5tert.tex",ncolumns=1)

kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbMtogtert) &
                           !is.na(V2primtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) + factor(Avidprimtert[kp]) + factor(ADCCprimtert[kp]) + factor(NAbMtogtert[kp]) 
                + factor(V2primtert[kp]) + factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[5:6])%*%solve(fit$var[c(5:6),c(5:6)])%*%fit$coef[5:6]
pval3 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[7:8])%*%solve(fit$var[c(7:8),c(7:8)])%*%fit$coef[7:8]
pval4 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[9:10])%*%solve(fit$var[c(9:10),c(9:10)])%*%fit$coef[9:10]
pval5 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[11:12])%*%solve(fit$var[c(11:12),c(11:12)])%*%fit$coef[11:12]
pval6 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(c(pval1,pval1,pval2,pval2,pval3,pval3,pval4,pval4,pval5,pval5,pval6,pval6,round(x[,4],4)[13:15]),4))
qval <- rep(round(p.adjust(c(pval1,pval2,pval3,pval4,pval5,pval6),method="BH"),4),each=2)
x <- x[1:12,]
x <- cbind(x,qval)

temp <- t(fit$coef[1:12])%*%solve(fit$var[c(1:12),c(1:12)])%*%fit$coef[1:12]
pval.global <- 1 - pchisq(temp, 12)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidmvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCmvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[7,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[8,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbmvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[9,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[10,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[11,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog6tert.tex",ncolumns=6,sep=" & ")
write(x[12,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSmvtog6tert.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog6tert.tex",ncolumns=1)

kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbprimtert) &
                           !is.na(V2togtert) & !is.na(CD4ICSprimtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) + factor(Avidprimtert[kp]) + factor(ADCCprimtert[kp]) + factor(NAbprimtert[kp]) 
                + factor(V2togtert[kp]) + factor(CD4ICSprimtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[5:6])%*%solve(fit$var[c(5:6),c(5:6)])%*%fit$coef[5:6]
pval3 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[7:8])%*%solve(fit$var[c(7:8),c(7:8)])%*%fit$coef[7:8]
pval4 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[9:10])%*%solve(fit$var[c(9:10),c(9:10)])%*%fit$coef[9:10]
pval5 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[11:12])%*%solve(fit$var[c(11:12),c(11:12)])%*%fit$coef[11:12]
pval6 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(c(pval1,pval1,pval2,pval2,pval3,pval3,pval4,pval4,pval5,pval5,pval6,pval6,round(x[,4],4)[13:15]),4))
qval <- rep(round(p.adjust(c(pval1,pval2,pval3,pval4,pval5,pval6),method="BH"),4),each=2)
x <- x[1:12,]
x <- cbind(x,qval)

temp <- t(fit$coef[1:12])%*%solve(fit$var[c(1:12),c(1:12)])%*%fit$coef[1:12]
pval.global <- 1 - pchisq(temp, 12)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidmvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCmvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[7,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[8,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbmvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[9,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[10,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[11,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog7tert.tex",ncolumns=6,sep=" & ")
write(x[12,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSmvtog7tert.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog7tert.tex",ncolumns=1)


kp <- dat$trt=="VACCINE" & !is.na(IgAprimtert) & !is.na(Avidprimtert) & !is.na(ADCCprimtert) & !is.na(NAbprimtert) &
                           !is.na(V2primtert) & !is.na(CD4lumtogtert) 
fit <- cch(Surv(flrtime[kp],flrstatus[kp]) ~ factor(IgAprimtert[kp]) + factor(Avidprimtert[kp]) + factor(ADCCprimtert[kp]) + factor(NAbprimtert[kp]) 
                + factor(V2primtert[kp]) + factor(CD4lumtogtert[kp]) + sex[kp] + risk.medium[kp] + risk.high[kp],
                subcoh=in.subcohort[kp],stratum=stratuminds[kp],
                id=dat$pin[kp],cohort.size=cohortstratasizes,method="II.Borgan")
x <- summary(fit)$coeff
temp <- t(fit$coef[1:2])%*%solve(fit$var[c(1:2),c(1:2)])%*%fit$coef[1:2]
pval1 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[3:4])%*%solve(fit$var[c(3:4),c(3:4)])%*%fit$coef[3:4]
pval2 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[5:6])%*%solve(fit$var[c(5:6),c(5:6)])%*%fit$coef[5:6]
pval3 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[7:8])%*%solve(fit$var[c(7:8),c(7:8)])%*%fit$coef[7:8]
pval4 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[9:10])%*%solve(fit$var[c(9:10),c(9:10)])%*%fit$coef[9:10]
pval5 <- 1 - pchisq(temp, 2)
temp <- t(fit$coef[11:12])%*%solve(fit$var[c(11:12),c(11:12)])%*%fit$coef[11:12]
pval6 <- 1 - pchisq(temp, 2)
x <- cbind(round(exp(x[,1]),2),round(exp(x[,1] - x[,2]*1.96),2),
round(exp(x[,1] + x[,2]*1.96),2),round(x[,4],4),round(c(pval1,pval1,pval2,pval2,pval3,pval3,pval4,pval4,pval5,pval5,pval6,pval6,round(x[,4],4)[13:15]),4))
qval <- rep(round(p.adjust(c(pval1,pval2,pval3,pval4,pval5,pval6),method="BH"),4),each=2)
x <- x[1:12,]
x <- cbind(x,qval)

temp <- t(fit$coef[1:12])%*%solve(fit$var[c(1:12),c(1:12)])%*%fit$coef[1:12]
pval.global <- 1 - pchisq(temp, 12)

write(x[1,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1IgAmvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[2,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2IgAmvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[3,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1Avidmvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[4,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2Avidmvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[5,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1ADCCmvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[6,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2ADCCmvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[7,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1NAbmvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[8,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2NAbmvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[9,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1V2mvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[10,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2V2mvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[11,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table1CD4ICSmvtog8tert.tex",ncolumns=6,sep=" & ")
write(x[12,],file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/table2CD4ICSmvtog8tert.tex",ncolumns=6,sep=" & ")
write(round(pval.global,4),file="/trials/vaccine/thai_trial_rv144/case_control_analysis/tables/pvalmvtog8tert.tex",ncolumns=1)

################################
# Plots
################################

##### Need to weight non-cases appropriately for cum inc plot:
## Vaccine recipients only

kp <- dat$trt=="VACCINE" & !is.na(IgAprim)
wt <- ifelse(flrstatus[kp]==1,1,ifelse(flrstatus[kp]==0 & dat$dem_sex[kp]=="Female" & dat$perprot[kp]=="Yes" & dat$vaccno[kp]==4,
1/alpha1.0,
ifelse(flrstatus[kp]==0 & dat$dem_sex[kp]=="Male" & dat$perprot[kp]=="Yes" & dat$vaccno[kp]==4,
1/alpha2.0,
ifelse(flrstatus[kp]==0 & dat$dem_sex[kp]=="Female" & dat$perprot[kp]=="No" & dat$vaccno[kp]==4,
1/alpha3.0,
ifelse(flrstatus[kp]==0 & dat$dem_sex[kp]=="Male" & dat$perprot[kp]=="No" & dat$vaccno[kp]==4,
1/alpha4.0,1/alpha5.0)))))

##############################################################
# KM plots by antibody quartiles, weighted KM estimates  

scl <- 365.25/12
postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincIgAtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~IgAprimtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7) 
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Variable IgA Binding M-B gD-")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincAvidtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~Avidprimtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Variable Avidity A244 gD- D11")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincADCCtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~ADCCprimtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Variable ADCC luc 92TH023 gD-")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincNAbtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~NAbprimtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Variable NAb TZM-bl/A3R5 AUC-MB")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincV2tert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~V2primtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Variable V2 gp70 V1-V2")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincCD4ICStert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~CD4ICSprimtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Variable CD4+ ICS 92TH023")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincIgAtogtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~IgAtogtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Toggle Variable IgA A244 gD- D11")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincAvidtogtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~Avidtogtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Toggle Variable Avidity MN gD-")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincADCCtogtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~ADCCtogtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Toggle Variable ADCC gp120 A244 gD-")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincNAbAtogtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~NAbAtogtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Toggle Variable NAb TZM-bl AUC-MB-3")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincNAbStogtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~NAbStogtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Toggle Variable NAb TZM-bl-1")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincNAbMtogtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~NAbMtogtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Toggle Variable NAb A3R5 AUC-MB-2")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincV2togtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~V2togtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Toggle Variable V2 cyclic 42aa")
dev.off()

postscript("/trials/vaccine/thai_trial_rv144/case_control_analysis/figures/figure1cumincCD4lumtogtert.ps")
par(mfrow=c(1,1),las=1,cex.axis=1.5,cex.lab=1.5,cex.main=1.3,oma=c(0,3,0,0))
KM <- survfit(Surv(flrtime[kp]/scl, flrstatus[kp])~CD4lumtogtert[kp], type="kaplan-meier",weights=wt)
plot(KM,fun="event",col=c("black","red","blue"),lty=c(1,2,3),xlab="Months since the week 26 visit",
ylab="",lwd=4.7)
box()
mtext("Probability of acquiring HIV infection",side=2,cex=1.4,outer=T,las=3,line=1)
legend(x="topleft",legend=c("Low","Medium","High"),col=c("black","red","blue"),lty=c(1,2,3),cex=1.5,lwd=4.7)
title("Cumulative HIV Incidence by Category of Primary Toggle Variable CD4+ Luminex cytokines")
dev.off()
}