library(Kmisc)
library(COMPASS)
library(ggplot2)
library(gridExtra)
library(data.table.extras)
library(osDesign)
library(pROC)
library(plyr)
library(reshape2)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
source("Response/compute_PI.R") # ignore
source("computePI_bc.R") # ignore
fit<-readRDS("./manuscript/test/RV144_CC.rds") # don't need this bc I have my own COMPASS fit
q =  c(1, 1.2, 2)
logPI <- computePI(q , fit$data$n_s, fit$fit$categories)
PI <- exp(logPI)
PI_corrected <- computePI_bc(q, fit$data$n_s, fit$data$n_u, fit$fit$categories)

colnames(PI) = q
colnames(PI_corrected) = q


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# TODO: get these below from Box folder; need to calculate p-values
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
source("./Correlates/my.tps.function.R")
source("Correlates/primaryinferences_cch.r")
source("Correlates/primaryinferences_main_090111.r")

assign("depthtrigger", 100, data.table:::.global)

## Merge in infection status post-hoc (which is case and which is control)
## has case-control
orig <- read.csv("./Correlates/PTID_map_case_control.csv", header=T, sep=",")
# 
# to_merge <- dat[ c("pin", "infect")]
# to_merge <- merge(to_merge, orig, by.x="pin", by.y="PTID_orig", all.y=TRUE, sort = FALSE)
# to_merge <- to_merge[ !(names(to_merge) %in% c("X", "pin")) ]
# to_merge <- unique(to_merge)
# 
# meta <- fit$data$meta
# meta <- merge(meta, to_merge, by="PTID", all.x=TRUE, sort = FALSE)
# meta$infect <- as.character(meta$infect)
# meta$infect[ is.na(meta$infect) ] <- "PLACEBO"
# meta$infect[which(meta$infect == "No")] <- "NON-INFECTED"
# meta$infect[which(meta$infect == "Yes")] <- "INFECTED"
# meta$infect <- factor(meta$infect, levels=c("NON-INFECTED", "INFECTED", "PLACEBO"))
# 
# fit$data$meta <- meta

I <- dim(fit$data$n_s)[1] #subjects
K <- dim(fit$data$n_s)[2] #categories
T <- dim(fit$fit$alpha_s)[1] #
K1 <- K-1
M <- dim(fit$data$categories)[2]-1
Cnames <- colnames(fit$data$categories)[1:M]
## get the vaccination and infection status for individuals in the n_s/n_u
## note: all individuals in n_s are in the metadata
indiv <- rownames(fit$data$n_s) #individual PTIDs
meta <- fit$data$meta
vacc <- unique(meta[ meta$PTID %in% indiv, c("PTID", "vaccine")])
vacc <- vacc[ match(indiv, vacc$PTID), ]
rownames(vacc) <- 1:nrow(vacc)

infect <- unique(meta[ meta$PTID %in% indiv, c("PTID", "infect")])
infect <- infect[ match(indiv, infect$PTID), ]
rownames(infect) <- 1:nrow(infect)


#######################
## get explanatory variables for regression
rPTID =  unique(orig[,2:3])
colnames(rPTID) = c("PTID", "PTID_orig")

## meta$PTID <- swap( meta$PTID, rPTID$PTID_orig, rPTID$PTID )
sel_outcome <- NULL #reorder n_s/n_u according to dat$pin
sel_c <- NULL # which dat$pin maps with indiv
PTID_sub <- NULL

for (ii in 1:length(dat$pin)) {
  tpp = which(indiv == rPTID[which(rPTID[,2]==dat$pin[ii]),1])
  if(length(tpp) >1){
    cat(ii)
  }
  if (length(tpp)>0) {
    sel_c = c(sel_c, ii)
    sel_outcome = c(sel_outcome, tpp)
    PTID_sub<- c(PTID_sub, rPTID[rPTID[,1] == indiv[tpp],2])
  }
}

## already have this information: Lynn is merging flow with metadata here
flrstatus = flrstatus[sel_c]
IgAprim = IgAprim[sel_c]
V2prim = V2prim[sel_c]
NAbprim = NAbprim[sel_c]

Avidprim <- Avidprim[sel_c]
ADCCprim <- ADCCprim[sel_c]

CD4ICSprim  <- CD4ICSprim[sel_c] 

sex = sex[sel_c]
risk.medium = risk.medium[sel_c]
risk.high = risk.high[sel_c]
stratuminds = stratuminds[sel_c]

PFS <- PolyfunctionalityScore(fit)
FS <- FunctionalityScore(fit)

fit23 = readRDS("./model-fit/COMPASS-fit-RV144_23.rds") ## Lynn's original COMPASS fit, replace with new one

##### original 6 primary variables ###
pvalues <- NULL
pnames <- NULL
OR <- NULL
#outcome = scale(as.vector(fit$fit$mean_gamma[sel_outcome,23]))
outcome = scale(as.vector(PFS[sel_outcome])) # can replace this with PFS, FS, mean gamma as needed
fit.tps <- try(tps(flrstatus ~ outcome + IgAprim +V2prim+NAbprim +Avidprim+ADCCprim+ 
                     sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE) 
for(ii in 2:7){
  x <- as.matrix(cbind(fit.tps$coef[ii], round(exp(fit.tps$coef[ii]),3),round(exp(fit.tps$coef[ii] - sqrt(fit.tps$covm[ii,ii])*1.96),3),
                       round(exp(fit.tps$coef[ii] + sqrt(fit.tps$covm[ii,ii])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[ii]/sqrt(fit.tps$covm[ii,ii])))),1.0),4)))
  colnames(x)<-c("Coef","OR","CI.low","CI.up","p-value")
  pvalues <- c(pvalues, x[5])
  OR <- c(OR, x[2])
  pnames <- c(pnames, names(fit.tps$coef)[ii])
  
}

names(pvalues) <- pnames
qvalues <- p.adjust(pvalues, method = "fdr", n = length(pvalues))

#####################################
#########volcano plot ###############
#####################################
Cate <- NULL
Cnames[which(Cnames == 'CD154')] = "CD40L"
dd <- fit$fit$categories[1:23,1:6]
colnames(dd)[which(colnames(dd) == "CD154")] = "CD40L"
for(k in 1:23){
  tmp = Cnames[which(dd[k,] == 1)]
  Cate <- c(Cate, paste(tmp, "+", sep="", collapse = ""))
}

#log fold
#logfoldchange <- log(fit$data$n_s[,1:23] + 10^(-20)) - log(fit$data$counts_s)-
#                 log(fit$data$n_u[,1:23] + 10^(-20)) + log(fit$data$counts_u)
# difference
logfoldchange <- fit$data$n_s[,1:23]/fit$data$counts_s - 
                 fit$data$n_u[,1:23]/fit$data$counts_u

colnames(logfoldchange) = Cate
DFF.fold <- melt(logfoldchange)
colnames(DFF.fold)[3] = "logFC"

PrResponse <- fit$fit$mean_gamma[,1:23]
colnames(PrResponse) = Cate
DFF.Pr <- melt(PrResponse)
colnames(DFF.Pr)[3] = "PrResponse"

DFF <- cbind(DFF.fold, DFF.Pr$PrResponse)
colnames(DFF)[4] <- colnames(DFF.Pr)[3]

ggplot(data=DFF, aes(x=logFC, y=PrResponse)) +
  geom_point(alpha=0.4, size=1.75) + facet_wrap(~Var2, scales = "free")+
  xlab("ps-pu") + ylab("Prob. Response")+
  theme(axis.text.x = element_text(angle=90))

# "IFNg+CD40L+"
selcat <- which(Cate %in% c("TNFa+IL2+CD40L+","IFNg+CD40L+" ))

#selcat <- 7:23
logfoldchange <- fit$data$n_s[,selcat]/fit$data$counts_s - 
  fit$data$n_u[,selcat]/fit$data$counts_u
colnames(logfoldchange) <- Cate[selcat]
DFF.fold <- melt(logfoldchange); colnames(DFF.fold)[3] = "ps-pu"
PrResponse <- fit$fit$mean_gamma[,selcat]
colnames(PrResponse) = Cate[selcat]
DFF.Pr <- melt(PrResponse)
colnames(DFF.Pr)[3] = "PrResponse"
DFF1 <- merge(DFF.fold, DFF.Pr, by = c("Var1", "Var2"))

selns <- fit$data$n_s[,selcat]
colnames(selns)<- Cate[selcat]
DFF.ns <- melt(selns); colnames(DFF.ns)[3] <- "ns"
selnu <- fit$data$n_u[,selcat]
colnames(selnu) <- Cate[selcat]
DFF.nu <- melt(selnu); colnames(DFF.nu)[3] <- "nu"
DFF2 <- merge(DFF.ns, DFF.nu, by = c("Var1", "Var2"))

DFF <- merge(DFF1, DFF2, by = c("Var1", "Var2"))
DFF <- cbind(DFF, "sqrt(ns)" = sqrt(DFF$ns), "sqrt(nu)" = sqrt(DFF$nu))
#####cts color scheme
ggplot(data=DFF, aes(x=`ps-pu`, y=PrResponse, color = `sqrt(ns+nu)`)) +
  geom_point(alpha=0.4, size=1.75) + facet_wrap(~Var2, scales = "free")+
  scale_gradient()+
  #scale_colour_gradient2( mid = "grey", high = "purple")+
  xlab("ps-pu") + ylab("Prob. Response") + theme_bw()+
  theme(axis.text.x = element_text(angle=90)) 

#### color manual make nu as factors
DFF_posi <- subset(DFF, `ps-pu` > 0)
DFF1 <- cbind(DFF_posi, "col" = DFF_posi$nu)
DFF1$col <- as.factor(DFF1$col)
labels <- names(table(DFF_posi$nu))
colfunc <- colorRampPalette(c("blue","yellow","red"))
breaks <- colfunc(length(labels))[1:length(labels)]
png(file="Response/Supplementary/figs/volcano.png",res = 300, height = 1600, width = 2200)
ggplot(DFF1) + 
  geom_point(aes(x=`ps-pu`, y=PrResponse,size=nu), alpha=1) +
  scale_size_continuous(trans="sqrt",guide="legend") +
  ylab("Response probability") + scale_x_continuous(limits = c(0, 0.00078))+
  facet_wrap(~Var2, scales = "free") +theme_bw()+
  theme(axis.text.x = element_text(angle=90))

####### discretize ############
################################################
DFF_posi <- subset(DFF, `ps-pu` > 0)
bincol <- function(x,low,medium,high) {
  breaks <- function(x) pretty(range(x), n = 1*nclass.Sturges(x), min.n = 1)  
  colfunc <- colorRampPalette(c(low, medium, high))
  binned <- cut(x,breaks(x), include.lowest=TRUE)  
  res <- colfunc(length(unique(as.integer(binned))))[as.integer(binned)]
  names(res) <- as.character(binned)
  res
}

var_plot = DFF_posi$`sqrt(nu)`
labels <- unique(names(bincol(var_plot,"blue","yellow","red")))
breaks <- unique(bincol(var_plot,"blue","yellow","red"))

values <- unlist(lapply(labels, function(x) {
  tmp <- strsplit((strsplit(unlist(strsplit(x, ",")), " "))[[1]], "")[[1]]
  paste0(unlist(tmp[2:length(tmp)]), collapse = "")
}))

o <- order(as.numeric(values))
breaks <- breaks[o]
labels <- labels[o]

ggplot(DFF_posi) + 
  geom_point(aes(x=`ps-pu`, y=PrResponse,
                 colour=bincol(`sqrt(nu)`,"blue","yellow","red"),size=`sqrt(nu)`), alpha=1) +
  scale_color_identity("sqrt(nu)", labels=labels, 
                       breaks=breaks, guide="legend") +
  facet_wrap(~Var2, scales = "free") +theme_bw()+
  theme(axis.text.x = element_text(angle=90))

################### ns + nu ########
labels <- unique(names(bincol(DFF$`ns+nu`,"blue","yellow","red")))
breaks <- unique(bincol(DFF$`ns+nu`,"blue","yellow","red"))

values <- unlist(lapply(labels, function(x) {
  tmp <- strsplit((strsplit(unlist(strsplit(x, ",")), " "))[[1]], "")[[1]]
  paste0(unlist(tmp[2:length(tmp)]), collapse = "")
}))

o <- order(as.numeric(values))
breaks <- breaks[o]
labels <- labels[o]

ggplot(DFF) + 
  geom_point(aes(x=`ps-pu`, y=PrResponse,
                 colour=bincol(`ns+nu`,"blue","yellow","red")), alpha=0.6, size=1.85) +
  scale_color_identity("ns+nu", labels=labels, 
                       breaks=breaks, guide="legend") +
  facet_wrap(~Var2, scales = "free") +theme_bw()+
  theme(axis.text.x = element_text(angle=90))

######################################
##### boxplot for PI, FS, and PFS#####
######################################
DFF <- NULL
diff <- rbind(FS, vacc[,2], "Functionality score")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)
diff <- rbind(PFS, vacc[,2],"Polyfunctionality score")

rownames(diff) <- NULL
DFF <- cbind(DFF, diff)
diff <- rbind(PI[,"1"], vacc[,2], "PI")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)

diff <- rbind(PI_corrected[,"1"], vacc[,2], "PI_corrected")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)


diff <- rbind(PI[,"1"], vacc[,2], "PI, q = 1")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)

diff <- rbind(PI[,"1.2"], vacc[,2], "PI, q = 1.2")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)

diff <- rbind(PI[,"2"], vacc[,2], "PI, q = 2")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)

diff <- rbind(PI[,"2"], vacc[,2], "PI, q = 2")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)

diff<- rbind(FS_MIMOSA[match(names(FS), names(FS_MIMOSA))], vacc[,2], "FS_MIMOSA")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)

diff<- rbind(PFS_MIMOSA[match(names(PFS), names(PFS_MIMOSA))], vacc[,2], "PFS_MIMOSA")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)

diff<- rbind(UMIMOSA$`TNFa&IFNg&!IL4&IL2&CD154&!IL17a`[match(names(FS), UMIMOSA$PTID)], vacc[,2],"5degree_MIMOSA")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)

diff<- rbind(fit23$mean_gamma[,1], vacc[,2], "5degree_COMPASSuni")
rownames(diff) <- NULL
DFF<- cbind(DFF, diff)

diff<- rbind(fit$fit$mean_gamma[,23], vacc[,2], "5degree_COMPASS")
rownames(diff) <- NULL
DFF <- cbind(DFF, diff)

colnames(DFF) <- NULL
DFF <- data.frame(t(DFF))
colnames(DFF) <- c("Score","Status","names")
DFF$Score <- as.numeric(as.character(DFF$Score))
levels(DFF$Status) = levels(vacc[,2])
DFF$names <- factor(DFF$names, 
                    levels=c("Functionality score", "Polyfunctionality score", "PI", "PI_corrected"))

DFF$names <- factor(DFF$names, c("PI, q = 1", "PI, q = 1.2", 
                                 "PI, q = 2","Functionality score","Polyfunctionality score" ,"FS_MIMOSA", 
                                 "PFS_MIMOSA", "5degree_MIMOSA","5degree_COMPASSuni","5degree_COMPASS"))

scale <- 300 / 72
#tiff(file="Response/Figures/boxplots_RV144.tiff", res=300, width=680*scale, height=480*scale)
#png(file="Response/Supplementary/figs/boxplots_RV144.png", res=300, width=2200, height=1500)
png(file="Response/Supplementary/figs/boxplots_RV144_v2.png", res=300, width=2200, height=2400)
ggplot(DFF, aes(Status, Score))+ geom_boxplot(outlier.size = 0)+
  geom_jitter(aes(color=Status))+facet_wrap(~names, scales = "free_y")+
  scale_colour_grey(start = 0.2, end = 0.45)+theme_bw()+ theme(legend.position="none")


subDFF <- subset(DFF, names %in% c("5degree_MIMOSA", "5degree_COMPASSuni","5degree_COMPASS"))
subDFF <- subset(DFF, names == "Polyfunctionality score")
wilcox.test(subDFF$Score[which(subDFF$Status == "VACCINE")],
            subDFF$Score[which(subDFF$Status == "PLACEBO")],alternative = "two.sided")

##############################
##### ROC (PFS, FS, PI) ##################
library(pROC)
library(reshape2)

DFF <- data.table( PFS_COMPASS = unname(PFS), FS_COMPASS = unname(FS), `PI, q = 1` = PI[,"1"], `PI, q = 1.2` = PI[,"1.2"], 
                   `PI, q = 2` = PI[,"2"], FS_MIMOSA = FS_MIMOSA[match(names(FS), names(FS_MIMOSA))],
                   PFS_MIMOSA = PFS_MIMOSA[match(names(PFS), names(PFS_MIMOSA))],
                   InfectionStatus = vacc$vaccine, PTID = vacc$PTID)

DFF <- data.table( PFS = unname(PFS), FS = unname(FS), 
                   FS_MIMOSA = unname(FS_MIMOSA[match(names(FS), names(FS_MIMOSA))]),
                   PFS_MIMOSA = unname(PFS_MIMOSA[match(names(PFS), names(PFS_MIMOSA))]),
                   InfectionStatus = vacc$vaccine, PTID = vacc$PTID)
nm <- colnames(DFF)[1:(length(colnames(DFF))-2)]
rocs <- lapply(nm, function(x) {
  roc(DFF$InfectionStatus, DFF[[x]])
})

roc_dt <- lapply(rocs, function(x) {
  data.table(
    Sensitivity=x$sensitivities,
    Specificity=x$specificities
  )
})

names(roc_dt) <- nm

rocs_dt <- rbindlistn(roc_dt)

rocs_dt <- rocs_dt[ order(Sensitivity), ]

rocs_dt[, Methods :=  .Names ]

#png("Response/Figures/ROC_RV144.png", res=300, width=1600, height=1300)
ggplot(rocs_dt, aes(x=1-Specificity, y=Sensitivity, col=Methods)) +
  geom_line() +
  geom_abline(a=0, b=1, lty='dashed')+theme_bw()
#dev.off()

########### ROC (MIMOSA, COMPASS, Fisher's) #####
MMIMOSA <- readRDS("./MIMOSA_RV144/MIMOSAVariables.rds")
test <- MMIMOSA[ match(names(PFS),MMIMOSA$PTID), ]

################ 2x2 fisher's test ############################
n_s <- fit$data$n_s
n_u <- fit$data$n_u
kk <- c(1,2,4,5)
sel = which(rowSums(fit$data$categories[,kk])>0)
tmp = 1:K
tmp = tmp[-sel]
n.stim = cbind(rowSums(n_s[,tmp]),rowSums(n_s[,sel]))
n.unstim = cbind(rowSums(n_u[,tmp]),rowSums(n_u[,sel]))

Table = array(0,dim=c(2,2));
pvalue_F = array(0,dim=c(I,1));
ppp_F2 = pvalue_F
for (i in 1:I) {
  Table[1,1] = n.unstim[i,2]
  Table[2,1] = n.stim[i,2]
  Table[1,2] = n.unstim[i,1];
  Table[2,2] = n.stim[i,1];
  resu = fisher.test(Table,alternative = "less");
  pvalue_F[i] <- resu$p.value
}
ppp_F2 = p.adjust(pvalue_F, method = "fdr", n = I)
ppp_F2 = as.vector(ppp_F2)

#######fisher's 2xK ############
DFF <- NULL
Table = array(0,dim=c(2,K));
pvalue_F = array(0,dim=c(I,1));
ppp_F = pvalue_F
for (i in 1:I) {
  Table[1,] =fit$data$n_u[i,] 
  Table[2,] = fit$data$n_s[i,]
  res1 = fisher.test(Table,simulate.p.value = TRUE);
  pvalue_F[i]=res1$p.value 
}
ppp_F = p.adjust(pvalue_F, method = "fdr", n = I)
Mgamma_list = as.vector(ppp_F);


Ps <- rep(0,I)
for (tt in 1:T) {
  Ps <- Ps +1*(rowSums(fit$fit$gamma[,1:K1,tt])>2)
}
Ps <- Ps/T  #COMPASS


result_cc <- readRDS("./model-fit/MIMOSA-fit-RV144.rds") #MIMOSA

DFF <- data.table(MIMOSA = result_cc[,2], MultivariateMIMOSA = test$PrRespMultivMIMOSA, COMPASS = unname(Ps), 
                  `Fisher's 2xK Exact Test` = Mgamma_list,`Fisher's 2x2 Exact Test` = ppp_F2,PI = PI[,"1"],
                  PI_corrected = PI_corrected[,"1"],
                  InfectionStatus = vacc$vaccine, PTID = vacc$PTID)
nm <- colnames(DFF)[1:(length(colnames(DFF))-2)]
rocs <- lapply(nm, function(x) {
  roc(DFF$InfectionStatus, DFF[[x]])
})


roc_dt <- lapply(rocs, function(x) {
  data.table(
    Sensitivity=x$sensitivities,
    Specificity=x$specificities
  )
})

names(roc_dt) <- nm

rocs_dt <- rbindlistn(roc_dt)

rocs_dt <- rocs_dt[ order(Sensitivity), ]

rocs_dt[, Methods :=  .Names ]

png("Response/Supplementary/figs/ROC_RV144_MIMOSAPI_v2.png", res=300, width=2000, height=1300)
ggplot(rocs_dt, aes(x=1-Specificity, y=Sensitivity, col=Methods)) +
  geom_line() +
  geom_abline(a=0, b=1, lty='dashed')+theme_bw()
dev.off()
##############################
###### Logistic Regression ######
#################################
p_value <- NULL
pnames<-NULL
outcome = scale(as.vector(PI_corrected[sel_outcome,"1"]))
fit.tps <- try(tps(flrstatus ~ outcome + IgAprim + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE) 
x <- as.matrix(cbind(fit.tps$coef[2], round(exp(fit.tps$coef[2]),3),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),3),
                     round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0),4)))
colnames(x)<-c("Coef","OR","CI.low","CI.up","p-value")
p_value <- c(p_value, x[,5])
pnames <- c(pnames, "PI, q = 1")

outcome = scale(as.vector(PI[sel_outcome,"1.2"]))
fit.tps <- try(tps(flrstatus ~ outcome + IgAprim + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE) 
x <- as.matrix(cbind(fit.tps$coef[2], round(exp(fit.tps$coef[2]),3),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),3),
                     round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0),4)))
colnames(x)<-c("Coef","OR","CI.low","CI.up","p-value")
p_value <- c(p_value, x[,5])
pnames <- c(pnames, "PI, q = 1.2")

outcome = scale(as.vector(PI[sel_outcome,"2"]))
fit.tps <- try(tps(flrstatus ~ outcome + IgAprim + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE) 
x <- as.matrix(cbind(fit.tps$coef[2], round(exp(fit.tps$coef[2]),3),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),3),
                     round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0),4)))
colnames(x)<-c("Coef","OR","CI.low","CI.up","p-value")
p_value <- c(p_value, x[,5])
pnames <- c(pnames, "PI, q = 2")
names(p_value) = pnames
print(p_value)

######################################
########### luminex ##################
#####################################
conc= read.csv("./Correlates/usm_mcelrath_cytokines_luminex_est.csv.csv")
################################
lana = length(unique(conc$analyte))
uana = unique(conc$analyte)
posit_s = array(0,dim=c(I,lana))
posit_u = posit_s
est_s = posit_s
est_u = posit_s

for ( kk in 1:lana) {
  visit800_s = NULL; visit800_u = NULL; 
  for (ii in 1:length(indiv)) {
    visit800_s =c(visit800_s, 
                  which(conc$visit==800 & conc$ptid == rPTID[which(rPTID[,1]==indiv[ii]),2]  & conc$stimulation=="92TH023Env" & conc$analyte == uana[kk]))
    visit800_u =c(visit800_u, 
                  which(conc$visit==800 & conc$ptid == rPTID[which(rPTID[,1]==indiv[ii]),2]  & conc$stimulation=="negctrl" & conc$analyte == uana[kk]))
  }
  posit_s[,kk] = visit800_s
  posit_u[,kk] = visit800_u
  est_s[,kk] = exp(conc$est.log.conc[posit_s[,kk]])
  est_u[,kk] = exp(conc$est.log.conc[posit_u[,kk]])
}

diff_est=((est_s)-(est_u))
colnames(diff_est)=uana
diff_est = cbind(diff_est,rowMeans(diff_est))
colnames(diff_est)[13] = "avg"

colnames(est_s) = uana 
colnames(est_u) = uana 

######Positivity call for luminex
min_conc = c(50, 75, 50, 50, 40, 5, 15, 20, 10, 100, 30, 10)
min_fold = c(3, 5, 4, 3, 3, 3, 3, 5, 3, 3, 3, 4)
max_back = c(50, 40, 40, 120, 40, 5, 20, 20, 5, 250, 20, 10)

pos = array(0, dim=c(I,lana))
for (ss in 1:lana) {
  t1 = (est_u[,ss])<=max_back[ss]
  t2 = ((est_s[,ss])-(est_u[,ss]))>=min_conc[ss]
  t3 = (est_s[,ss]/est_u[,ss])>= min_fold[ss]
  pos[which(t1==TRUE & t2==TRUE & t3==TRUE),ss] = 1;
}
pos_all = rowSums(pos)/lana # positivity indicator

#####common cytokines ####
common_cyto = intersect(Cnames,uana)
pos1 = rowSums(pos[,which(uana %in% common_cyto)])/length(common_cyto)# positivity indicator

fitcommon <- readRDS("./model-fit/model-fit-RV144_CD4_CD154_IL17a.rds")
FS_common <- FunctionalityScore(fitcommon)
PFS_common <- PolyfunctionalityScore(fitcommon)
DFF = cbind(pos1, FS_common, PFS_common, vacc[,2])
DFF = data.frame(DFF)
colnames(DFF) = c("Luminex","FS","PFS","Status")
DFF$Status = as.factor(DFF$Status)
levels(DFF$Status) = levels(vacc[,2])
p11<-ggplot(subset(DFF, Status =="VACCINE"), aes(x=FS, y=Luminex)) +
  geom_point(alpha=I(0.8)) +    
  geom_smooth(method="lm",colour = getPalette(22)[4], fill = getPalette(30)[30]) +xlab("Functionality score") + ylab("Proportion of expressed cytokines \n(Multiplex bead array)") +
  theme_bw()#+theme(axis.title = element_text(size=17), legend.text = element_text(size=12),legend.title = element_text(size=13))

p22<-ggplot(subset(DFF, Status =="VACCINE"), aes(x=PFS, y=Luminex)) +
  geom_point(alpha=I(0.8)) +    
  geom_smooth(method="lm",colour = getPalette(22)[4], fill = getPalette(30)[30]) +xlab("Polyfunctionality score") + ylab("Proportion of expressed cytokines \n(Multiplex bead array)") +
  theme_bw()#+theme(axis.title = element_text(size=17), legend.text = element_text(size=12),legend.title = element_text(size=13))

#tiff("Response/Figures/Luminx_FSPFS_sameCyto.tiff", res=300, width=1700, height=1000)
grid.arrange(p11, p22, ncol=2)
#dev.off()
# spearman coefficient
tmp <- subset(DFF,Status == "VACCINE")
cor.test(tmp$Luminex, tmp$FS, alternative = "two.sided", method = "spearm")
lm_FS = (lm(Luminex ~ FS, data=tmp))
sqrt(summary(lm_FS)$r.squared)
confint(lm_FS, "FS", level = 0.95)
sqrt(summary(lm(Luminex ~ PFS, data=tmp))$r.squared)

##########################################
########## Correlation with Luminex #####
#########################################
p_value <- NULL
pnames <- NULL
outcome = scale(as.vector(pos_all[sel_outcome]))
fit.tps <- try(tps(flrstatus ~ outcome + IgAprim + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE) 
x <- as.matrix(cbind(fit.tps$coef[2], round(exp(fit.tps$coef[2]),3),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),3),
                     round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0),4)))
colnames(x)<-c("Coef","OR","CI.low","CI.up","p-value")
p_value <- c(p_value, x[,5])
pnames <- c(pnames, "Luminex")

outcome = scale(as.vector(pos_all[sel_outcome]))
fit.tps <- try(tps(flrstatus ~ outcome + FS[sel_outcome] + IgAprim + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE) 
x <- as.matrix(cbind(fit.tps$coef[2], round(exp(fit.tps$coef[2]),3),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),3),
                     round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0),4)))
y <- as.matrix(cbind(fit.tps$coef[3], round(exp(fit.tps$coef[3]),3),round(exp(fit.tps$coef[3] - sqrt(fit.tps$covm[3,3])*1.96),3),
                     round(exp(fit.tps$coef[3] + sqrt(fit.tps$covm[3,3])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[3]/sqrt(fit.tps$covm[3,3])))),1.0),4)))
colnames(x)<-c("Coef","OR","CI.low","CI.up","p-value")
colnames(y)<-c("Coef","OR","CI.low","CI.up","p-value")

p_value <- c(p_value, x[,5])
pnames <- c(pnames, "Luminex (Luminex + FS)")

p_value <- c(p_value, y[,5])
pnames <- c(pnames, "FS (Luminex + FS)")

outcome = scale(as.vector(pos_all[sel_outcome]))
fit.tps <- try(tps(flrstatus ~ outcome + PFS[sel_outcome] + IgAprim + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE) 
x <- as.matrix(cbind(fit.tps$coef[2], round(exp(fit.tps$coef[2]),3),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),3),
                     round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0),4)))
y <- as.matrix(cbind(fit.tps$coef[3], round(exp(fit.tps$coef[3]),3),round(exp(fit.tps$coef[3] - sqrt(fit.tps$covm[3,3])*1.96),3),
                     round(exp(fit.tps$coef[3] + sqrt(fit.tps$covm[3,3])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[3]/sqrt(fit.tps$covm[3,3])))),1.0),4)))
#z <- as.matrix(cbind(fit.tps$coef[8], round(exp(fit.tps$coef[8]),3),round(exp(fit.tps$coef[8] - sqrt(fit.tps$covm[8,8])*1.96),3),
#                     round(exp(fit.tps$coef[8] + sqrt(fit.tps$covm[8,8])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[8]/sqrt(fit.tps$covm[8,8])))),1.0),4)))

colnames(x)<-c("Coef","OR","CI.low","CI.up","p-value")
colnames(y)<-c("Coef","OR","CI.low","CI.up","p-value")

p_value <- c(p_value, x[,5])
pnames <- c(pnames, "Luminex (Luminex + PFS)")

p_value <- c(p_value, y[,5])
pnames <- c(pnames, "PFS (Luminex + PFS)")

outcome = scale(as.vector(pos_all[sel_outcome]))
fit.tps <- try(tps(flrstatus ~ outcome + fit$fit$mean_gamma[sel_outcome,23] + IgAprim + sex + risk.medium + risk.high,nn0=nn0, nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE) 
x <- as.matrix(cbind(fit.tps$coef[2], round(exp(fit.tps$coef[2]),3),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),3),
                     round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0),4)))
y <- as.matrix(cbind(fit.tps$coef[3], round(exp(fit.tps$coef[3]),3),round(exp(fit.tps$coef[3] - sqrt(fit.tps$covm[3,3])*1.96),3),
                     round(exp(fit.tps$coef[3] + sqrt(fit.tps$covm[3,3])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[3]/sqrt(fit.tps$covm[3,3])))),1.0),4)))
colnames(x)<-c("Coef","OR","CI.low","CI.up","p-value")
colnames(y)<-c("Coef","OR","CI.low","CI.up","p-value")

p_value <- c(p_value, x[,5])
pnames <- c(pnames, "Luminex (Luminex + 5degree)")

p_value <- c(p_value, y[,5])
pnames <- c(pnames, "5degree (Luminex + 5degree)")

outcome = scale(as.vector(pos_all[sel_outcome]))
fit.tps <- try(tps(flrstatus ~ outcome + FS[sel_outcome] +
                     PFS[sel_outcome] + IgAprim + sex + risk.medium + risk.high,nn0=nn0, 
                   nn1=nn1,group=stratuminds,method="PL",cohort=TRUE),silent=TRUE) 
x <- as.matrix(cbind(fit.tps$coef[2], round(exp(fit.tps$coef[2]),3),round(exp(fit.tps$coef[2] - sqrt(fit.tps$covm[2,2])*1.96),3),
                     round(exp(fit.tps$coef[2] + sqrt(fit.tps$covm[2,2])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[2]/sqrt(fit.tps$covm[2,2])))),1.0),4)))
y <- as.matrix(cbind(fit.tps$coef[3], round(exp(fit.tps$coef[3]),3),round(exp(fit.tps$coef[3] - sqrt(fit.tps$covm[3,3])*1.96),3),
                     round(exp(fit.tps$coef[3] + sqrt(fit.tps$covm[3,3])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[3]/sqrt(fit.tps$covm[3,3])))),1.0),4)))
z <- as.matrix(cbind(fit.tps$coef[4], round(exp(fit.tps$coef[4]),3),round(exp(fit.tps$coef[4] - sqrt(fit.tps$covm[4,4])*1.96),3),
                     round(exp(fit.tps$coef[4] + sqrt(fit.tps$covm[4,4])*1.96),3),round(min(2*(1-pnorm(abs(fit.tps$coef[4]/sqrt(fit.tps$covm[4,4])))),1.0),4)))

colnames(x)<-c("Coef","OR","CI.low","CI.up","p-value")
colnames(y)<-c("Coef","OR","CI.low","CI.up","p-value")
colnames(z)<-c("Coef","OR","CI.low","CI.up","p-value")

p_value <- c(p_value, x[,5])
pnames <- c(pnames, "Luminex (Luminex + FS + PFS)")

p_value <- c(p_value, y[,5])
pnames <- c(pnames, "FS (Luminex + FS + PFS)")

p_value <- c(p_value, z[,5])
pnames <- c(pnames, "PFS (Luminex + FS + PFS)")
names(p_value) = pnames


