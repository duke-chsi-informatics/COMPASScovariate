
####################
### Revised 9/1/11 to have both tps and cch versions calculate sampling weights with reference to *case-cohort eligible* study population


library(osDesign)
library(survival)
library(xtable)
# Read in the master RV144 data-set:
dat1 <- read.csv("~/CHSI/PackageDev/COMPASScovariate/incubator/Frances/RV144/Regression/rv144_master_wk26.csv",header=T)

# Windows version
#dat1 <- read.csv("./Analysis/rv144_master_wk26.csv",header=T)

n1 <- nrow(dat1)

# Read in the immune response data
dat2 <- read.csv("~/CHSI/PackageDev/COMPASScovariate/incubator/Frances/RV144/Regression/rv144_data_wk26_correlates_primary_scaled_single_imp.csv",header=T)

#Windows version
#dat2 <- read.csv("T:/vaccine/thai_trial_rv144/case_control_analysis/pdata/compiled/rv144_data_wk26_correlates_primary_scaled_single_imp.csv",header=T)

n2 <- nrow(dat2)

# Merge the data-sets:
dat <- merge(dat1,dat2,by.y="pin")

# Set the flrtime as the time since the Week 26 visit
flrtime <- dat$time_to_infect - dat$time_to_wk26

## Remove the infected subjects with a negative failure time, as they likely were infected
## before the immune response was measured
#remv <- c(1:nrow(dat))[flrtime <= 0]
#dat <- dat[-remv,]
n <- nrow(dat)

# Make short names for convenience:
IgAprim  <- dat$primary_iga_score_tomaras_wk26
Avidprim <- dat$primary_avidity_a244_gdneg_delta11_alam_wk26
ADCCprim <- dat$primary_adcc_luc_92th023_evans_wk26
NAbprim <- dat$primary_nab_tzmbl_a3r5_auc_afrims_siriraj_montefiori_wk26
V2prim <- dat$primary_v2_gp70_v1v2_zollapazner_wk26
CD4ICSprim  <- dat$primary_cd4_ics_any_mcelrath_wk26 # this contains imputed values from single imputation
CD4ICSprim_imp <- dat$primary_cd4_ics_any_mcelrath_wk26_imp
IgAtog <- dat$toggle_iga_gdneg_delta11_tomaras_wk26
Avidtog <- dat$toggle_avidity_mn_gdneg_alam_wk26
ADCCtog <- dat$toggle_adcc_gp120coated_A244_ferrari_wk26
NAbAtog <- dat$toggle_nab_tzmbl_aucmb_afrims_wk26
NAbStog <- dat$toggle_nab_tzmbl_aucmb_siriraj_wk26
NAbMtog <- dat$toggle_nab_tzmbl_a3r5_montefiori_wk26
V2tog <- dat$toggle_v2_cyclic_peptide_42aa_wk26
CD4lumtog_imp <- dat$toggle_luminex_cytokines_mcelrath_wk26_imp
CD4lumtog <- dat$toggle_luminex_cytokines_mcelrath_wk26
flrstatus <- ifelse(dat$infect=="Yes",1,0)
sex <- ifelse(dat$dem_sex=="Female",1,0)
risk.cat <- ifelse(dat$BRA_risk=="Low",0,ifelse(dat$BRA_risk=="Medium",1,2))
risk.medium <- ifelse(dat$BRA_risk=="Medium",1,0)
risk.high <- ifelse(dat$BRA_risk=="High",1,0)

# Load categorical versions
IgAprimtert  <- dat$primary_iga_score_tomaras_wk26_cat
Avidprimtert <- dat$primary_avidity_a244_gdneg_delta11_alam_wk26_cat
ADCCprimtert <- dat$primary_adcc_luc_92th023_evans_wk26_cat
NAbprimtert <- dat$primary_nab_tzmbl_a3r5_auc_afrims_siriraj_montefiori_wk26_cat
V2primtert <- dat$primary_v2_gp70_v1v2_zollapazner_wk26_cat
CD4ICSprimtert  <- dat$primary_cd4_ics_any_mcelrath_wk26_cat
IgAtogtert <- dat$toggle_iga_gdneg_delta11_tomaras_wk26_cat
Avidtogtert <- dat$toggle_avidity_mn_gdneg_alam_wk26_cat
ADCCtogtert <- dat$toggle_adcc_gp120coated_A244_ferrari_wk26_cat
NAbAtogtert <- dat$toggle_nab_tzmbl_aucmb_afrims_wk26_cat
NAbStogtert <- dat$toggle_nab_tzmbl_aucmb_siriraj_wk26_cat
NAbMtogtert <- dat$toggle_nab_tzmbl_a3r5_montefiori_wk26_cat
V2togtert <- dat$toggle_v2_cyclic_peptide_42aa_wk26_cat
CD4lumtogtert <- dat$toggle_luminex_cytokines_mcelrath_wk26_cat


options(contrasts=c(factor="contr.treatment",ordered="contr.poly"))


################################
# Create dichotomized variables:

IgAprimdichot <- ifelse(IgAprimtert==3,1,0)
Avidprimdichot <- ifelse(Avidprimtert==3,1,0)
ADCCprimdichot <- ifelse(ADCCprimtert==3,1,0)
NAbprimdichot <- ifelse(NAbprimtert==3,1,0)
V2primdichot <- ifelse(V2primtert==3,1,0)
CD4ICSprimdichot <- ifelse(CD4ICSprimtert==3,1,0)
IgAtogdichot <- ifelse(IgAtogtert==3,1,0)
Avidtogdichot <- ifelse(Avidtogtert==3,1,0)
ADCCtogdichot <- ifelse(ADCCtogtert==3,1,0)
NAbAtogdichot <- ifelse(NAbAtogtert==3,1,0)
NAbStogdichot <- ifelse(NAbStogtert==3,1,0)
NAbMtogdichot <- ifelse(NAbMtogtert==3,1,0)
V2togdichot <- ifelse(V2togtert==3,1,0)
CD4lumtogdichot <- ifelse(CD4lumtogtert==3,1,0)


##############################
# two-phase design inference / estimation set up - added by YH

# Define phase I strata, excluded the empty 6th stratum: Male/no/<4
# nn0: vector indicating the number of controls for each of the 5 Phase I strata
nn0 <- c(length(dat1$pin[dat1$infect=="No" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Female" & dat1$perprot=="Yes" & dat1$vaccno==4]),
    length(dat1$pin[dat1$infect=="No" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Male" & dat1$perprot=="Yes" & dat1$vaccno==4]),
    length(dat1$pin[dat1$infect=="No" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Female" & dat1$perprot=="No" & dat1$vaccno==4]),
    length(dat1$pin[dat1$infect=="No" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Male" & dat1$perprot=="No" & dat1$vaccno==4]),
    length(dat1$pin[dat1$infect=="No" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Female" & dat1$perprot=="No" & dat1$vaccno<4]))

# nn1: vector indicating the number of cases for each of the 5 Phase I strata
nn1 <- c(length(dat1$pin[dat1$infect=="Yes" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Female" & dat1$perprot=="Yes" & dat1$vaccno==4]),
    length(dat1$pin[dat1$infect=="Yes" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Male" & dat1$perprot=="Yes" & dat1$vaccno==4]),
    length(dat1$pin[dat1$infect=="Yes" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Female" & dat1$perprot=="No" & dat1$vaccno==4]),
    length(dat1$pin[dat1$infect=="Yes" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Male" & dat1$perprot=="No" & dat1$vaccno==4]),
    length(dat1$pin[dat1$infect=="Yes" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Female" & dat1$perprot=="No" & dat1$vaccno<4]))
############################# -- added by YH


#############################
# two-phase design inference collapsing over per-protocol status and number of vaccinations strata

# Define phase I strata
# nn0.collapsed: vector indicating the number of controls for each of the 2 Phase I strata
nn0.collapsed <- c(length(dat1$pin[dat1$infect=="No" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Female"]),
    length(dat1$pin[dat1$infect=="No" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Male"]))

# nn1.collapsed: vector indicating the number of cases for each of the 2 Phase I strata
nn1.collapsed <- c(length(dat1$pin[dat1$infect=="Yes" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Female"]),
    length(dat1$pin[dat1$infect=="Yes" & dat1$cc_cohort==1 & dat1$trt=="VACCINE" & dat1$dem_sex=="Male"]))
#############################



##############################
# case-cohort inference / estimation set up

# Define inverse sampling weights
# Note all cases have weight 1

alpha1 <- length(unique(dat$pin[dat$trt=="VACCINE" & dat$dem_sex=="Female" & dat$perprot=="Yes" & dat$vaccno==4]))/
                length(dat1$pin[dat1$trt=="VACCINE" & dat1$cc_cohort==1 & dat1$dem_sex=="Female" & dat1$perprot=="Yes" & dat1$vaccno==4])

alpha2 <- length(unique(dat$pin[dat$trt=="VACCINE" & dat$dem_sex=="Male" & dat$perprot=="Yes" & dat$vaccno==4]))/
                length(dat1$pin[dat1$trt=="VACCINE" & dat1$cc_cohort==1 & dat1$dem_sex=="Male" & dat1$perprot=="Yes" & dat1$vaccno==4])

alpha3 <- length(unique(dat$pin[dat$trt=="VACCINE" & dat$dem_sex=="Female" & dat$perprot=="No" & dat$vaccno==4]))/
                length(dat1$pin[dat1$trt=="VACCINE" & dat1$cc_cohort==1 & dat1$dem_sex=="Female" & dat1$perprot=="No" & dat1$vaccno==4])

alpha4 <- length(unique(dat$pin[dat$trt=="VACCINE" & dat$dem_sex=="Male" & dat$perprot=="No" & dat$vaccno==4]))/
                length(dat1$pin[dat1$trt=="VACCINE" & dat1$cc_cohort==1 & dat1$dem_sex=="Male" & dat1$perprot=="No" & dat1$vaccno==4])

alpha5 <- length(unique(dat$pin[dat$trt=="VACCINE" & dat$dem_sex=="Female" & dat$perprot=="No" & dat$vaccno<4]))/
                length(dat1$pin[dat1$trt=="VACCINE" & dat1$cc_cohort==1 & dat1$dem_sex=="Female" & dat1$perprot=="No" & dat1$vaccno<4])

alpha6 <- length(unique(dat$pin[dat$trt=="VACCINE" & dat$dem_sex=="Male" & dat$perprot=="No" & dat$vaccno<4]))/
                length(dat1$pin[dat1$trt=="VACCINE" & dat1$cc_cohort==1 & dat1$dem_sex=="Male" & dat1$perprot=="No" & dat1$vaccno<4])

# alpha6 is zero- so eliminate that category

stratuminds <- ifelse(dat$dem_sex=="Female" & dat$perprot=="Yes" & dat$vaccno==4,1,
               ifelse(dat$dem_sex=="Male" & dat$perprot=="Yes" & dat$vaccno==4,2,
               ifelse(dat$dem_sex=="Female" & dat$perprot=="No" & dat$vaccno==4,3,
               ifelse(dat$dem_sex=="Male" & dat$perprot=="No" & dat$vaccno==4,4,5))))

cohortstratasizes <- floor(table(stratuminds[dat$trt=="VACCINE"])*c(1/alpha1,1/alpha2,1/alpha3,1/alpha4,1/alpha5))
in.subcohort <- dat$infect=="No"
#############################


#############################
# case-cohort inference/estimation collapsed over per-protocol and number of vaccinations strata

# Define inverse sampling weights
# Note all cases have weight 1

alpha1.collapsed <- length(unique(dat$pin[dat$trt=="VACCINE" & dat$dem_sex=="Female"]))/
                length(dat1$pin[dat1$trt=="VACCINE" & dat1$cc_cohort==1 & dat1$dem_sex=="Female"])

alpha2.collapsed <- length(unique(dat$pin[dat$trt=="VACCINE" & dat$dem_sex=="Male"]))/
                length(dat1$pin[dat1$trt=="VACCINE" & dat1$cc_cohort==1 & dat1$dem_sex=="Male"])

stratuminds.collapsed <- ifelse(dat$dem_sex=="Female",1,2)

cohortstratasizes.collapsed <- floor(table(stratuminds.collapsed[dat$trt=="VACCINE"])*c(1/alpha1,1/alpha2))
#############################


#source("/trials/vaccine/thai_trial_rv144/case_control_analysis/code/primaryinferences_tps.r")

#source("/trials/vaccine/thai_trial_rv144/case_control_analysis/code/primaryinferences_cch.r")

#source("/trials/vaccine/thai_trial_rv144/case_control_analysis/code/modelselection_main.r")

#source("/trials/vaccine/thai_trial_rv144/case_control_analysis/code/interactionplots.r")

#source("/trials/vaccine/thai_trial_rv144/case_control_analysis/code/interactionmodels.r")

#source("/trials/vaccine/thai_trial_rv144/case_control_analysis/code/moreinteractionplots.r")

