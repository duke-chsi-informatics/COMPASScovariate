
#--------------------------------------------------------
# Read in multiple COMPASS seeded runs
#--------------------------------------------------------
readInObjs <- function(data_dir, pattern) {
  files <- list.files(data_dir, pattern=pattern)
  seeds <- str_extract(files, "seed(.*)\\.rds",group=1)

  Objs <- vector(mode="list", length=length(files))
  # read in fold runs
  for (i in 1:length(files)) {

    orig_path <- paste0(data_dir, files[i])

    Objs[[i]] <- readRDS(orig_path)

  }
  names(Objs) <- seeds

  return(Objs)
}




#--------------------------------------------------------
# Fix COMPASS-covariate output into a COMPASS-like object ----
# @param: uncleaned covar obj, cleaned orig obj, covariate matrix
#--------------------------------------------------------
COVtoCOMPASSObj <- function(covarObj, origObj, covMat) {
  covarNew <- list(fit=covarObj$COMPASScovarobj, data = origObj$data)
  # add covariate matrix
  covarNew$data$X <- covMat
  # make sample ID rownames in meta
  rownames(covarNew$data$meta) <- covarNew$data$meta[,covarNew$data$individual_id]
  # make beta gate rownames
  rownames(covarNew$fit$beta) <- colnames(covarNew$fit$mean_gamma)[1:(ncol(covarNew$fit$mean_gamma)-1)]
  # make beta covar colnames
  colnames(covarNew$fit$beta) <- c("Intercept", colnames(covMat))

  return(covarNew)
}

#--------------------------------------------------------
# Get Polyfunctional Score for a covariate  or discrete run ----
# @param: cleaned COMPASS obj, name of ID column, name of cohort column
#--------------------------------------------------------
# getPFSDF <- function(Obj, ID, Cohort) {
#   PFSdf <- if(Obj$fit$model=="discrete") {
#                   data.frame(PFSOrig=PolyfunctionalityScore(Obj))}
#   else {data.frame(PFSCovar=PolyfunctionalityScore.COMPASSResult(Obj))}
#
#     PFSdf %>%
#     rownames_to_column(var="ID") %>%
#     pivot_longer(cols = -ID,
#                  names_to = "COMPASS.Type",
#                  values_to = "PFS") %>%
#     left_join(Obj$data$meta %>% select(matches(ID), matches(Cohort)) %>%
#                 `colnames<-`(c("ID","Cohort")))
#
# }

#TODO: fix this version 2 potentially
#--------------------------------------------------------
# Get Polyfunctional Score for a covariate  or discrete run ----
#--------------------------------------------------------
getPFSMetaDF <- function(Obj) {
  PFSdf <- if(Obj$fit$model=="discrete") {
    data.frame(PFSOrig=PolyfunctionalityScore(Obj),
               FSOrig=FunctionalityScore(Obj))}
  else {data.frame(PFSCovar=PolyfunctionalityScore.COMPASSResult(Obj),
                   FSCovar=FunctionalityScore.COMPASSResult(Obj))}

if(Obj$fit$model=="discrete") {
  return(PFSdf %>%
    bind_cols(Obj$data$meta) %>%
      as.data.frame())
}
  else {
    return(PFSdf %>%
      bind_cols(Obj$data$X) %>%
      bind_cols(Obj$data$meta) %>%
        as.data.frame())}

}

#--------------------------------------------------------
# Put multiple condition PFS scores in dataframe ----
# @params: original Obj, list of covar obj, vector of names, list of cov matrices
#--------------------------------------------------------
# aggPFS <- function(origObj, covarObjList, covarNames, covarMatList) {
#   origPFS <- getPFSDF(origObj)
#   PFSs <- do.call(rbind, lapply(1:length(covarObjList),
#                                 function(i) getPFSDF(covarObjList[[i]]) %>%
#                                   mutate(Type=covarNames[i]))) %>%
#     rbind(getPFSDF(origObj) %>% mutate(Type="Original"))
# }


#--------------------------------------------------------
# Plot paired boxplots of multiple condition PFS ----
# @params: original Obj, list of covar obj, vector of names, list of cov matrices
#--------------------------------------------------------
# boxplotPFS <- function(origObj, covarObjList, covarNames, covarMatList) {
#   PFSs <- aggPFS(origObj, covarObjList, covarNames, covarMatList)
#
#   PFSs %>% filter(Type!="Original") %>%
#     # finangling because we want original to be in all facets
#     mutate(TypeFacet = Type) %>%
#     bind_rows(expand_grid(PFSs %>%
#                             filter(Type=="Original"),
#                           covarNames) %>%
#                 rename(TypeFacet=covarNames)) %>%
#     mutate(Type = factor(Type,
#                          levels=c("Original", covarNames))) %>%
#     ggplot(aes(x=Type, y=PFS)) + geom_boxplot() +
#     geom_line(aes(group=as.factor(ID), color=Cohort),alpha=0.5) +
#     facet_wrap(Cohort~TypeFacet, scales = "free_x", nrow=2)
# }

#########################################################
# Train-Test Evals ----
# 1. We first split and save train-test pairs(meta, covariate matrices, COMPASS train objects)
#    We have to predict mean-gamma for test observations using beta from train and covariates from test.
# 2. We then train a regression of cohort based on age, sex, and PFS.
# 3. We finally calculate ROC for test observations.
#########################################################

# inverse logit function
invlogit <- function(x) {
  exp(x)/(1+exp(x))
}

#--------------------------------------------------------
# Multiply the surface marker observations for each patient ----
# by mean beta to get XB, then take inverse logit to get mean gamma
# @param: COMPASSCov object, beta dataframe (can be from another object)
#--------------------------------------------------------

getPredGamma <- function(Obj, betaDF) {
  #XB <- sweep(sampledBeta[,(ncol(sampledBeta)/2 +1):(ncol(sampledBeta))], MARGIN=1, testMeta[,2], "*")
  X <- Obj$data$X
  XB <- cbind(rep(1, nrow(X)),as.matrix(X))%*%t(betaDF)
  pred <- invlogit(XB)
  Obj$fit$mean_gamma <- cbind(pred, rep(1, nrow(pred)))
  return(Obj)
}

#--------------------------------------------------------
# Use mean gamma to get PFS score (note mean gamma is missing null gate) ----
# @param:
#--------------------------------------------------------
# PolyfunctionalityScore_gamMean <- function(x,  mean_gamma) {
#
#   # remove the 0-value category from degree
#   degree <- head(x$categories[, "Counts"],-1)
#   #degree <- x$categories[, "Counts"]
#   n <- ncol(x$categories) - 1
#   y <- mean_gamma
#   pfs = apply(y, 1, function(row) {
#     ## (2 / (n+1)) is a factor that normalized the score between 0 and 1
#     sum(row * degree / choose(n, degree)) / n * (2 / (n + 1))
#   })
#   return(pfs)
# }

#--------------------------------------------------------
# Scale and center test covs according to training mean/sd ----
# @param: raw train covar matrix, raw test covar matrix
#--------------------------------------------------------
scaleCenterTest <- function(trainCov, testRawCov) {
  meanCovs <- colMeans(trainCov) %>% as.matrix() %>% t()
  sdCovs <- apply(trainCov, 2, sd)

  print(paste("Mean", meanCovs))
  print(paste("SD", sdCovs))

  sweep(testRawCov, 2, meanCovs, FUN="-") %>%
    sweep(2, sdCovs, FUN="/")
}


#--------------------------------------------------------
# Use PCA train loadings to get test PCs ----
# @param: raw train covar matrix, raw test covar matrix
#--------------------------------------------------------
PCATest <- function(trainCov, testRawCov) {
  p <- prcomp(trainCov, center=FALSE)
  as.matrix(testRawCov)%*%as.matrix(p$rotation %>% as.data.frame() %>%
                                      select(PC1:PC10))

}

#--------------------------------------------------------
# Put together test regression input ----
# @params: fold list, fold number, full meta DF (minus placebos), test covar DF, train COMPASScovarobj
# @return: test meta DF with PFS
#--------------------------------------------------------
# testMetawPFS <- function(foldList, foldNum, fullMeta, testCovar, COMPASStrainobj, ID, Cohort, summaryM="mean") {
#   # mean beta
#   if(summaryM=="mean") {
#     meanBeta <- apply(COMPASStrainobj$fit$beta, c(1,2), mean)
#   }
#   else if (summaryM=="median") {
#     meanBeta <- apply(COMPASStrainobj$fit$beta, c(1,2), median)
#   }
#
#   else if (summaryM=="mode") {
#     Mode <- function(x) {
#       ux <- unique(x)
#       ux[which.max(tabulate(match(x, ux)))]
#     }
#
#     meanBeta <- apply(COMPASStrainobj$fit$beta, c(1,2), Mode)
#   }
#
#   # test observation surface markers and meta information
#   #testMeta <- fullMeta[-foldList[[foldNum]],]
#
#
#   testGamma <- getPredTestGamma(as.matrix(testCovar), meanBeta, testMeta, ID)
#   testPFS <- PolyfunctionalityScore_gamMean(COMPASStrainobj$data, testGamma)
#   testMeta <- testMeta %>%
#     rename(Cohort=Cohort, ID=ID) %>%
#     # filter out placebos
#     filter(Cohort!="Placebo") %>%
#     inner_join(testPFS %>% as.data.frame() %>% rownames_to_column(var=ID) %>%
#                  `colnames<-`(c("ID", "PFS"))) %>%
#     select(Cohort, PFS)
#
# }

#--------------------------------------------------------
# Train regression, then return model of binary response for COMPASS-covar ----
# @params: clean COMPASScovarobj, ID column name, Cohort column name
# @return: glm object
#--------------------------------------------------------
getTrainReg <- function(COMPASStrainobj, ID, Cohort) {

  # calculate PFS
  PFSDFTrain <- getPFSDF(COMPASStrainobj, ID, Cohort)

  # add PFS to train meta
  # COMPASStrainobj$data$meta <- COMPASStrainobj$data$meta %>%
  #   inner_join(PFSDFTrain)

  CovsTrain <- COMPASStrainobj$data$meta %>%
    rename(Cohort=Cohort, ID=ID) %>%
    inner_join(PFSDFTrain) %>%
    select(Cohort, PFS) %>%
    # get rid of placebos
    filter(Cohort!="Placebo") %>%
    mutate(Cohort_Bin = if_else(Cohort=="Case", 1, 0))


  ModTrain <- glm(Cohort_Bin~ PFS, #Age + Sex +
                  data=CovsTrain,
                  family="binomial")

  return(ModTrain)

}

#--------------------------------------------------------
# Train regression, then predict probability of binary response for COMPASS-covar ----
# @params: train or test, fold list, fold number, full meta DF, test covar DF, COMPASScovarobj, ID field, summary measure for beta
# @return: dataframe of test obs' probabilities of response
#--------------------------------------------------------
getResponseProbs <- function(type="train", foldList, foldNum, fullMeta, testCovar, COMPASStrainobj, ID, Cohort, summaryM) {

  ModTrain <- getTrainReg(COMPASStrainobj, ID, Cohort)

  if (type=="train") {
    Preds <- predict(ModTrain, type="response")
  }
  else {
  MetaTest <- testMetawPFS(foldList, foldNum, fullMeta, testCovar, COMPASStrainobj, ID, Cohort,summaryM)
  Preds <- predict(ModTrain, newdata=MetaTest, type="response")
  }

  return(Preds)

}

#--------------------------------------------------------
# get ROCs ----
# @param: fold list, fold number, full meta, test predicted probs
# assumes only binary cohort vars in fullmeta
#--------------------------------------------------------

getROCs <-function(idxs, fullMeta, preds) {
  filtMeta <- fullMeta[idxs,]
  roc(filtMeta %>%
        dplyr::mutate(Cohort=factor(Cohort,levels=c("Case","Null"))) %>%
        dplyr::mutate(Cohort_Bin=(Cohort=="Case")*1) %>%
        .$Cohort_Bin, preds %>% as.vector())
}

#--------------------------------------------------------
# AUC and data characteristics for ith fold ----
# @param:  i (fold number), uncleaned original and covariate obj lists, train and test covariate mat lists, meta lists, beta summary measure
# @param cont: name of ID column, name of cohort column
# @return:  data frame of i, # infected/uninfected, ROC
#--------------------------------------------------------
fold_result <- function(i, foldsList, origCOMPASSList, covarCOMPASSList, trainCovarMats, testCovarMats, metaDFs, summaryM, ID, Cohort) {

  # Objects of ith place in list
  COMPASSObj <- origCOMPASSList[[i]]$COMPASSobj
  COMPASScovarObj <- covarCOMPASSList[[i]]
  trainCovarMat <- trainCovarMats[[i]]
  testCovarMat <- testCovarMats[[i]]

  metaDF <- metaDFs[[i]] %>%
    rename(ID=ID, Cohort=Cohort)

  # remove placebos from fold list and metaDF
  PlaceboIDs <- metaDF %>%
    filter(Cohort=="Placebo") %>% .$ID

  for(i in 1:length(foldsList)) {
    foldsList[[i]]<-setdiff(foldsList[[i]],which(ID %in% PlaceboIDs))
  }



  # train and test indices
  trainIdx <- foldsList[[i]]
  testIdx <- setdiff(1:nrow(metaDF), trainIdx)

  # cleaned covariate COMPASS object
  COMPASScovarObjclean <- COVtoCOMPASSObj(COMPASScovarObj, COMPASSObj, trainCovarMat)

  # get test meta DF plus PFS
  testMetaDF <- testMetawPFS(foldsList,i, metaDF, testCovarMat, COMPASScovarObjclean, ID, summaryM)
  # get train meta DF plus PFS
  PFSDForig <- getPFSDF(COMPASSObj, ID, Cohort)
  trainMetaDF <- COMPASSObj$data$meta %>%
    rename(Cohort=Cohort, ID=ID) %>%
    inner_join(PFSDForig) %>%
    select(Cohort, PFS) %>%
    # get rid of placebos
    filter(Cohort!="Placebo") %>%
    mutate(Cohort_Bin = if_else(Cohort=="Case", 1, 0))

  # train regression model on fold covariate PFS and get probability of hospitalization
  trainMod <- getTrainReg(COMPASScovarObjclean,ID, Cohort)
  # trainMod <- getTrainReg(trainMetaDF)
  # train regression model on fold original PFS and probability of hospitalization
  trainModOG <- glm(Cohort_Bin~ PFS,
                    data=trainMetaDF,
                    family="binomial")
  # PFS p-value, train + test probabilities, beta for PFS
  PFSpval <- broom::tidy(trainMod) %>% filter(term=="PFS") %>% .$p.value
  PFSpvalOG <- broom::tidy(trainModOG) %>% filter(term=="PFS") %>% .$p.value
  testProbs <- getResponseProbs("test",foldsList, i, metaDF, testCovarMat, COMPASScovarObjclean, ID, Cohort, summaryM)
  trainProbs <- getResponseProbs("train",foldsList, i, metaDF, trainCovarMat, COMPASScovarObjclean,ID, Cohort, summaryM)

  # train and test indices without the placebo
  trainIdx <- setdiff(foldsList[[i]], which(metaDF %>% rename(Cohort=Cohort) %>%
                                              .$Cohort=="Placebo"))
  testIdx <- setdiff(setdiff(1:nrow(metaDF), trainIdx),which(metaDF %>% rename(Cohort=Cohort) %>%
                                                               .$Cohort=="Placebo"))

  testROC <- getROCs(testIdx, metaDF %>% rename(Cohort=Cohort), testProbs)
  trainROC <- getROCs(trainIdx, metaDF %>% rename(Cohort=Cohort), trainProbs)

  return(data.frame(i=i,
                    numTestInfect = length(testROC$cases),
                    numTestUninfect = length(testROC$controls),
                    numTrainInfect = length(trainROC$cases),
                    numTrainUninfect = length(trainROC$controls),
                    testROC = as.numeric(testROC$auc),
                    trainROC = as.numeric(trainROC$auc),
                    estBeta = trainMod$coefficients["PFS"],
                    PFSpval=PFSpval,
                    PFSpvalOG=PFSpvalOG))
}

#--------------------------------------------------------
# Get table of mean/median ROCs and p-values from above table ----
#--------------------------------------------------------
avgRegVals <- function(resultsDF, covarType) {
  data.frame(Median=apply(resultsDF[,(ncol(resultsDF)-3):ncol(resultsDF)], 2, median),
             Mean = apply(resultsDF[,(ncol(resultsDF)-3):ncol(resultsDF)], 2, mean)) %>%
    kable(caption=paste0("Median and Mean ROCs/P-values for ", covarType, " Covariates")) %>%
    kable_classic_2
}


