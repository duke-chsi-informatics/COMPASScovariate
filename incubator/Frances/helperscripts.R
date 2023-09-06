
#--------------------------------------------------------
# Read in multiple COMPASS seeded runs
# @param: data directory path, pattern
# @return: list of uncleaned COMPASS/COMPASS-covar objects
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
# @return: cleaned COMPASS-covar object with X added
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
# @param: clean COMPASS or COMPASS-covariate object
# @return: DF of PFS/FS scores, meta data, and covariates (X)
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
# Remove Placebo rows from COMPASS (Cov) Object
#--------------------------------------------------------
removePlacebos <- function(Obj, Cohort, PlaceboMark) {
  placebo_idx = which(Obj$data$meta[,Cohort]==PlaceboMark)
  nrows = nrow(Obj$data$meta)
  # ignore beta, alpha u, alpha s
  subdf_names = c("data$meta", "data$X", "data$n_s", "data$n_u",
                  "fit$mean_gamma")
  subvec_names = c("data$counts_u", "data$counts_s", "fit$A_alphas", "fit$A_alphau", "fit$A_gamma")
  submat_names = c("fit$gamma")

  for(name in subdf_names) {
    components <- strsplit(name, "\\$")[[1]]  # Split into components
    parent <- components[1]  # Parent sublist name
    child <- components[2]   # Child data frame name

    # Remove the placebo rows from the specified data frame
    Obj[[parent]][[child]] <- Obj[[parent]][[child]][setdiff(1:nrows, placebo_idx), ]
  }

  for(name in subvec_names) {
    components <- strsplit(name, "\\$")[[1]]  # Split into components
    parent <- components[1]  # Parent sublist name
    child <- components[2]   # Child vec name

    # Remove the placebo entries from the specified vec
    Obj[[parent]][[child]] <- Obj[[parent]][[child]][setdiff(1:nrows, placebo_idx)]
  }

  for(name in submat_names) {
    components <- strsplit(name, "\\$")[[1]]  # Split into components
    parent <- components[1]  # Parent sublist name
    child <- components[2]   # Child vec name

    # Remove the placebo rows from the specified matrix
    Obj[[parent]][[child]] <- Obj[[parent]][[child]][setdiff(1:nrows, placebo_idx),,]
  }

  return(Obj)
}



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
# @param: cleaned COMPASS-Cov object, beta dataframe (can be from another object)
# @return: cleaned COMPASS-cov object with new predicted mean gamma
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
# Scale and center test covs according to training mean/sd ----
# @param: raw train covar matrix, raw test covar matrix
# @return: test covar matrix with scaled and centered columns according to training
#--------------------------------------------------------
scaleCenterTest <- function(trainCov, testRawCov) {
  meanCovs <- colMeans(trainCov) %>% as.matrix() %>% t()
  sdCovs <- apply(trainCov, 2, sd)

  sweep(testRawCov, 2, meanCovs, FUN="-") %>%
    sweep(2, sdCovs, FUN="/")
}


#--------------------------------------------------------
# Use PCA train loadings to get test PCs ----
# @param: raw train covar matrix, raw test covar matrix
# @return: the 10-PC representation of the test data
#--------------------------------------------------------
PCATest <- function(trainCov, testRawCov) {
  p <- prcomp(trainCov, center=FALSE)
  as.matrix(testRawCov)%*%as.matrix(p$rotation %>% as.data.frame() %>%
                                      select(PC1:PC10))

}



#--------------------------------------------------------
# Fit train models and test them with the test COMPASS data
# @param: cleaned COMPASS-covar and COMPASS train + test objects,
#         vectors of covariate names to use for COMPASS-cov and COMPASS PFS regression,
#         strings denoting column name of cohort binary designation and what the positive value is,
#         whether to use glm or lasso
# @return: list including the train covariate and original regression models,
#          ROCs of train and test covariate and original data
#--------------------------------------------------------

rocs_Calc <- function(newTrainCovarObj, newTestCovarObj, newTrainOrigObj, newTestOrigObj, covarNamesCov, covarNamesOrig, cohortName, cohortPosVal, model="glm") {
  # calculate PFS/FS and append to meta/cov information for all objects
  trainCovarX <- getPFSMetaDF(newTrainCovarObj)[,c(covarNamesCov),drop=FALSE]
  testCovarX <- getPFSMetaDF(newTestCovarObj)[,c(covarNamesCov), drop=FALSE]
  trainOrigX <- getPFSMetaDF(newTrainOrigObj)[,c(covarNamesOrig), drop=FALSE]
  testOrigX <- getPFSMetaDF(newTestOrigObj)[,c(covarNamesOrig), drop=FALSE]

  # use the train matrices to make GLM or LASS0 models for the cohort membership
  trainY <-data.frame(Y=(newTrainCovarObj$data$meta[,cohortName]==cohortPosVal)*1)
  testY <-data.frame(Y=(newTestCovarObj$data$meta[,cohortName]==cohortPosVal)*1)

  if(model=="glm") {
    trainCovarMod <- glm(Y ~ ., data=data.frame(trainY, trainCovarX), family="binomial")  # covariate
    trainOrigMod <- glm(Y ~ .,  data=data.frame(trainY, trainOrigX), family="binomial")    # original
  }
  else if (model=="lasso") {
    cvCovfit <- cv.glmnet(as.matrix(trainCovarX), as.matrix(trainY),
                          family = "binomial", type.measure = "class")
    trainCovarMod <- glmnet(x=as.matrix(trainCovarX),                        # covariate
                            y=as.matrix(trainY),
                            family = "binomial",
                            lambda=cvCovfit$lambda.min)

    cvOGfit <- cv.glmnet(as.matrix(trainOrigX), as.matrix(trainY), family = "binomial", type.measure = "class")
    trainOrigMod <- glmnet(x=as.matrix(trainOrigX),                         # original
                           y=as.matrix(trainY),
                           family = "binomial",
                           lambda=cvOGfit$lambda.min)
  }

  # predict the membership of train and test datasets
  if (model=="glm") {
    trainCovPreds <- predict(trainCovarMod, type="response")
    testCovPreds <- predict(trainCovarMod, newdata=data.frame(Y=testY,testCovarX), type="response")

    trainOrigPreds <- predict(trainOrigMod, type="response")
    testOrigPreds <- predict(trainOrigMod, newdata=data.frame(Y=testY,testOrigX), type="response")
  }
  else if(model=="lasso") {
    trainCovPreds <- predict(trainCovarMod, newx=as.matrix(trainCovarX), type="response")
    testCovPreds <- predict(trainCovarMod, newx=as.matrix(testCovarX), type="response")

    trainOrigPreds <- predict(trainOrigMod, newx=as.matrix(trainOrigX),type="response")
    testOrigPreds <- predict(trainOrigMod, newx=as.matrix(testOrigX), type="response")
  }

  # get ROCs
  testOrigROC <- roc(testY$Y, as.vector(testOrigPreds))
  trainOrigROC <- roc(trainY$Y, as.vector(trainOrigPreds))
  testCovROC <- roc(testY$Y, as.vector(testCovPreds))
  trainCovROC <- roc(trainY$Y, as.vector(trainCovPreds))

  return(list(covariate_model = trainCovarMod, original_model=trainOrigMod,
              testOrigROC = testOrigROC, trainOrigROC=trainOrigROC,
              testCovROC = testCovROC, trainCovROC = trainCovROC))

}

#--------------------------------------------------------
# For the MCMC-specific preds
#--------------------------------------------------------
rocs_Calc_Avg <- function(newTrainCovarObj, newTestCovarObj, newTrainOrigObj, newTestOrigObj, covarNamesCov, covarNamesOrig, cohortName, cohortPosVal, nSamps, model="glm") {
  # train model outside here
  trainCovarX <- getPFSMetaDF(newTrainCovarObj)[,c(covarNamesCov), drop=FALSE]
  trainY <- data.frame(Y=(newTrainCovarObj$data$meta[,cohortName]==cohortPosVal)*1)
  testY <-data.frame(Y=(newTestCovarObj$data$meta[,cohortName]==cohortPosVal)*1)

  if(model=="glm") {
    trainCovarMod <- glm(Y ~ ., data=data.frame(trainY, trainCovarX), family="binomial")
    trainCovPreds <- predict(trainCovarMod, type="response")
  }
  else if (model=="lasso") {
    cvCovfit <- cv.glmnet(as.matrix(trainCovarX), as.matrix(trainY),
                          family = "binomial", type.measure = "class")
    trainCovarMod <- glmnet(x=as.matrix(trainCovarX),
                            y=as.matrix(trainY),
                            family = "binomial",
                            lambda=cvCovfit$lambda.min)

    trainCovPreds <- predict(trainCovarMod, newx=as.matrix(trainCovarX), type="response")
  }

  probsPos <- one_Iter(newTrainCovarObj, newTestCovarObj, trainCovarMod, covarNamesCov, cohortName, cohortPosVal, nSamps, model="glm")
  testCovROC <- roc(testY$Y, rowMeans(probsPos))

  return(list(probsPos= probsPos, testCovROC=testCovROC))

}


one_Iter <- function(newTrainCovarObj, newTestCovarObj, trainCovarMod, covarNamesCov, cohortName, cohortPosVal, nSamps, model="glm") {
  covarBetas <- newTrainCovarObj$fit$beta
  sampleIdx <- sample(1:dim(covarBetas)[3], nSamps)

  # calculate mean gamma, then PFS using MCMC beta
  probsPos <- do.call(cbind,lapply(sampleIdx, function(i) {
    predTestCovarObj <- getPredGamma(newTestCovarObj, covarBetas[,,i])
    testCovarX <- getPFSMetaDF(predTestCovarObj)[,c(covarNamesCov), drop=FALSE]

    # use the train matrices to make GLM or LASS0 models for the cohort membership
    testY <-data.frame(Y=(newTestCovarObj$data$meta[,cohortName]==cohortPosVal)*1)

    if (model=="glm") {
      testCovPreds <- as.matrix(predict(trainCovarMod, newdata=data.frame(Y=testY,testCovarX), type="response"))
    }
    else if(model=="lasso") {
      testCovPreds <- predict(trainCovarMod, newx=as.matrix(testCovarX), type="response")

    }
  }))

  return(probsPos %>% `colnames<-`(sampleIdx))
}


#--------------------------------------------------------
# Output as dataframe: with COMPASS test run
#--------------------------------------------------------

COMPASSTestRun.Results <- function(foldResult, model="glm") {
  if (model=="glm") {
    return(data.frame(numTrainInfected = length(foldResult$trainCovROC$cases),
                      numTrainHealthy = length(foldResult$trainCovROC$controls),
                      numTestInfected = length(foldResult$testCovROC$cases),
                      numTestHealthy = length(foldResult$testCovROC$controls),
                      TrainCovarROC = foldResult$trainCovROC$auc,
                      TestCovarROC = foldResult$testCovROC$auc,
                      TrainOrigROC = foldResult$trainOrigROC$auc,
                      TestOrigROC = foldResult$testOrigROC$auc))
  }
  else return(data.frame(numTrainInfected = length(foldResult$trainCovROC$cases),
                         numTrainHealthy = length(foldResult$trainCovROC$controls),
                         numTestInfected = length(foldResult$testCovROC$cases),
                         numTestHealthy = length(foldResult$testCovROC$controls),
                         TrainCovarROC = foldResult$trainCovROC$auc,
                         TestCovarROC = foldResult$testCovROC$auc,
                         TrainOrigROC = foldResult$trainOrigROC$auc,
                         TestOrigROC = foldResult$testOrigROC$auc))
}

#--------------------------------------------------------
# Output as dataframe: with Post Preds
#--------------------------------------------------------

COMPASSPostPred.Results <- function(foldResult, model="glm") {
  if (model=="glm") {
    return(data.frame(numTestInfected = length(foldResult$testCovROC$cases),
                      numTestHealthy = length(foldResult$testCovROC$controls),
                      TestCovarROC = foldResult$testCovROC$auc,
                      nSamps = ncol(foldResult$probsPos)))
  }
  else return(data.frame(numTestInfected = length(foldResult$testCovROC$cases),
                         numTestHealthy = length(foldResult$testCovROC$controls),
                         TestCovarROC = foldResult$testCovROC$auc,
                         nSamps = ncol(foldResult$probsPos)))
}

