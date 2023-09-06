#########################################################
# Comparison with PI ----
#########################################################

#--------------------------------------------------------
# Function to count the number of "!" characters in a column name ----
#--------------------------------------------------------

count_exclamation_marks <- function(column_name) {
  str_count(column_name, "!")
}

#--------------------------------------------------------
# Function for getting agg percentage matrix from count matrix
# @param: original counts dataframe, number of cytokines tested
# @returns: percentage matrix (each column corresponds to 0, 1, 2... cytokines expressed)
#--------------------------------------------------------

Percent.Matrix <- function(df, numCyt) {
  # initialize new aggregated counts dataframe
  aggregated_counts <- data.frame(matrix(0,
                                         nrow=nrow(df),
                                         ncol=numCyt+1))
  # Group the columns based on the number of "!" characters and sum their counts
  for (i in 1:ncol(df)) {
    num_exclamation_marks <- count_exclamation_marks(colnames(df)[i])
    aggregated_counts[,numCyt-num_exclamation_marks+1] <- aggregated_counts[,numCyt-num_exclamation_marks+1] + df[, i]
  }

  # Convert the count matrix to a percentage matrix by row. Note there are extra columns of 0s.
  percentage_matrix <- t(apply(aggregated_counts %>% as.matrix(), 1, function(row) row / sum(row) * 100))

  return(percentage_matrix)
}


#--------------------------------------------------------
# Function for getting PI score -----
# @params: vector of percentages (adding up to 1), number of cytokines tested, q
# @returns: PI score
#--------------------------------------------------------
PI.Score <- function(perc_vector, cytNum, q) {
  score=0
  # summing scores for each #cytokines-expressed category
  for (i in 1:length(perc_vector)) {
    score=score + perc_vector[i]*((i-1)/cytNum)^q
  }
  return(score)
}

rocs_Calc_PI <- function(newTrainOrigObj, newTestOrigObj, covarNames, cohortName, cohortPosVal, q, model="glm") {

  # calculate PI and append to meta information for original objects
  trainOrigX <- data.frame(newTrainOrigObj$data$meta,
                           PI=apply(Percent.Matrix(newTrainOrigObj$data$n_s, numCyt),
                                    1, function(row) PI.Score(row, numCyt, q)))[,c(covarNames), drop=FALSE]
  testOrigX <- data.frame(newTestOrigObj$data$meta,
                          PI=apply(Percent.Matrix(newTestOrigObj$data$n_s, numCyt),
                                   1, function(row) PI.Score(row, numCyt, q)))[,c(covarNames), drop=FALSE]

  # use the train matrices to make GLM or LASS0 models for the cohort membership
  trainY <-data.frame(Y=(newTrainOrigObj$data$meta[,cohortName]==cohortPosVal)*1)
  testY <-data.frame(Y=(newTestOrigObj$data$meta[,cohortName]==cohortPosVal)*1)

  if(model=="glm") {
    trainOrigMod <- glm(Y ~ .,  data=data.frame(trainY, trainOrigX), family="binomial")    # original
  }
  else if (model=="lasso") {
    cvOGfit <- cv.glmnet(as.matrix(trainOrigX), as.matrix(trainY), family = "binomial", type.measure = "class")
    trainOrigMod <- glmnet(x=as.matrix(trainOrigX),                         # original
                           y=as.matrix(trainY),
                           family = "binomial",
                           lambda=cvOGfit$lambda.min)
  }

  # predict the membership of train and test datasets
  if (model=="glm") {
    trainOrigPreds <- predict(trainOrigMod, type="response")
    testOrigPreds <- predict(trainOrigMod, newdata=data.frame(Y=testY,testOrigX), type="response")
  }
  else if(model=="lasso") {
    trainOrigPreds <- predict(trainOrigMod, newx=as.matrix(trainOrigX),type="response")
    testOrigPreds <- predict(trainOrigMod, newx=as.matrix(testOrigX), type="response")
  }

  # get ROCs
  testOrigROC <- roc(testY$Y, as.vector(testOrigPreds))
  trainOrigROC <- roc(trainY$Y, as.vector(trainOrigPreds))

  return(list(model=trainOrigMod,
              testOrigROC = testOrigROC, trainOrigROC=trainOrigROC))

}

#--------------------------------------------------------
# Output as dataframe: PI
#--------------------------------------------------------
PIRun.Results <- function(foldResult, model="glm") {
  if (model=="glm") {
    return(data.frame(numTrainInfected = length(foldResult$trainOrigROC$cases),
                      numTrainHealthy = length(foldResult$trainOrigROC$controls),
                      numTestInfected = length(foldResult$testOrigROC$cases),
                      numTestHealthy = length(foldResult$testOrigROC$controls),
                      TrainOrigROC = foldResult$trainOrigROC$auc,
                      TestOrigROC = foldResult$testOrigROC$auc))
  }
  else return(data.frame(numTrainInfected = length(foldResult$trainOrigROC$cases),
                         numTrainHealthy = length(foldResult$trainOrigROC$controls),
                         numTestInfected = length(foldResult$testOrigROC$cases),
                         numTestHealthy = length(foldResult$testOrigROC$controls),
                         TrainOrigROC = foldResult$trainOrigROC$auc,
                         TestOrigROC = foldResult$testOrigROC$auc))
}



#--------------------------------------------------------
# Put together test regression input ----
# @params: fold list, fold number, full meta DF, full covar DF, # cytokines, q, ID col name (will change cbind to join in future)
# @return: test meta DF with PFS
#--------------------------------------------------------
testMetawPI <- function(foldList, foldNum, fullMeta, Ns.full, numCyt, q, ID) {

  # test observation surface markers and meta information
  testMeta <- fullMeta[-foldList[[foldNum]],]

  testPI <- apply(Percent.Matrix(Ns.full[-foldList[[foldNum]],], numCyt), 1, function(row) PI.Score(row, numCyt, q))
  testMeta <- testMeta %>%
    cbind(testPI %>% as.data.frame()  %>%
            `colnames<-`(c("PI")))
}


#--------------------------------------------------------
# Put together train regression input ----
# @params: fold list, fold number, full meta DF, full covar DF, # cytokines, q, ID col name (will change cbind to join in future)
# @return: train meta DF with PFS
# Requires positive cases to be denoted as "Case" in column "Cohort", negatives as "Null"
#--------------------------------------------------------
trainMetawPI <- function(foldList, foldNum, fullMeta, Ns.full, trainCov, numCyt, q, ID) {

  # train observation surface markers and meta information
  trainMeta <- fullMeta[foldList[[foldNum]],]

  trainPI <- apply(Percent.Matrix(Ns.full[foldList[[foldNum]],], numCyt), 1, function(row) PI.Score(row, numCyt, q))
  trainMeta <- trainMeta %>%
    cbind(trainPI %>% as.data.frame()  %>%
            `colnames<-`(c("PI")))

  CovsTrain <- trainMeta %>%
    select(Cohort, matches(ID), PI) %>%
    mutate(Cohort_Bin = case_when(Cohort=="Case" ~ 1,
                                  Cohort=="Null" ~0,
                                  TRUE ~ NA)) %>%
    # filter out placebos
    filter(!is.na(Cohort_Bin)) %>%
    # add covariates
    left_join(trainCov %>% rownames_to_column(var=ID))

  X=CovsTrain %>% select(-c(Cohort, matches(ID))) %>% as.matrix()

  return(X)
}

#--------------------------------------------------------
# Train regression, then return model of binary response for COMPASS-covar ----
# @params: fold list, fold number, full meta DF, full covar DF, COMPASScovarobj
# @return: glm object
#--------------------------------------------------------
getTrainRegPI <- function(foldList, foldNum, fullMeta, Ns.full, trainCov, numCyt, q, ID) {

  # calculate PI/train observation surface markers and meta information
  fulldata=trainMetawPI(foldList, foldNum, fullMeta, Ns.full, trainCov, numCyt, q, ID)
  X=fulldata[,colnames(fulldata)!="Cohort_Bin"]
  Y=fulldata[,"Cohort_Bin"]
  # initial fit to get best lambda
  cvfit <- cv.glmnet(as.matrix(X), Y, family = "binomial", type.measure = "class")
  ModTrain <- glmnet(x=X,
                     y=Y,
                     family = "binomial",
                     lambda=cvfit$lambda.1se)


  return(ModTrain)

}

#--------------------------------------------------------
# Train regression, then predict probability of binary response for PI ----
# @params: train or test, fold list, fold number, full meta DF, full covar DF, COMPASScovarobj
# @return: dataframe of test obs' probabilities of response
#--------------------------------------------------------
getResponseProbsPI <- function(type="train", foldList, foldNum, fullMeta, Ns.full, trainCov, testCov, numCyt, q, ID) {

  ModTrain <- getTrainRegPI(foldList, foldNum, fullMeta, Ns.full, trainCov, numCyt, q, ID)

  if (type=="train") {
    fulldata <- trainMetawPI(foldList, foldNum, fullMeta, Ns.full, trainCov, numCyt, q, ID)
    Preds <- predict(ModTrain,
                     newx=fulldata[,colnames(fulldata)!="Cohort_Bin"], type="response")
  }
  else {
    MetaTest <- testMetawPI(foldList, foldNum, fullMeta, Ns.full, numCyt, q, ID) %>%
      select(Cohort, matches(ID), PI) %>%
      mutate(Cohort_Bin = case_when(Cohort=="Case" ~ 1,
                                    Cohort=="Null" ~0,
                                    TRUE ~ NA)) %>%
      # filter out placebos
      filter(!is.na(Cohort_Bin)) %>%
      # add covariates
      left_join(testCov %>% rownames_to_column(var=ID)) %>%
      select(-c(Cohort,Cohort_Bin,matches(ID))) %>%
      as.matrix()
    Preds <- predict(ModTrain, newx=MetaTest, type="response")
  }

  return(Preds)

}

#--------------------------------------------------------
# AUC and data characteristics for ith fold ----
# @param:  i (fold number), uncleaned original and covariate obj lists
# @return:  data frame of i, # infected/uninfected, ROC
#--------------------------------------------------------
fold_resultPI <- function(i, foldsList, fullMeta, Ns.full,trainCov, testCov, numCyt,  q, ID) {


  # get test meta DF plus PFS
  testMetaDF <- testMetawPI(foldsList,i, fullMeta, Ns.full, numCyt, q, ID)

  # train regression model on fold covariate PFS and get probability of hospitalization
  trainMod <- getTrainRegPI(foldsList, i, fullMeta, Ns.full, trainCov, numCyt, q, ID)
  # TODO: add the following selected vars to table?
  as.matrix(coef(trainMod))[,1] ->test
  names(test[test!=0])

  # PFS p-value, train + test probabilities
  #PFSpval <- broom::tidy(trainMod) %>% filter(term=="PI") %>% .$p.value
  testProbs <- getResponseProbsPI("test",foldsList, i, fullMeta, Ns.full,trainCov, testCov, numCyt, q, ID)
  trainProbs <- getResponseProbsPI("train",foldsList, i, fullMeta, Ns.full, trainCov, testCov,numCyt, q, ID)

  # train and test indices without the placebo
  trainIdx <- setdiff(foldsList[[i]], which(fullMeta$Cohort=="Placebo"))
  testIdx <- setdiff(setdiff(1:nrow(fullMeta), trainIdx),which(fullMeta$Cohort=="Placebo"))

  testROC <- getROCs(testIdx, fullMeta, testProbs)
  trainROC <- getROCs(trainIdx, fullMeta, trainProbs)

  return(data.frame(i=i,
                    numTestInfect = length(testROC$cases),
                    numTestUninfect = length(testROC$controls),
                    numTrainInfect = length(trainROC$cases),
                    numTrainUninfect = length(trainROC$controls),
                    testROC = as.numeric(testROC$auc),
                    trainROC = as.numeric(trainROC$auc)
                    #PIpval=PFSpval
  ))
}

#########################################################
# LASSO with Log-fold Change ----
#########################################################

#--------------------------------------------------------
# Train regression, then predict probability of binary response for PI ----
# @params: train or test, fold list, fold number, full meta DF, full covar DF, COMPASScovarobj
# @return: dataframe of test obs' probabilities of response
#--------------------------------------------------------
getResponseProbsGates <- function(type="train", foldList, foldNum, fullMeta, Ns.full, Nu.full, numCyt, ID) {
  logfoldChange <- log2(Ns.full/(Nu.full+0.001))

  trainMeta <- fullMeta[foldList[[foldNum]],]
  testMeta <- fullMeta[-foldList[[foldNum]],]
  trainLFC <- logfoldChange[foldList[[foldNum]],]
  testLFC <- logfoldChange[-foldList[[foldNum]],]

  fulldata <- trainMeta %>%
    `rownames<-`(c(trainMeta$`SAMPLE ID`)) %>%
    select(Cohort) %>%
    cbind(trainLFC) %>%
    mutate(Cohort_Bin=if_else(Cohort=="Hospitalized", 1, 0)) %>%
    select(-Cohort)

  # training the LASSO on the fold training data
  X=fulldata[,colnames(fulldata)!="Cohort_Bin"]
  Y=fulldata[,"Cohort_Bin"]
  # initial fit to get best lambda
  cvfit <- cv.glmnet(as.matrix(X), Y, family = "binomial", type.measure = "class")
  ModTrain <- glmnet(x=X,
                     y=Y,
                     family = "binomial",
                     lambda=cvfit$lambda.1se)

  if (type=="train") {
    newdata <- fulldata %>% as.matrix()
    Preds <- predict(ModTrain,
                     newx=newdata[,colnames(newdata)!="Cohort_Bin"], type="response")
  }
  else {
    MetaTest <-  testMeta %>%
      `rownames<-`(c(testMeta$`SAMPLE ID`)) %>%
      select(Cohort) %>%
      cbind(testLFC) %>%
      mutate(Cohort_Bin=if_else(Cohort=="Hospitalized", 1, 0)) %>%
      select(-c(Cohort, Cohort_Bin)) %>%
      as.matrix()
    Preds <- predict(ModTrain, newx=MetaTest, type="response")
  }

  return(Preds)

}

#--------------------------------------------------------
# AUC and data characteristics for ith fold
# @param:  i (fold number), uncleaned original and covariate obj lists
# @return:  data frame of i, # infected/uninfected, ROC
#--------------------------------------------------------
fold_resultGates <- function(i, foldsList, fullMeta, Ns.full, Nu.full, numCyt, ID) {
  # train and test indices
  trainIdx <- foldsList[[i]]
  testIdx <- setdiff(1:nrow(fullMeta), trainIdx)

  #PFSpval <- broom::tidy(trainMod) %>% filter(term=="PI") %>% .$p.value
  testProbs <- getResponseProbsGates("test",foldsList, i, fullMeta, Ns.full, Nu.full, numCyt, ID)
  trainProbs <- getResponseProbsGates("train",foldsList, i, fullMeta, Ns.full, Nu.full, numCyt, ID)


  testROC <- getROCs(testIdx, fullMeta, testProbs)
  trainROC <- getROCs(trainIdx, fullMeta, trainProbs)

  return(data.frame(i=i,
                    numTestInfect = length(testROC$cases),
                    numTestUninfect = length(testROC$controls),
                    numTrainInfect = length(trainROC$cases),
                    numTrainUninfect = length(trainROC$controls),
                    testROC = as.numeric(testROC$auc),
                    trainROC = as.numeric(trainROC$auc)
                    #PIpval=PFSpval
  ))
}
