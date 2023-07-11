
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
# Original Heatmap
# @param: cr COMPASSResult with fit$mean_gamma, fit$categories, data$meta
#--------------------------------------------------------
# Like plot.COMPASSResult, but using the ComplexHeatmap package
plot.COMPASSResult.ComplexHeatmap <- function(cr,
                                              row_annotation = NULL,
                                              cytokine_order_for_annotation = NULL,
                                              dichotomize_by_cytokine = NULL,
                                              dichotomize_by_cytokine_color = NULL,
                                              row_annotation_colors = NULL,
                                              staircase_cytokine_annotation = TRUE,
                                              row_gap = unit(0, "in")) {
  library(tidyverse)
  library(ComplexHeatmap)
  library(COMPASS)
  library(RColorBrewer)

  if(!is.null(row_annotation) & is.null(row_annotation_colors)) { stop("row_annotation_colors must not be NULL if row_annotation is not NULL")}

  mean_gamma <- cr$fit$mean_gamma
  cats <- as.data.frame(cr$fit$categories[, -ncol(cr$fit$categories),drop=FALSE]) # drop the "Counts" column
  numMarkers <- ncol(cats)
  rownames(cats) <- colnames(mean_gamma)
  meta <- cr$data$meta# %>% `rownames<-`(.$`SAMPLE ID`) #TODO: change later assuming rownames are already IDs

  # Filter the subsets to those where the average mean_gamma is greater than the threshold (default in heatmap and this function is 0.01)
  compassSubsetsFiltered <- names(which(apply(mean_gamma, 2, function(x) { mean(x, na.rm = TRUE) }) > 0.01))
  # And remove the subset with 0 positive markers
  compassSubsetsFiltered <- compassSubsetsFiltered[lengths(regmatches(compassSubsetsFiltered, gregexpr("!", compassSubsetsFiltered))) != numMarkers]

  # Subset the cats rows to compassSubsetsFiltered, and put the columns in the order of cytokine_order_for_annotation
  cats <- cats[compassSubsetsFiltered,]
  if(!is.null(cytokine_order_for_annotation)) {
    cats <- cats[, cytokine_order_for_annotation]
  }

  # Before plotting, put the categories data frame rows in the desired order (columns were already re-ordered above)
  # This is essentially the same code I added to the pheatmap function
  cats <- cats[rev(do.call(order, cats)),,drop=FALSE]
  # And order the cats df rows by degrees (number of cytokines in subset)
  ckr<-apply(cats,1,function(x)sum(as.numeric(as.character(x))))
  cats = cats[order(ckr),]
  if(!is.null(dichotomize_by_cytokine)) {
    # And then dichotomize the cats df rows so that all subsets containing the cytokine in dichotomize_by_cytokine (e.g. "IFNg") appear last
    cats <- cats[order(cats[,dichotomize_by_cytokine]),]
  }

  # Now re-order the columns of mean_gamma to match the rows of cats
  mean_gamma <- mean_gamma[,rownames(cats)]

  # If dichotomizing and coloring by a cytokine, specify which cells in the cats matrix should be colored differently
  if(!is.null(dichotomize_by_cytokine) & !is.null(dichotomize_by_cytokine_color)) {
    current_cats_rownames <- rownames(cats)
    cats <- cats %>%
      mutate_all(~ ifelse(. == 1 & !!as.symbol(dichotomize_by_cytokine) > 0, 2, .))
    rownames(cats) <- current_cats_rownames
  }

  # Now re-order the rows of mean_gamma by FunctionalityScore
  FS_order <- order(FunctionalityScore(mean_gamma, n = numMarkers), decreasing=T)
  mean_gamma <- mean_gamma[FS_order,]
  # And make sure to update the metadata rows
  meta <- meta[FS_order,]
  stopifnot(all.equal(rownames(mean_gamma), rownames(meta)))

  # Now additionally order the rows of meta and then mean_gamma based on row_annotation
  suppressWarnings(if(!is.null(row_annotation)) {
    if(length(row_annotation) == 1) {
      meta <- meta %>% arrange(!! rlang::sym(row_annotation))
    } else if(length(row_annotation) > 1) {
      meta <- meta %>% arrange(!!! rlang::syms(row_annotation))
    }
    mean_gamma <- mean_gamma[match(rownames(meta), rownames(mean_gamma)),]
  })

  ht_opt$simple_anno_size = unit(2.5, "mm")
  heatmap_main <- Heatmap(mean_gamma,
                          cluster_rows = FALSE,
                          show_row_dend = FALSE,
                          cluster_columns = FALSE, cluster_column_slices = FALSE,
                          row_split = suppressWarnings(if(!is.null(row_annotation)) { meta %>% dplyr::select(all_of(row_annotation)) } else { NULL }),
                          border = "black",
                          show_column_names = FALSE,
                          show_row_names = FALSE, row_title = NULL,
                          right_annotation = suppressWarnings(if(!is.null(row_annotation)) {
                            HeatmapAnnotation(df = meta %>% dplyr::select(all_of(row_annotation)),
                                              col = row_annotation_colors,
                                              which = "row", border = F,
                                              show_annotation_name = F,
                                              annotation_legend_param = list(border=F,
                                                                             title_gp = gpar(fontsize = 13, fontface = "bold"),
                                                                             labels_gp = gpar(fontsize = 12)))
                          } else {
                            NULL
                          }),
                          heatmap_legend_param = list(border=F, fontsize = 9),
                          name="Response", use_raster = F,
                          col = colorRampPalette(brewer.pal(9,"Purples"))(20),
                          height = unit(5, "in"),
                          row_gap = row_gap)
  # draw(heatmap_main)

  heatmap_cats <- Heatmap(cats %>% dplyr::select(rev(everything())) %>% t(),
                          cluster_rows = FALSE,
                          show_row_dend = FALSE,
                          cluster_columns = FALSE,
                          border = T,
                          show_column_names = FALSE,
                          show_row_names = TRUE,
                          row_names_gp = gpar(fontsize=14),
                          row_names_side = "left",
                          row_title = NULL,
                          col = if(!is.null(dichotomize_by_cytokine) & !is.null(dichotomize_by_cytokine_color)) {
                            c("white", "black", dichotomize_by_cytokine_color) } else { c("white", "black") },
                          show_heatmap_legend = FALSE,
                          rect_gp = gpar(col = "white", lwd = 1),
                          height = unit(1.2, "in"))
  # draw(heatmap_cats)

  # draw(heatmap_main %v% heatmap_cats, gap = unit(0.1, "in"), merge_legends = T) # heatmap_legend_side = "left"
  # ht_opt(RESET = TRUE)

  heatmap_main %v% heatmap_cats
}


#--------------------------------------------------------
# Covariate Heatmap: add horizontal (beta) and vertical (covar value) annots
# @param: cr COMPASSResult with fit$mean_gamma, fit$categories, data$meta
#--------------------------------------------------------
# Like plot.COMPASSResult, but using the ComplexHeatmap package
plot.COMPASSCovarResult.ComplexHeatmap <- function(cr,
                                              covar_to_plot = NULL,
                                              row_annotation = NULL,
                                              cytokine_order_for_annotation = NULL,
                                              dichotomize_by_cytokine = NULL,
                                              dichotomize_by_cytokine_color = NULL,
                                              row_annotation_colors = NULL,
                                              staircase_cytokine_annotation = TRUE,
                                              row_gap = unit(0, "in")) {
  library(tidyverse)
  library(ComplexHeatmap)
  library(COMPASS)
  library(RColorBrewer)

  if(!is.null(row_annotation) & is.null(row_annotation_colors)) { stop("row_annotation_colors must not be NULL if row_annotation is not NULL")}

  mean_gamma <- cr$fit$mean_gamma
  cats <- as.data.frame(cr$fit$categories[, -ncol(cr$fit$categories),drop=FALSE]) # drop the "Counts" column
  numMarkers <- ncol(cats)
  rownames(cats) <- colnames(mean_gamma)
  meta <- cr$data$meta #%>% `rownames<-`(.$`SAMPLE ID`) #TODO: change later assuming rownames are already IDs
  # add beta
  mean_beta <- rowMeans(cr$fit$beta, dims=2)


  # Filter the subsets to those where the average mean_gamma is greater than the threshold (default in heatmap and this function is 0.01)
  compassSubsetsFiltered <- names(which(apply(mean_gamma, 2, function(x) { mean(x, na.rm = TRUE) }) > 0.01))
  # And remove the subset with 0 positive markers
  compassSubsetsFiltered <- compassSubsetsFiltered[lengths(regmatches(compassSubsetsFiltered, gregexpr("!", compassSubsetsFiltered))) != numMarkers]

  # Subset the cats rows to compassSubsetsFiltered, and put the columns in the order of cytokine_order_for_annotation
  cats <- cats[compassSubsetsFiltered,]
  if(!is.null(cytokine_order_for_annotation)) {
    cats <- cats[, cytokine_order_for_annotation]
  }

  # Before plotting, put the categories data frame rows in the desired order (columns were already re-ordered above)
  # This is essentially the same code I added to the pheatmap function
  cats <- cats[rev(do.call(order, cats)),,drop=FALSE]
  # And order the cats df rows by degrees (number of cytokines in subset)
  ckr<-apply(cats,1,function(x)sum(as.numeric(as.character(x))))
  cats = cats[order(ckr),]
  if(!is.null(dichotomize_by_cytokine)) {
    # And then dichotomize the cats df rows so that all subsets containing the cytokine in dichotomize_by_cytokine (e.g. "IFNg") appear last
    cats <- cats[order(cats[,dichotomize_by_cytokine]),]
  }

  # Now re-order the columns of mean_gamma and beta to match the rows of cats
  mean_gamma <- mean_gamma[,rownames(cats)]
  mean_beta_to_plot <- as.matrix(mean_beta[rownames(cats),covar_to_plot])
  colnames(mean_beta_to_plot) <- covar_to_plot

  # If dichotomizing and coloring by a cytokine, specify which cells in the cats matrix should be colored differently
  if(!is.null(dichotomize_by_cytokine) & !is.null(dichotomize_by_cytokine_color)) {
    current_cats_rownames <- rownames(cats)
    cats <- cats %>%
      mutate_all(~ ifelse(. == 1 & !!as.symbol(dichotomize_by_cytokine) > 0, 2, .))
    rownames(cats) <- current_cats_rownames
  }

  # Now re-order the rows of mean_gamma by FunctionalityScore
  FS_order <- order(FunctionalityScore(mean_gamma, n = numMarkers), decreasing=T)
  mean_gamma <- mean_gamma[FS_order,]
  # And make sure to update the metadata rows
  meta <- meta[FS_order,]
  stopifnot(all.equal(rownames(mean_gamma), rownames(meta)))

  # Now additionally order the rows of meta and then mean_gamma based on row_annotation
  suppressWarnings(if(!is.null(row_annotation)) {
    if(length(row_annotation) == 1) {
      meta <- meta %>% arrange(!! rlang::sym(row_annotation))
    } else if(length(row_annotation) > 1) {
      meta <- meta %>% arrange(!!! rlang::syms(row_annotation))
    }
    mean_gamma <- mean_gamma[match(rownames(meta), rownames(mean_gamma)),]
  })

  ht_opt$simple_anno_size = unit(2.5, "mm")
  heatmap_main <- Heatmap(mean_gamma,
                          cluster_rows = FALSE,
                          show_row_dend = FALSE,
                          cluster_columns = FALSE, cluster_column_slices = FALSE,
                          row_split = suppressWarnings(if(!is.null(row_annotation)) { meta %>% dplyr::select(all_of(row_annotation)) } else { NULL }),
                          border = "black",
                          show_column_names = FALSE,
                          show_row_names = FALSE, row_title = NULL,
                          right_annotation = suppressWarnings(if(!is.null(row_annotation)) {
                            HeatmapAnnotation(df = meta %>% dplyr::select(all_of(row_annotation)),
                                              col = row_annotation_colors,
                                              which = "row", border = F,
                                              show_annotation_name = F,
                                              annotation_legend_param = list(border=F,
                                                                             title_gp = gpar(fontsize = 13, fontface = "bold"),
                                                                             labels_gp = gpar(fontsize = 12)))
                          } else {
                            NULL
                          }),
                          bottom_annotation = suppressWarnings(if(!is.null(covar_to_plot)) {
                            HeatmapAnnotation(covar_to_plot = mean_beta_to_plot, #%>% as.data.frame(),
                                              col = list(covar_to_plot =colorRamp2(c(floor(min(mean_beta_to_plot)), 0, ceiling(max(mean_beta_to_plot))),
                                                               c("blue", "white", "red"))),
                                                                which="column",
                                                                show_legend = T,
                                                                #`Effects`=anno_barplot(mean_beta_to_plot),
                                                                show_annotation_name = T,
                                              annotation_name_side = "left",
                                              annotation_name_gp = gpar(fontsize=6),
                                              annotation_legend_param = list(border=F,
                                                                             title_gp = gpar(fontsize = 9, fontface = "bold"),
                                                                             labels_gp = gpar(fontsize = 9)))
                            } else {
                              NULL
                            }),
                          heatmap_legend_param = list(border=F, fontsize = 9),
                          name="Response", use_raster = F,
                          col = colorRampPalette(brewer.pal(9,"Purples"))(20),
                          height = unit(5, "in"),
                          row_gap = row_gap)
  # draw(heatmap_main)

  heatmap_cats <- Heatmap(cats %>% dplyr::select(rev(everything())) %>% t(),
                          cluster_rows = FALSE,
                          show_row_dend = FALSE,
                          cluster_columns = FALSE,
                          border = T,
                          show_column_names = FALSE,
                          show_row_names = TRUE,
                          row_names_gp = gpar(fontsize=14),
                          row_names_side = "left",
                          row_title = NULL,
                          col = if(!is.null(dichotomize_by_cytokine) & !is.null(dichotomize_by_cytokine_color)) {
                            c("white", "black", dichotomize_by_cytokine_color) } else { c("white", "black") },
                          show_heatmap_legend = FALSE,
                          rect_gp = gpar(col = "white", lwd = 1),
                          height = unit(1.2, "in"))
  # draw(heatmap_cats)

  # draw(heatmap_main %v% heatmap_cats, gap = unit(0.1, "in"), merge_legends = T) # heatmap_legend_side = "left"
  # ht_opt(RESET = TRUE)

  heatmap_main %v% heatmap_cats
}

#--------------------------------------------------------
# Fix COMPASS-covariate output into a COMPASS-like object
#--------------------------------------------------------
COVtoCOMPASSObj <- function(covarObj, origObj, covMat) {
  covarNew <- list(fit=covarObj$COMPASScovarobj, data = origObj$COMPASSobj$data)
  # make sample ID rownames in meta
  rownames(covarNew$data$meta) <- covarNew$data$meta[,covarNew$data$individual_id]
  # make beta gate rownames
  rownames(covarNew$fit$beta) <- colnames(covarNew$fit$mean_gamma)[1:(ncol(covarNew$fit$mean_gamma)-1)]
  # make beta covar colnames
  colnames(covarNew$fit$beta) <- c("Intercept", colnames(covMat))

  return(covarNew)
}

#--------------------------------------------------------
# Get Polyfunctional Score for a covariate  or discrete run
#--------------------------------------------------------
getPFSDF <- function(Obj) {
  PFSdf <- if(Obj$fit$model=="discrete") {
                  data.frame(PFSOrig=PolyfunctionalityScore(Obj))}
  else {data.frame(PFSCovar=PolyfunctionalityScore.COMPASSResult(Obj))}

    PFSdf %>%
    rownames_to_column(var="ID") %>%
    pivot_longer(cols = -ID,
                 names_to = "COMPASS.Type",
                 values_to = "PFS") %>%
    left_join(Obj$data$meta %>% select(`SAMPLE ID`, Cohort) %>%
                rename(ID=`SAMPLE ID`))

}

#########################################################
# Train-Test Evals
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
# Multiply the surface marker observations for each test patient
# by mean beta to get XB, then take inverse logit to get mean gamma
#--------------------------------------------------------

getPredTestGamma <- function(testCovarDF, meanBetaDF, metaDF, sampName) {
  #XB <- sweep(sampledBeta[,(ncol(sampledBeta)/2 +1):(ncol(sampledBeta))], MARGIN=1, testMeta[,2], "*")
  XB <- cbind(rep(1, nrow(testCovarDF)),testCovarDF)%*%t(meanBetaDF)
  pred <- invlogit(XB)
  rownames(pred) <- metaDF[,sampName]
  return(pred)
}

#--------------------------------------------------------
# Use mean gamma to get PFS score
#--------------------------------------------------------
PolyfunctionalityScore_gamMean <- function(x,  mean_gamma) {

  # remove the 0-value category from degree
  degree <- head(x$categories[, "Counts"],-1)
  n <- ncol(x$categories) - 1
  y <- mean_gamma
  pfs = apply(y, 1, function(row) {
    ## (2 / (n+1)) is a factor that normalized the score between 0 and 1
    sum(row * degree / choose(n, degree)) / n * (2 / (n + 1))
  })
  return(pfs)
}

#--------------------------------------------------------
# Scale and center test covs according to training mean/sd
#--------------------------------------------------------
scaleCenterTest <- function(trainCov, testRawCov) {
  meanCovs <- colMeans(trainCov) %>% as.matrix() %>% t()
  sdCovs <- apply(trainCov, 2, sd)

  sweep(testRawCov, 2, meanCovs, FUN="-") %>%
    sweep(2, sdCovs, FUN="/")
}


#--------------------------------------------------------
# Put together test regression input
# @params: fold list, fold number, full meta DF, full covar DF, train COMPASScovarobj
# @return: test meta DF with PFS
#--------------------------------------------------------
testMetawPFS <- function(foldList, foldNum, fullMeta, fullCovar, COMPASStrainobj, ID) {
  # mean beta
  meanBeta <- apply(COMPASStrainobj$fit$beta, c(1,2), mean)

  # test observation surface markers and meta information
  testMeta <- fullMeta[-foldList[[foldNum]],]
  # calculate correct surface test covariates
  testSurface <- fullCovar[-foldList[[foldNum]],]


  testGamma <- getPredTestGamma(testSurface, meanBeta, testMeta, ID)
  testPFS <- PolyfunctionalityScore_gamMean(COMPASStrainobj$data, testGamma)
  testMeta <- testMeta %>%
    inner_join(testPFS %>% as.data.frame() %>% rownames_to_column(var=ID) %>%
                 `colnames<-`(c(ID, "PFS")))
}

#--------------------------------------------------------
# Train regression, then return model of binary response for COMPASS-covar
# @params: fold list, fold number, full meta DF, full covar DF, COMPASScovarobj
# @return: glm object
#--------------------------------------------------------
getTrainReg <- function(foldList, foldNum, fullMeta, fullCovar, COMPASStrainobj, ID) {

  # calculate PFS
  PFSDFTrain <- getPFSDF(COMPASStrainobj)

  # add PFS to train meta
  COMPASStrainobj$data$meta <- COMPASStrainobj$data$meta %>%
    inner_join(PFSDFTrain %>% rename("SAMPLE ID"=ID))

  CovsTrain <- COMPASStrainobj$data$meta %>%
    select(Cohort, Age, Sex, PFS) %>%
    mutate(Cohort_Bin = if_else(Cohort=="Hospitalized", 1, 0))

  ModTrain <- glm(Cohort_Bin~ PFS, #Age + Sex +
                  data=CovsTrain,
                  family="binomial")

  return(ModTrain)

}

#--------------------------------------------------------
# Train regression, then predict probability of binary response for COMPASS-covar
# @params: train or test, fold list, fold number, full meta DF, full covar DF, COMPASScovarobj
# @return: dataframe of test obs' probabilities of response
#--------------------------------------------------------
getResponseProbs <- function(type="train", foldList, foldNum, fullMeta, fullCovar, COMPASStrainobj, ID) {

  ModTrain <- getTrainReg(foldList, foldNum, fullMeta, fullCovar, COMPASStrainobj, ID)

  if (type=="train") {
    Preds <- predict(ModTrain, type="response")
  }
  else {
  MetaTest <- testMetawPFS(foldList, foldNum, fullMeta, fullCovar, COMPASStrainobj, ID)
  Preds <- predict(ModTrain, newdata=MetaTest, type="response")
  }

  return(Preds)

}

#--------------------------------------------------------
# get ROCs
# @param: fold list, fold number, full meta, test predicted probs
#--------------------------------------------------------

getROCs <-function(idxs, fullMeta, preds) {
  filtMeta <- fullMeta[idxs,]
  roc(filtMeta %>%
        mutate(Cohort_Bin=(Cohort=="Hospitalized")*1) %>%
        .$Cohort_Bin, preds)
}

#--------------------------------------------------------
# AUC and data characteristics for ith fold
# @param:  i (fold number)
# @return:  data frame of i, # infected/uninfected, ROC
#--------------------------------------------------------
fold_result <- function(i, foldsList, origCOMPASSList, covarCOMPASSList, covarMats, metaDFs) {

  # Objects of ith place in list
  COMPASSObj <- origCOMPASSList[[i]]
  COMPASScovarObj <- covarCOMPASSList[[i]]
  covarMat <- covarMats[[i]]
  metaDF <- metaDFs[[i]]

  # train and test indices
  trainIdx <- foldsList[[i]]
  testIdx <- setdiff(1:nrow(metaDF), trainIdx)

  # cleaned covariate COMPASS object
  COMPASScovarObjclean <- COVtoCOMPASSObj(COMPASScovarObj, COMPASSObj, covarMat)

  # get test meta DF plus PFS
  testMetaDF <- testMetawPFS(foldsList,i, metaDF, covarMat, COMPASScovarObjclean, "SAMPLE ID")
  # get train meta DF plus PFS
  PFSDForig <- getPFSDF(COMPASSObj$COMPASSobj)
  trainMetaDF <- COMPASSObj$COMPASSobj$data$meta %>%
    inner_join(PFSDForig %>% rename("SAMPLE ID"=ID)) %>%
    select(Cohort, Age, Sex, PFS) %>%
    mutate(Cohort_Bin = if_else(Cohort=="Hospitalized", 1, 0))

  # train regression model on fold covariate PFS and get probability of hospitalization
  trainMod <- getTrainReg(foldsList, i, metaDF, covarMat, COMPASScovarObjclean,"SAMPLE ID")
  # train regression model on fold original PFS and probability of hospitalization
  trainModOG <- glm(Cohort_Bin~ PFS,
                    data=trainMetaDF,
                    family="binomial")
  # PFS p-value, train + test probabilities
  PFSpval <- broom::tidy(trainMod) %>% filter(term=="PFS") %>% .$p.value
  PFSpvalOG <- broom::tidy(trainModOG) %>% filter(term=="PFS") %>% .$p.value
  testProbs <- getResponseProbs("test",foldsList, i, metaDF, covarMat, COMPASScovarObjclean,"SAMPLE ID")
  trainProbs <- getResponseProbs("train",foldsList, i, metaDF, covarMat, COMPASScovarObjclean,"SAMPLE ID")


  testROC <- getROCs(testIdx, metaDF, testProbs)
  trainROC <- getROCs(trainIdx, metaDF, trainProbs)

  return(data.frame(i=i,
                    numTestInfect = length(testROC$cases),
                    numTestUninfect = length(testROC$controls),
                    numTrainInfect = length(trainROC$cases),
                    numTrainUninfect = length(trainROC$controls),
                    testROC = as.numeric(testROC$auc),
                    trainROC = as.numeric(trainROC$auc),
                    PFSpval=PFSpval,
                    PFSpvalOG=PFSpvalOG))
}

#--------------------------------------------------------
# Get table of mean/median ROCs and p-values from above table
#--------------------------------------------------------
avgRegVals <- function(resultsDF, covarType) {
  data.frame(Median=apply(resultsDF[,(ncol(resultsDF)-3):ncol(resultsDF)], 2, median),
             Mean = apply(resultsDF[,(ncol(resultsDF)-3):ncol(resultsDF)], 2, mean)) %>%
    kable(caption=paste0("Median and Mean ROCs/P-values for ", covarType, " Covariates")) %>%
    kable_classic_2
}


