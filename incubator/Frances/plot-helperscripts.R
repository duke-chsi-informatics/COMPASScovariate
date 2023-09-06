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
# Save heatmap
# @param: cleaned covariate obj,
#      path to file, row annotation column name, named list of colors
#--------------------------------------------------------
plotHeatM <- function(Obj, pathtofile, row_annot, row_annot_colors) {

  mean_beta <- rowMeans(Obj$fit$beta, dim=2)

  #--------------------------------------------------------
  # saving covar heatmap
  #--------------------------------------------------------

  png(file=pathtofile,
      width=500, height = 750)
  covarplot <- plot.COMPASSCovarResult.ComplexHeatmap(Obj,
                                                      covar_to_plot = colnames(mean_beta),
                                                      row_annotation = row_annot,
                                                      row_annotation_colors = row_annot_colors)

  draw(covarplot)
  dev.off()
}
