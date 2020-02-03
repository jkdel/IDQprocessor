#' Plot the assay plate layout.
#'
#' @param eset A Biobase::ExpressionSet.
#' @param highlight An optional expression used to highlight certain samples.
#' @param fill Name of the column used to distinguish samples and wells.
#' @param well_col Name of the column containing the well position.
#' @param batch_ind The name of the column containing the information to use for
#'   determining the batches. Defaults to "Plate Production No.).
#' @return An invisible `ggplot()`.
#' @example
#' \notrun{
#' plot_layout(es, `Sample Identification` == "MMU08_EDTA_96dpi")
#' }
#' @import ggplot2
#' @export
plot_layout <- function(eset,
                        highlight,
                        fill = "Sample Type",
                        well_col = "Well Position",
                        batch_ind = "Plate Production No.") {
  pDat <- Biobase::pData(eset)
  pDat$PCol <- factor(((pDat[[well_col]]-1) %% 12)+1)
  pDat$PRow <- factor(((pDat[[well_col]]-1) %/% 12)+1)
  pDat$PRow <- factor(pDat[["PRow"]], levels = rev(levels(pDat[["PRow"]])))
  pDat$batch <- batches(eset, batch_ind)
  if (!missing(highlight)) {
    pDat$high <- ifelse(eval(substitute(highlight), pDat),"red",NA)
  }
  g <- ggplot(pDat, aes(PCol, PRow, fill = pDat[[fill]])) +
    scale_fill_discrete(fill) +
    scale_color_identity() +
    coord_equal() +
    facet_wrap(. ~ batch) +
    theme_void() +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_legend(ncol = 5)) +
    theme(legend.position = "top")
  if (missing(highlight)) print(g + geom_tile(size = 1.2, width = 0.9, height = 0.9))
  else print(g + geom_tile(size = 1.2, width = 0.9, height = 0.9, color = pDat[["high"]]))
  invisible(g)
}

#' Produce a PCA plot.
#'
#' @param eset A Biobase::ExpressionSet.
#' @param axes Numeric vector of length two specifying which components to draw.
#' @param scaling Logical specifying if scaling should be performed, or function
#'   to use to perform scaling. Defaults to TRUE.
#' @param group The variable to use to group observations.
#' @param ellipse Logical specifying if a normal ellipse should be drawn for
#'   each group.
#' @param ellipse_prob Numerical, size of the ellipse under normal distribution.
#' @param annotations The variable to use to name the points.
#' @return An invisible `ggplot()`.
#' @seealso `ggbiplot::ggbiplot()`
#' @import ggplot2
#' @export
plot_pca <- function(eset,
                     axes = 1:2,
                     scaling = T,
                     group,
                     ellipse = F,
                     ellipse_prob=0.68,
                     annotations) {
  if (length(axes)!=2) stop("`axes` should be a numeric vector of length two.")
  dat <- t(Biobase::exprs(eset))
  sel <- stats::complete.cases(dat)
  if (any(dim(dat[sel, ]) != dim(dat)))
    warning("Incomplete data. Only complete cases were used.", call. = F)
  del <- which(lapply(as.data.frame(dat[sel,]), stats::sd)==0)
  if (length(del)>0) {
    dat <- dat[,-del]
    # sel <- sel[-del]
    warning("Following features were removed because of unit variance: ", names(del), call. = F)
  }
  if (is.function(scaling)) {
    dat <- scaling(dat[sel, ])
    pca <- stats::prcomp(dat, scale. = FALSE)
  } else if (is.logical(scaling)) {
    pca <- stats::prcomp(dat[sel, ], scale. = scaling)
  } else {
    stop("`scaling` appears to be neither logical nor a function.", call. = F)
  }
  perc <- round(pca$sdev[axes]^2*100/sum(pca$sdev^2))
  ppca <- stats::predict(pca)
  cols <- colnames(ppca)[axes]
  ppca <- stats::setNames(as.data.frame(ppca[,axes]),c("x","y"))
  g <- ggplot() +
    coord_equal() +
    labs(x = paste0(cols[[1]]," (",perc[[1]],"%)"), y = paste0(cols[[2]]," (",perc[[2]],"%)"))
  if (!missing(group)) {
    ppca$gp <- eval(substitute(group), Biobase::pData(eset))[sel]
    g <- g + geom_point(data = ppca, aes(x, y, color = gp))
  } else g <- g + geom_point(data = ppca, aes(x, y))
  if (ellipse && length(unique(ppca$gp))>1) {
    # mostly taken from
    # https://github.com/vqv/ggbiplot/blob/master/R/ggbiplot.r
    # fd999409dfc8b52ec8d3ddb6341a14eafa7c7812
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- do.call(rbind.data.frame,lapply(split(ppca[,1:2],ppca$gp), function(x) {
      if(nrow(x) <= 2) return(NULL)
      sigma <- stats::var(x)
      mu <- sapply(x,base::mean)
      ed <- base::sqrt(stats::qchisq(ellipse_prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
    }))
    names(ell)[1:2] <- c("x", "y")
    ell$gp <- rep(unique(ppca$gp),each=100)
    g <- g + geom_path(data = ell, aes(x, y, color = gp, group = gp))
  }
  if (!missing(annotations)) {
    ppca$annot <- eval(substitute(annotations), Biobase::pData(eset))[sel]
    g <- g + geom_text(data = ppca, aes(x, y, label = annot), hjust = 0, nudge_x = 0.05)
  }
  print(g)
  invisible(g)
}

#' Plot features over processing time.
#'
#' @param eset A Biobase::ExpressionSet.
#' @param features Character vector of feature names for which to plot.
#' @param ind_qcs Name of the column to use to identify QCs (any sample
#'   measured at least twice)
#' @param measurement_col The name of the column containing the measurement
#'   time (àefaults to "Measurement Time").
#' @param ref_index The index of the subcolumn to use for ordering of
#'   measurement times (typically 1, and therefore default value).
#' @return An invisible `ggplot()`.
#' @import ggplot2
#' @export
plot_qcs <- function(eset, features, ind_qcs = "Sample Identification",
                     measurement_col = "Measurement Time", ref_index = 1) {
  eset$proc <- processing_order(eset, measurement_col, ref_index)
  # eset$batch <- factor(batches(eset, batch_ind))
  qcs <- names(which(table(eset[[ind_qcs]])>1))
  eset <- eset[features,]
  dat <- cbind.data.frame(Biobase::pData(eset)[,c("proc",ind_qcs)], t(Biobase::exprs(eset)))
  dat <- tidyr::pivot_longer(dat, -1:-2, names_to="Feature", values_to = "Expression")
  dat$QC <- ifelse(dat[[ind_qcs]] %in% qcs, dat[[ind_qcs]], NA)
  # bt <- sapply(split(dat,dat$batch),function(x) {
  #   m <- max(x$proc)
  #   ifelse(m==max(dat$proc), NA, m)
  # })
  g <- ggplot(dat[is.na(dat$QC),], aes(proc, Expression)) +
    geom_point(color="black",alpha=0.2,pch=19) +
    geom_point(data=dat[!is.na(dat$QC),], aes(color=QC),pch=19) +
    geom_line(data=dat[!is.na(dat$QC),], aes(color=QC, group=QC)) +
    # geom_vline(xintercept = bt+0.5) +
    facet_wrap(. ~ Feature, scales = "free_y") +
    guides(color=guide_legend(direction = "horizontal")) +
    theme(legend.position = "top")
  print(g)
  invisible(dat)
}

#' Boxplot of each feature.
#'
#' @param eset A Biobase::ExpressionSet.
#' @param by_batch Logical, `TRUE` if separate facets should be used for each
#'   bacth.
#' @param batch_ind Name of the column from where to derive batch information.
#' @return An invisible `ggplot()`.
#' @import ggplot2
#' @export
plot_features <- function(eset, by_batch = T, batch_ind = "Plate Production No.") {
  dat <- t(Biobase::exprs(eset))
  dat <- cbind.data.frame(dat, batch = batches(eset))
  dat <- tidyr::pivot_longer(dat, -batch, names_to = "Feature", values_to = "Expression")
  g <- ggplot(dat, aes(Feature, Expression, group=Feature)) +
    geom_boxplot()
  if (by_batch) g <- g + facet_wrap(. ~ batch)
  print(g)
  invisible(g)
}
