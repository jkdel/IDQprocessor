#' Log transfrom the data.
#'
#' This log transformation is robust again negative, zero and missing values.
#'
#' @param eset A Biobase::ExpressionSet.
#' @param log_base The log base to use. Typically 2.
#' @return A `Biobase::ExpressionSet()`.
#' @export
log_transform <- function(eset, log_base=2) {
  dat <- t(Biobase::exprs(eset))
  mins <- apply(dat, 2, function(x) min(abs(x[x!=0]),na.rm=T)/10)
  dat <- sweep(dat, 2, mins, function(y, min.val) {
    log((y + sqrt(y^2 + min.val^2))/2, base=log_base)
  })
  Biobase::exprs(eset) <- t(dat)
  return(eset)
}

# #' Perform quantile normalization on an ExpressionSet.
# #'
# #' @param eset A Biobase::ExpressionSet.
# #' @return A `Biobase::ExpressionSet()`.
# #' @export
# quantile_norm <- function(eset) {
#   exprs(eset) <- preprocessCore::normalize.quantiles(exprs(eset), F)
#   return(eset)
# }

#' Perform QC-RLSC for batch correction of chromatographic signal
#'
#' @param eset A `Biobase::ExpressionSet`.
#' @param qc_ind Either a vector corresponding to the `eset` samples, where `1`
#'   indicates a QC and `2` indicates a samples. Or the name of the column from
#'   where to derive this information. Repeatedly measured samples will be used
#'   as QC. If more than one sample was repeatedly mesured, only the first
#'   repeated sample identified will be used.
#' @param degree The degree of polynomials to be used.
#' @param span The parameter \ifelse{latex}{\out{$\alpha$}}{\ifelse{html}{\out{&alpha;}}{alpha}} which controls the LOESS smoothing.
#' @param measurement_col The name of the column from where to derive processing
#'   order (either using the MetIDQ Mesurement Time or by a vector of integers).
#' @param ref_index The index of the subcolumn to use for ordering of
#'   measurement times (typically 1, and therefore default value). Will only be
#'   used if the column corresponding to `measurement_column` is not a vector of
#'   integers.
#' @return A `Biobase::ExpressionSet`.
#' @export
#' @references Dunn et al. Nature Protocols 6, 1060-1083 (2011)
qc_rlsc <- function(eset, qc_ind = "Sample Identification", degree = 2,
                    span =  0.75, measurement_col = "Measurement Time",
                    ref_index = 1) {
  if (is.numeric(eset[[measurement_col]])) proc <- eset[[measurement_col]]
  else proc <- processing_order(eset, measurement_col, ref_index)
  eset <- eset[,order(proc)]
  or <- base::order(proc[base::order(proc)])
  tab <- as.data.frame(t(Biobase::exprs(eset)))
  if (is.numeric(eset[[qc_ind]])) {
    if (base::table(eset[[qc_ind]]==1)["TRUE"]<2) stop("There are either not enough QC samples or the `qc_ind` column is not accurately formated. Please refer to ?qc_rlsc.", call. = F)
    colv <- eset[[qc_ind]]
  } else {
    res <- base::which(table(eset[[qc_ind]])>2)
    if (length(res)>1) warning("There are more than one sample repeatedly measured. Will use the first one as reference QC.", call. = F)
    colv <- base::ifelse(eset[[qc_ind]]==names(res)[[1]], 1, 2)
  }
  if (or[base::which(colv==1)][[1]]!=1 | rev(or[base::which(colv==1)])[[1]]!=length(or))
    warning("The first and/or last sample is not a QC. Some expression levels might not be adjusted and set to NA!", call. = F)
  if (length(base::which(colv==1)*10)<length(or))
    warning("It seems that there are only very few QCs, normalization might not be accurate!", call. = F)
  tab <- base::sapply(tab, function(x) {
    x1 <- x[!is.na(x)]
    or1 <- or[!is.na(x)]
    colv1 <- colv[!is.na(x)]
    lo <- stats::loess(x1[which(colv1 == 1)] ~ or1[which(colv1 == 1)],
                       span = span,
                       degree = degree,
                       control = stats::loess.control(surface = "direct"))
		ap <- stats::approx(x = or1[which(colv1 == 1)], y = lo$fitted, xout = or)
		return(x / ap$y)
  })
  Biobase::exprs(eset) <- t(tab)
  return(eset)
}
