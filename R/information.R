#' Get order in which the samples were processed.
#'
#' @description
#' This function parses the complex measurement time column from MetIDQ to get
#'   samples processing order.
#'
#' @param eset A Biobase::ExpressionSet.
#' @param measurement_col The name of the column containing the measurement
#'   time (àefaults to "Measurement Time").
#' @param ref_index The index of the subcolumn to use for ordering of
#'   measurement times (typically 1, and therefore default value).
#' @return An integer vector of sample order.
#' @export
processing_order <- function(eset,
                             measurement_col = "Measurement Time",
                             ref_index = 1) {
  max_PTs <-
    max(sapply(strsplit(eset[[measurement_col]], split = " \\| "), length))
  op <- tidyr::separate(data.frame(MT=eset[[measurement_col]]),
                        "MT",
                        paste0("MT", 1:max_PTs),
                        sep = "\\|",
                        fill = "right")
  op <- as.data.frame(lapply(op, function(x) as.POSIXct(gsub("\\.", "-", x))))
  op <- apply(op, 1, function(x) {
    x <- sort(x,decreasing=T)
    length(x) <- 3
  x})
  return(rank(op[ref_index,]))
}

#' Get batch information from plate production number.
#'
#' @description
#' For large experiments with many samples, more than one IDQ Kit is often
#'   necessary. However, the results are often combined into one file. This
#'   function extracts the batch inforlation using the plate production number
#'   or any desired column. NAs are considered one batch.
#'
#' @param eset A Biobase::ExpressionSet.
#' @param batch_ind The name of the column containing the information to use for
#'   determining the batches. Defaults to "Plate Production No.).
#' @return An integer vector of batch numbers.
#' @export
batches <- function(eset, batch_ind = "Plate Production No.") {
  return(as.numeric(addNA(factor(eset[[batch_ind]]))))
}

# add annotations

#' Get the matrix of LOD values for each metabolite with respect to batches.
#'
#' @param eset A Biobase::ExpressionSet.
#' @return A matrix of LOD values of the same size as `exprs(eset)`.
#' @export
LOD_matrix <- function(eset) {
  pDat <- Biobase::pData(eset)
  pDat <- pDat[, base::startsWith(colnames(pDat), "PBC")]
  LOD_matrix <- apply(pDat, 1, function(x) {
    apply(Biobase::fData(eset)[, x[!is.na(x)], drop = F], 1, max)
  })
  return(LOD_matrix)
}

#' Get the proportions of values below LOD for each metabolite.
#'
#' @param eset A Biobase::ExpressionSet.
#' @return A numeric vector of the propotions of values below LOD.
#' @export
prop_below_LOD <- function(eset) {
  LOD_matrix <- LOD_matrix(eset)
  compMat <- Biobase::exprs(eset) < LOD_matrix
  res <- rowSums(compMat, na.rm = T) /
    (rep(ncol(compMat), nrow(compMat)) - rowSums(is.na(compMat)))
  return(res)
}
