#' Aggregate ExpressionSet values according to an indicator column
#'
#' @param eset A `Biobase::ExpressionSet`.
#' @param ind_dupes The name of the column of pData(eset) to use to find
#'   duplicate values to aggregate.
#' @param fun Aggregatation function. Should be able to handle both numeric and
#'   character vectors. Aggregation is performed sample wise for both the
#'   expression and the phenotype data, i.e. respectively column and row-wise.
#' @return A `Biobase::ExpressionSet`.
#' @export
aggregate_eset <- function(eset, ind_dupes, fun) {
  fd <- stats::aggregate(t(Biobase::exprs(eset)), by = list(eset[[ind_dupes]]), fun)
  sd <- stats::aggregate(Biobase::pData(eset), by = list(eset[[ind_dupes]]), fun)
  return(Biobase::ExpressionSet(
    assayData = t(fd[,-1]),
    phenoData = methods::as(sd[,-1], "AnnotatedDataFrame"),
    featureData = Biobase::featureData(eset),
    experimentData = Biobase::experimentData(eset)
  ))
}
