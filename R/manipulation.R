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

#' Add summary features to a `Biobase::ExpressionSet`
#'
#' @param eset A `Biobase::ExpressionSet`.
#' @param group_ind The name of the column of dData(eset) to use to find
#'   values to aggregate.
#' @param fun Aggregatation function to apply to each column. By default `sum`.
#' @return A `Biobase::ExpressionSet`.
#' @export
add_summaries <- function(eset, group_ind, fun=sum) {
  groups <- base::split(base::rownames(Biobase::fData(eset)),Biobase::fData(eset)[[group_ind]])
  groups <- groups[sapply(groups,length)>1]
  summaries <- t(base::sapply(groups, function(x) apply(Biobase::exprs(eset[x,]),2,fun)))
  fd <- Biobase::fData(eset)
  dd <- as.data.frame(base::matrix(nrow = length(groups),ncol = ncol(fd)))
  colnames(dd) <- colnames(fd)
  rownames(dd) <- names(groups)
  dd[,group_ind] <- names(groups)
  fd <- base::rbind(fd, dd)
  return(Biobase::ExpressionSet(
    assayData = base::rbind(Biobase::exprs(eset), summaries),
    phenoData = Biobase::phenoData(eset),
    featureData = methods::as(fd, "AnnotatedDataFrame"),
    experimentData = Biobase::experimentData(eset)
  ))
}
