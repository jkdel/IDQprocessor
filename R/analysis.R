#' Test for significance of enrichment using the Kolmogorow-Smirnow Test.
#'
#' @param p All P-values associated with the factor of interest .
#' @param set Indices of the p values associated with the pathway of interest.
#' @param nperm Number of permutations to get a p-value by gene permutations.
#'
#' @return A vector of p-values
#' @export
subramanian <- function(p, set, nperm = 1000) {
  ks.stat <- ks.test(p[set], p[-set], alternative = "less", exact = F)$statistic
  m <- length(p)
  n <- length(set)
  ks.perm <- sapply(1:nperm, function(i) {
    setperm <- sample(1:m, n, replace = FALSE)
    ks.test(p[setperm], p[-setperm], alternative = "less", exact = F)$statistic
  })
  pperm <- 1 - sum(ks.stat < ks.perm) / nperm
  return(pperm)
}

#' Perform metabolite set enrichment analysis (MSEA)
#'
#' @param metabolites A `data.frame` with the first column containing metabolite
#'   IDs and the second column the associated P-values.
#' @param pathways A `data.frame` with the first column containing a pathway ID
#'   and a second column containing a metabolite ID.
#' @param p.adj Method for P-value adjustment. Default: FDR. `NULL` for none.
#' @param min.set.size Minimal size of a set, should be at least 2.
#' @param fun Function to use to test for enrichment. Default is `subramanian`
#'   to use the Kolmogorow-Smirnow test. You can provide any function taking a
#'   vector of p-values as first argument and the indices of the
#'   genes/metabolites of interest within these. Further arguments are allowed.
#' @param ... Further arguments that will be passed to `fun`.
#'
#' @return A `data.frame` containing pathway ID, number of metabolites found in
#'   the pathway, P-value (significantly enriched) and adjusted P-value.
#' @export
msea <- function(metabolites, pathways, p.adj="fdr", min.set.size=5, fun=subramanian, ...) {
  pathways <- pathways[pathways[,2] %in% metabolites[,1],]
  res <- base::do.call(base::rbind.data.frame,
                       future.apply::future_lapply(unique(pathways[,1]),
                                                   function(pathway_ID) {
                                                     indices <- base::which(metabolites[,1] %in% pathways[pathways[,1]==pathway_ID,2])
                                                     set_size <- base::length(indices)
                                                     if (set_size>=min.set.size) {
                                                       p <- fun(metabolites[,2], indices, ...)
                                                     } else {p <- NA}
                                                     c(pathway_ID, set_size, p)
                                                   }))
  colnames(res) <- c("Pathway","Set.size","P")
  if (!is.null(p.adj)) res$P.adj <- stats::p.adjust(res$P, method = p.adj)
  res[order(res$P),]
}
