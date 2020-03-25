#' Test for significance of enrichment using the Kolmogorov-Smirnow test.
#'
#' @param p All P-values, coefficients or similar values associated with the factor of interest .
#' @param set Indices of the values associated with the pathway of interest in p.
#' @param alternative See `?ks.test`. Typically, for p-value should be "less".
#'   For coefficents, or variable importance, where higher values are associated
#'   with more importance variables, it should be "greater".
#' @param nperm Number of permutations to get a P-value of pathway significance.
#'
#' @return A vector of p-values
#' @importFrom stats ks.test
#' @export
subramanian <- function(p, set, alternative="less", nperm = 1000) {
  ks.stat <- ks.test(p[set], p[-set], alternative = alternative, exact = F)$statistic
  m <- length(p)
  n <- length(set)
  ks.perm <- sapply(1:nperm, function(i) {
    setperm <- sample(1:m, n, replace = FALSE)
    ks.test(p[setperm], p[-setperm], alternative = alternative, exact = F)$statistic
  })
  pperm <- 1 - sum(ks.stat < ks.perm) / nperm
  return(pperm)
}

#' Perform metabolite set enrichment analysis (MSEA)
#'
#' @param metabolites A `data.frame` with the first column containing metabolite
#'   IDs and the second column the associated P-values, coefficients or similar values.
#' @param pathways A `data.frame` with the first column containing a metabolite ID
#'   and a second column containing a pathway ID, as obtained by `pathway_list()`.
#' @param min.set.size Minimal size of a set.
#' @param fun Function to call to get the enrichment statistic. Should take a
#'   numeric vector of P-values or similar (e.g. VIP score,
#'   absolute coefficients) as first argument and the indices of the values
#'   belonging to the set of interest as second argument. Further arguments can
#'   be passed by `...`.
#' @param ... Further arguments to pass to the function calculating the
#'   enrichment statistic (e.g. number of permutations to perform).
#' @param p.adj Method for P-value adjustment. Default: FDR. `NULL` for none.
#'
#' @return A `data.frame` containing pathway ID, number of metabolites found in
#'   the pathway, P-value (significantly enriched) and adjusted P-value.
#' @export
msea <- function(metabolites, pathways, min.set.size=5, fun=subramanian, ..., p.adj="fdr") {
  res <- pathways[pathways[,1] %in% metabolites[,1],-1]
  res <- res[!duplicated(res),]
  funarg <- list(...)
  res <- cbind(res, furrr::future_pmap_dfr(unname(res), function(pid, ...) {
    indices <- which(metabolites[,1] %in% pathways[pathways[,2]==pid,1])
    set_size <- length(indices)
    if (set_size>=min.set.size) {
      p <- do.call(fun, c(list(metabolites[,2], indices), funarg))
    } else {p <- NA}
    data.frame(set.size=set_size, p.value=p)
  }))
  res <- res[res$set.size>=min.set.size,]
  if (!is.null(p.adj)) res$p.adj <- stats::p.adjust(res$p.value, method = p.adj)
  res[order(res$p.value),]
}

#' @title Provides a list of pathways associated to each feature by tidying the
#'   annotations of an eset.
#'
#' @param met_df A `data.frame` of feature informations, typically `fData(eset)`.
#' @param met_var The name of the variable containing the feature IDs. If equal
#'   to "rownames", the row names are used as IDs.
#' @param path_var The name of the variable containing the pathway IDs.
#' @param path_desc The name of the variable containing the pathway descriptors
#'   (optional).
#' @param path_sep Character to use to split up the different pathways and
#'   descriptions associated with each feature.
#'
#' @return A `data.frame` with one row per pathway associated with each feature.
#' @export
pathway_list <- function(met_df, met_var, path_var, path_desc=NULL, path_sep="; ") {
  if (met_var=="rownames") met_df[[met_var]] <- rownames(met_df)
  pths <- strsplit(met_df[[path_var]], path_sep, fixed=T)
  df <- data.frame(Metabolite=rep(met_df[[met_var]],times=sapply(pths, length)),
                   Pathway=unlist(pths), stringsAsFactors = F)
  colnames(df) <- c(met_var, path_var)
  if (!is.null(path_desc))
    df[[path_desc]] <- unlist(strsplit(met_df[[path_desc]], path_sep, fixed=T))
  df <- df[!is.na(df[[path_var]]),]
  message(length(unique(df[[met_var]]))," out of ", nrow(met_df), " features were associated with ", length(unique(df[[path_var]])), " distinct pathways.")
  df
}
