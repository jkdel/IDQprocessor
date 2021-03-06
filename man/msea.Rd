% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{msea}
\alias{msea}
\title{Perform metabolite set enrichment analysis (MSEA)}
\usage{
msea(
  metabolites,
  pathways,
  min.set.size = 5,
  fun = subramanian,
  ...,
  p.adj = "fdr"
)
}
\arguments{
\item{metabolites}{A `data.frame` with the first column containing metabolite
IDs and the second column the associated P-values.}

\item{pathways}{A `data.frame` with the first column containing a metabolite ID
and a second column containing a pathway ID, as obtained by `pathway_list()`.}

\item{min.set.size}{Minimal size of a set.}

\item{fun}{Function to call to get the enrichment statistic. Should take a
vector of P-values as first argument and the indices of the P-values
belonging to the set of interest as second argument. Further arguments can
be passed by `...`.}

\item{...}{Further arguments to pass to the function calculating the
enrichment statistic (e.g. number of permutations to perform).}

\item{p.adj}{Method for P-value adjustment. Default: FDR. `NULL` for none.}
}
\value{
A `data.frame` containing pathway ID, number of metabolites found in
  the pathway, P-value (significantly enriched) and adjusted P-value.
}
\description{
Perform metabolite set enrichment analysis (MSEA)
}
