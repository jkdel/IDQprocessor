% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_features}
\alias{plot_features}
\title{Boxplot of each feature.}
\usage{
plot_features(eset, by_batch = T, batch_ind = "Plate Production No.", plot = T)
}
\arguments{
\item{eset}{A Biobase::ExpressionSet.}

\item{by_batch}{Logical, `TRUE` if separate facets should be used for each
bacth.}

\item{batch_ind}{Name of the column from where to derive batch information.}

\item{plot}{Logical. Should the plot be plotted right away?}
}
\value{
An invisible `ggplot()`.
}
\description{
Boxplot of each feature.
}
