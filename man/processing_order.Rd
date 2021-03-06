% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/information.R
\name{processing_order}
\alias{processing_order}
\title{Get order in which the samples were processed.}
\usage{
processing_order(eset, measurement_col = "Measurement Time", ref_index = 1)
}
\arguments{
\item{eset}{A Biobase::ExpressionSet.}

\item{measurement_col}{The name of the column containing the measurement
time (defaults to "Measurement Time").}

\item{ref_index}{The index of the subcolumn to use for ordering of
measurement times (typically 1, and therefore default value).}
}
\value{
An integer vector of sample order.
}
\description{
This function parses the complex measurement time column from MetIDQ to get
  samples processing order.
}
