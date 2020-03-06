#' Annotations for Biocrates IDQ p180.
#'
#' @description Annotation informations based on the SMPDB. Some metabolites
#'   determined by FIA correspond to multiple isomeres. In these cases either
#'   the first best match in the HMDB database was taken, or no database ID was
#'   provided. Last updated 2020-03-06.
#'
#' @example
#' \notrun{
#' library(tidyverse)
#' fData(mii) <- fData(mii) %>%
#'   mutate(Descriptor = rownames(.)) %>%
#'   left_join(p180_annot) %>%
#'   column_to_rownames("Descriptor")
#' }
#' @source \url{http://www.smpdb.ca/}
"p180_annot"
