#' Transform an excel file from MetIDQ into a Biobase::ExpressionSet
#'
#' @param path A string saying where to find the excel file to use as source.
#' @param col_before String giving the name of the last column of metadata.
#'   Typically it is called 'Measurement.Time'.
#' @param LOD_string String used to parse the LOD rownames.
#'   Defaults to "LOD \\(calc\\.\\) ([0-9]+)/([0-9]+) \\[µM\\]".
#' @param LOD_rep String to rename the LOD columns.
#' @param plate_code Name of the column containing the Plate car code,
#'   typically "Plate.Bar.Code".
#' @return A Biobase::ExpressionSet object.
#' @export
IDQ2eset <- function(path,
                     col_before = "Measurement Time",
                     LOD_string = "LOD \\(calc\\.\\) ([0-9]+)/([0-9]+) \\[\u00B5M\\]",
                     LOD_rep = "\\1-\\2",
                     plate_code = "Plate Bar Code") {
  # First row contains an uninteresting string
  op <- openxlsx::read.xlsx(path, startRow = 2, sep.names = " ")
  col_start <- which(colnames(op) == col_before)
  row_start <- max(which(grepl(LOD_string, op[, col_before])))
  fd <- op[-1:-row_start, -1:-col_start]
  fd <- as.data.frame(t(sapply(fd, as.numeric)))
  sd <- op[-1:-row_start, 1:(col_start)]
  if (plate_code %in% colnames(sd)) {
    max_BCs <- max(sapply(strsplit(sd[[plate_code]], split = " \\| "), length))
    sd <- tidyr::separate(sd,
                          plate_code,
                          paste0("PBC", 1:max_BCs),
                          sep = " \\| ",
                          fill = "right")
  }
  # appears unnecessary complex, but t.data.frame does strange things here
  md <- utils::type.convert(as.data.frame(t(op[1:row_start, -1:-col_start]),
                                          stringsAsFactors = F))
  colnames(md) <- sub(LOD_string, LOD_rep, op[1:row_start,col_start])
  colnames(fd) <- rownames(sd)
  eset <- Biobase::ExpressionSet(
    assayData = as.matrix(fd),
    phenoData = methods::as(sd, "AnnotatedDataFrame"),
    featureData = methods::as(md, "AnnotatedDataFrame")
  )
  return(eset)
}
