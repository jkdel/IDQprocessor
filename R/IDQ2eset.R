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
#' @param experimentData Optional Biobase::MIAME object.
#' @param try_annotate If `TRUE`, tries to annotate the metabolites using a
#'   custom annotation file.
#' @return A Biobase::ExpressionSet object.
#' @examples
#' \dontrun{
#' IDQ2eset("orig/TiHo Juni2019_nn_nc.xlsx",
#'          experimentData = MIAME(name = "Delarocque, Julien",
#'                                 lab = "Clinic for Horse, University of Veterinary Medicine Hannover, Foundation",
#'                                 contact = "julien.delarocque@tiho-hannover.de",
#'                                 title = "Equine Metabolome Seasonality in OGT",
#'                                 abstract = "19 horses were subjected to 5 OGTs at regular intervals over a one year period",
#'                                 samples = list(all = "Basal and 120 min EDTA Plasma samples"),
#'                                 normControls = list(intQC = "Own QC called 02-W11-OGT-P-30")))
#' }
#' @export
IDQ2eset <- function(path,
                     col_before = "Measurement.Time",
                     LOD_string = "LOD \\(calc\\.\\) ([0-9]+)/([0-9]+) \\[\u03BCM\\]",
                     LOD_rep = "\\1 - \\2",
                     plate_code = "Plate.Bar.Code",
                     experimentData = Biobase::MIAME(),
                     try_annotate = T) {
  # First row contains an uninteresting string
  op <- openxlsx::read.xlsx(path, startRow = 2)
  col_start <- which(colnames(op) == col_before)
  row_start <- max(which(grepl(LOD_string, op[, col_before])))
  fd <- op[-1:-row_start, -1:-col_start]
  fd <- as.data.frame(t(sapply(fd, as.numeric)))
  sd <- op[-1:-row_start, 1:(col_start - 1)]
  if (plate_code %in% colnames(sd)) {
    max_BCs <- max(sapply(strsplit(sd[[plate_code]], split = " \\| "), length))
    sd <- tidyr::separate(sd,
                          plate_code,
                          paste0("PBC", 1:max_BCs),
                          sep = "\\|",
                          fill = "right")
  }
  # appears unnecessary complex, but t.data.frame does strange things here
  md <- type.convert(as.data.frame(t(op[1:row_start, -1:-(col_start)]),
                                   stringsAsFactors = F))
  colnames(md) <- sub(LOD_string, LOD_rep, op[1:row_start,col_start])
  colnames(fd) <- rownames(sd)
  annotated <- F
  if (try_annotate & nrow(md)) {
    md <- merge(md,an,by.x="row.names",by.y="descriptor")
    rownames(md) <- rownames(fd)
    annotated <- T
  }
  eset <- Biobase::ExpressionSet(
    assayData = as.matrix(fd),
    phenoData = as(sd, "AnnotatedDataFrame"),
    featureData = as(md, "AnnotatedDataFrame"),
    experimentData = experimentData,
    annotation = ifelse(annotated,"manual_p180","")
  )
  return(eset)
}
