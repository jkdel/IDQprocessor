#' Plot the assay plate layout.
#'
#' @param eset A Biobase::ExpressionSet.
#' @param highlight An optional expression used to highlight certain samples.
#' @param fill Name of the column used to distinguish samples and wells.
#' @param well_col Name of the column containing the well position.
#' @param batch_ind The name of the column containing the information to use for
#'   determining the batches. Defaults to "Plate Production No.).
#' @return An invisible `ggplot()`.
#' @example
#' \notrun{
#' plot_layout(es, `Sample Identification` == "MMU08_EDTA_96dpi")
#' }
#' @import ggplot2
#' @export
plot_layout <- function(eset,
                        highlight,
                        fill = "Sample Type",
                        well_col = "Well Position",
                        batch_ind = "Plate Production No.") {
  pDat <- Biobase::pData(eset)
  pDat$PCol <- factor(((pDat[[well_col]]-1) %% 12)+1)
  pDat$PRow <- factor(((pDat[[well_col]]-1) %/% 12)+1)
  pDat$PRow <- factor(pDat[["PRow"]], levels = rev(levels(pDat[["PRow"]])))
  pDat$batch <- batches(eset, batch_ind)
  if (!missing(highlight)) {
    pDat$high <- ifelse(eval(substitute(highlight), pDat),"red",NA)
  }
  g <- ggplot(pDat, aes(PCol, PRow, fill = pDat[[fill]])) +
    scale_fill_discrete(fill) +
    scale_color_identity() +
    coord_equal() +
    facet_wrap(. ~ batch) +
    theme_void() +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_legend(ncol = 5)) +
    theme(legend.position = "top")
  if (missing(highlight)) print(g + geom_tile(size = 1.2, width = 0.9, height = 0.9))
  else print(g + geom_tile(size = 1.2, width = 0.9, height = 0.9, color = pDat[["high"]]))
  invisible(g)
}
