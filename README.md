# Introduction

When using Biocrates metabolomic panels, data is managed using the MetIDQ software. For data analysis in R is used excel files exported from MetIDQ without normalization.

This package provides:
* an import function that parses said excel file into a `Biobase::ExpressionSet`
* functions for preprocessing and visualization of omics data, that work with `ExpressionSet`s such as:
  * `plot_layout` to plot the arrangement of sample on the plate
  * `plot_pca` to plot a PCA with one grouping factor
  * `plot_qcs` to plot the expression levels of selected metabolites over processing time (to investigate any drift during sample measurement) while connecting the measurements of QC-samples
  * `plot_features` to produce boxplots of each feature
  * `aggregate_eset` to easily merge two `ExpressionSet`s
  * `qc_rlsc` to perform QC-RLSC normalization on a selected QC-sample
  * and many more...

# Details

You can learn more about the `Biobase::ExpressionSet` [here](https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf). This data structure, initially developed for transcriptomic assays, ensures that the measured feature expression (e.g. in this case metabolite concentrations) stays with the sample data and any additional information like the LOD for each metabolite and batch. Filtering an `ExpressionSet` by sample data is easy and prevents getting mixed up between several tables, potentially not ordered accurately.

# Limitations

The package was written for my personal use and is still work in progress (bugs are to be expected)! It is provided under a [CC-BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) licence.

# Installation

This R package can be installed via github:
```
devtools::install_github("jkdel/IDQprocessor")
```
