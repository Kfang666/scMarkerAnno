# scMarkerAnno: Single-Cell RNA-seq Marker Gene Identification and Annotation

![R-CMD-check](https://github.com/Kfang666/scMarkerAnno/workflows/R-CMD-check/badge.svg)
[![codecov](https://codecov.io/gh/Kfang666/scMarkerAnno/branch/main/graph/badge.svg)](https://codecov.io/gh/Kfang666/scMarkerAnno)

## Overview

`scMarkerAnno` is an R package designed to streamline the analysis of single-cell RNA-sequencing (scRNA-seq) data. 
It provides a comprehensive workflow for data loading, quality control, dimensionality reduction, 
clustering, marker gene identification, and cell type annotation.

## Features

- Automated loading and merging of multiple 10X Genomics datasets
- Flexible sample grouping and batch correction
- Hierarchical clustering for fine-grained cell type identification
- Marker gene detection with customizable thresholds
- Cell type annotation based on predefined marker genes
- Comprehensive visualization capabilities

## Installation

You can install the development version of `scMarkerAnno` from GitHub:

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install scMarkerAnno
devtools::install_github("Kfang666/scMarkerAnno")