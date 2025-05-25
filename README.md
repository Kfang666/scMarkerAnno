# scMarkerAnno: Single-Cell RNA-seq Marker Gene Identification and Annotation

![R-CMD-check](https://github.com/Kfang666/scMarkerAnno/workflows/R-CMD-check/badge.svg)(https://codecov.io/gh/Kfang666/scMarkerAnno)

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

## User Guide
# 加载包
library(scRecluster)
library(Seurat)
library(dplyr)
# 质控
params <- set_analysis_params(
  cell_threshold = 1500,        # 触发重聚类的最小细胞数
  resolution_setting = 0.7,     # 聚类分辨率
  min_pct = 0.25,               # 标记基因筛选的最小表达比例
  logfc_threshold = 0.5,        # 标记基因筛选的最小log2FC
  min_features = 200,           # 最小基因数过滤
  max_features = 6000,          # 最大基因数过滤
  max_mito_percent = 15,        # 最大线粒体基因比例
  nfeatures = 2500,             # 高变基因数量
  npcs = 40,                    # PCA主成分数量
  pca_analysis_dims = 1:30,     # PCA分析使用的维度
  clustering_dims = 1:25        # 聚类使用的维度
)
# 定义自定义基因列表（用于可视化）
custom_genes <- c(
  "CD3E", "CD4", "CD8A",        # T细胞标记
  "CD19", "MS4A1",              # B细胞标记
  "CD14", "LYZ",                # 单核细胞标记
  "CD56", "NCAM1",              # NK细胞标记
  "FOXP3", "IL2RA",             # 调节性T细胞标记
  "PDCD1", "CTLA4"              # 免疫检查点基因
)
# 运行完整分析流程（含注释等待）
results <- main_workflow(
  data_dir = "~/workspace/10X_data",           # 数据目录
  sample_groups_file = "~/sample_groups.txt",  # 样本分组文件
  analysis_params = params,                     # 分析参数
  output_dir = "~/scRecluster_results",         # 输出目录
  custom_genes = custom_genes,                  # 自定义基因列表
  wait_for_annotation = TRUE,                   # 等待手动注释
  dims_setting = list(
    pca_dims = 1:30,            # PCA分析维度
    clustering_dims = 1:25,     # 聚类维度
    umap_dims = 1:25            # UMAP维度
  )
)
# 保存最终结果
saveRDS(results, file.path("~/scRecluster_results", "final_annotated_seurat.rds"))
