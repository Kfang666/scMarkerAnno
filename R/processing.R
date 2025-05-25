# R/processing.R
#' Perform dimensionality reduction on a Seurat object
#'
#' @param object Seurat object
#' @param group_by Variable to group cells by
#' @param nfeatures Number of highly variable features to retain
#' @param npcs Number of principal components to compute
#' @param clustering_dims Dimensions to use for clustering and UMAP
#' @param pca_analysis_dims Dimensions to use for PCA-based analysis
#' @return Processed Seurat object
#' @export
perform_dimensionality_reduction <- function(
    object,
    group_by = "sample",
    nfeatures = 2000,
    npcs = 50,
    clustering_dims = 1:30,
    pca_analysis_dims = 1:30
) {
  # Normalize data
  object <- NormalizeData(object)
  
  # Find variable features
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures)
  
  # Scale data
  object <- ScaleData(object)
  
  # Perform PCA with specified number of components
  object <- RunPCA(object, npcs = npcs, verbose = FALSE)
  
  # Determine dimensionality for clustering
  if (is.null(clustering_dims)) {
    clustering_dims <- 1:min(30, npcs)
  }
  
  # Cluster cells using specified dimensions
  object <- FindNeighbors(object, dims = clustering_dims)
  object <- FindClusters(object, resolution = 0.5)
  
  # Run UMAP using specified dimensions
  object <- RunUMAP(object, dims = clustering_dims)
  
  return(object)
}