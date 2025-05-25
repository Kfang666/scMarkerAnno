# R/config.R
#' Set analysis parameters
#'
#' @param cell_threshold Minimum cell count for reclustering
#' @param pca_dims Dimensions to use for PCA calculation
#' @param resolution_setting Clustering resolution
#' @param parent_marker Separator for hierarchical clustering
#' @param min_pct Minimum percentage of cells expressing a gene
#' @param logfc_threshold Log fold change threshold for marker genes
#' @param umap_dims Dimensions to use for UMAP
#' @param min_features Minimum number of features per cell
#' @param max_features Maximum number of features per cell
#' @param max_mito_percent Maximum mitochondrial percentage per cell
#' @param nfeatures Number of highly variable features to retain
#' @param npcs Number of principal components to compute
#' @param clustering_dims Dimensions to use for clustering
#' @param pca_analysis_dims Dimensions to use for PCA-based analysis
#' @return A list of analysis parameters
#' @export
set_analysis_params <- function(
    cell_threshold = 2000,
    pca_dims = 1:30,
    resolution_setting = 1,
    parent_marker = "_",
    min_pct = 0.25,
    logfc_threshold = 0.585,
    umap_dims = 1:30,
    min_features = 200,
    max_features = 4000,
    max_mito_percent = 15,
    nfeatures = 2000,
    npcs = 50,
    clustering_dims = 1:30,
    pca_analysis_dims = 1:30
) {
  params <- list(
    cell_threshold = cell_threshold,
    pca_dims = pca_dims,
    resolution_setting = resolution_setting,
    parent_marker = parent_marker,
    min_pct = min_pct,
    logfc_threshold = logfc_threshold,
    umap_dims = umap_dims,
    min_features = min_features,
    max_features = max_features,
    max_mito_percent = max_mito_percent,
    nfeatures = nfeatures,
    npcs = npcs,
    clustering_dims = clustering_dims,
    pca_analysis_dims = pca_analysis_dims
  )
  
  # Validate parameters
  if (cell_threshold <= 0) stop("cell_threshold must be a positive integer")
  if (resolution_setting < 0) stop("resolution_setting must be non-negative")
  if (min_pct < 0 || min_pct > 1) stop("min_pct must be between 0 and 1")
  if (min_features <= 0 || max_features <= 0 || min_features >= max_features) {
    stop("Invalid min_features or max_features values")
  }
  if (max_mito_percent < 0 || max_mito_percent > 100) {
    stop("max_mito_percent must be between 0 and 100")
  }
  if (nfeatures <= 0) stop("nfeatures must be a positive integer")
  if (npcs <= 0) stop("npcs must be a positive integer")
  if (max(clustering_dims) > npcs) {
    stop("clustering_dims cannot include dimensions greater than npcs")
  }
  if (max(pca_analysis_dims) > npcs) {
    stop("pca_analysis_dims cannot include dimensions greater than npcs")
  }
  
  return(params)
}