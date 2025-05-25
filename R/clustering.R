# R/clustering.R
#' Perform hierarchical clustering on a Seurat object
#'
#' @param object Seurat object
#' @param cell_threshold Minimum number of cells to form a cluster
#' @param resolution Clustering resolution
#' @param parent_marker Separator for hierarchical cluster names
#' @param pca_analysis_dims Dimensions to use for PCA-based analysis
#' @return Seurat object with hierarchical clustering results
#' @export
perform_hierarchical_clustering <- function(
    object,
    cell_threshold = 2000,
    resolution = 1,
    parent_marker = "_",
    pca_analysis_dims = 1:30
) {
  # 初始聚类
  object <- FindNeighbors(object, dims = pca_analysis_dims)
  object <- FindClusters(object, resolution = resolution)
  
  # 初始化层次聚类结果
  object$hierarchical_cluster <- paste0("C", object$seurat_clusters)
  
  # 对每个聚类进行递归细分
  unique_clusters <- sort(unique(object$seurat_clusters))
  
  for (cluster in unique_clusters) {
    # 选择当前聚类的细胞
    cluster_cells <- Cells(object)[object$seurat_clusters == cluster]
    
    # 检查细胞数量是否超过阈值
    if (length(cluster_cells) >= cell_threshold) {
      # 创建子聚类
      sub_cluster_obj <- subset(object, cells = cluster_cells)
      
      # 重新进行降维和聚类
      sub_cluster_obj <- NormalizeData(sub_cluster_obj)
      sub_cluster_obj <- FindVariableFeatures(sub_cluster_obj)
      sub_cluster_obj <- ScaleData(sub_cluster_obj)
      sub_cluster_obj <- RunPCA(sub_cluster_obj, npcs = min(50, ncol(sub_cluster_obj)-1))
      
      # 使用指定的维度进行子聚类
      sub_cluster_obj <- FindNeighbors(sub_cluster_obj, dims = pca_analysis_dims)
      sub_cluster_obj <- FindClusters(sub_cluster_obj, resolution = resolution)
      
      # 更新主对象中的层次聚类标签
      for (sub_cluster in unique(sub_cluster_obj$seurat_clusters)) {
        sub_cells <- Cells(sub_cluster_obj)[sub_cluster_obj$seurat_clusters == sub_cluster]
        object$hierarchical_cluster[sub_cells] <- paste0(
          object$hierarchical_cluster[sub_cells[1]], 
          parent_marker, 
          sub_cluster
        )
      }
    }
  }
  
  return(object)
}