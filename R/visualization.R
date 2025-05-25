# R/visualization.R
#' Generate visualizations for a processed Seurat object
#'
#' @param object Seurat object
#' @param output_dir Directory to save visualizations
#' @param group_by Variables to group cells by for visualization
#' @param dims Dimensions to use for PCA and UMAP plots
#' @param custom_genes Custom gene list for dot plot visualization
#' @param dotplot_group_by Variable to group cells by for dot plot
#' @param dotplot_scale Scale factor for dot size in dot plot
#' @return NULL (saves plots to file)
#' @export
generate_visualizations <- function(
    object,
    output_dir = "./results",
    group_by = c("sample", "seurat_clusters", "cell_type"),
    dims = 1:30,
    custom_genes = NULL,
    dotplot_group_by = "cell_type",
    dotplot_scale = 8
) {
  # Create visualization directory
  vis_dir <- file.path(output_dir, "visualizations")
  if (!dir.exists(vis_dir)) {
    dir.create(vis_dir, recursive = TRUE)
  }
  
  # PCA plots
  for (group in group_by) {
    if (group %in% colnames(object@meta.data)) {
      pca_plot <- DimPlot(object, reduction = "pca", group.by = group, dims = dims[1:2])
      ggsave(file.path(vis_dir, paste0("pca_", group, ".png")), 
             pca_plot, width = 10, height = 8)
    }
  }
  
  # UMAP plots
  for (group in group_by) {
    if (group %in% colnames(object@meta.data)) {
      umap_plot <- DimPlot(object, reduction = "umap", group.by = group)
      ggsave(file.path(vis_dir, paste0("umap_", group, ".png")), 
             umap_plot, width = 10, height = 8)
    }
  }
  
  # Feature plots for top markers
  top_markers <- head(rownames(object@assays$RNA@var.features), 10)
  for (marker in top_markers) {
    if (marker %in% rownames(object)) {
      feature_plot <- FeaturePlot(object, features = marker)
      ggsave(file.path(vis_dir, paste0("feature_", marker, ".png")), 
             feature_plot, width = 10, height = 8)
    }
  }
  
  # Heatmap
  heatmap_plot <- DoHeatmap(object, features = top_markers, group.by = group_by[1])
  ggsave(file.path(vis_dir, "heatmap_top_markers.png"), 
         heatmap_plot, width = 12, height = 10)
  
  # Dot plot for custom genes (if provided)
  if (!is.null(custom_genes) && length(custom_genes) > 0) {
    # Filter genes that exist in the dataset
    valid_genes <- custom_genes[custom_genes %in% rownames(object)]
    
    if (length(valid_genes) > 0) {
      # Check if group variable exists
      if (dotplot_group_by %in% colnames(object@meta.data)) {
        dot_plot <- DotPlot(
          object, 
          features = valid_genes,
          group.by = dotplot_group_by,
          dot.scale = dotplot_scale,
          scale = TRUE
        ) +
          scale_color_gradientn(colors = c("#4575B4", "#FFFFBF", "#D73027")) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(file.path(vis_dir, "dotplot_custom_genes.png"), 
               dot_plot, width = 12, height = 10)
        
        message(paste("Dot plot saved for", length(valid_genes), "custom genes"))
      } else {
        warning(paste("Group variable", dotplot_group_by, "not found in metadata"))
      }
    } else {
      warning("None of the provided custom genes were found in the dataset")
    }
  }
  
  message(paste("Visualizations saved to:", vis_dir))
}