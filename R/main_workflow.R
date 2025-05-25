# R/main_workflow.R
#' Main workflow for single-cell RNA-seq analysis
#'
#' @param data_dir Directory containing 10X Genomics data
#' @param sample_groups_file Path to a file defining sample groups (optional)
#' @param analysis_params List of analysis parameters
#' @param output_dir Directory to save results
#' @param custom_genes Custom gene list for visualization
#' @return A processed Seurat object
#' @export
main_workflow <- function(
    data_dir = "~/workspace",
    sample_groups_file = NULL,
    analysis_params = set_analysis_params(),
    output_dir = "./results",
    custom_genes = NULL
) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load and preprocess data
  pbmc <- load_and_preprocess_data(
    data_dir = data_dir,
    min_features = analysis_params$min_features,
    max_features = analysis_params$max_features,
    max_mito_percent = analysis_params$max_mito_percent
  )
  
  # Quality control
  pbmc <- perform_quality_control(pbmc)
  
  # Get available samples
  available_samples <- unique(pbmc$sample)
  
  # Read sample groups from file (if provided)
  if (!is.null(sample_groups_file)) {
    sample_groups <- read_sample_groups(sample_groups_file)
    
    # Validate groups
    sample_groups <- validate_sample_groups(sample_groups, available_samples)
    
    if (length(sample_groups) == 0) {
      stop("No valid sample groups found after validation.")
    }
  } else {
    sample_groups <- NULL
  }
  
  # Apply sample groups (if any)
  if (!is.null(sample_groups)) {
    pbmc <- assign_sample_groups(pbmc, sample_groups)
    group_by_var <- "sample_group"
    
    # Save group assignments
    group_df <- data.frame(
      cell_id = colnames(pbmc),
      sample = pbmc$sample,
      group = pbmc$sample_group
    )
    write.csv(group_df, file.path(output_dir, "sample_group_assignments.csv"), row.names = FALSE)
    
    # Output group statistics
    group_counts <- table(pbmc$sample_group)
    message("Sample group statistics:")
    print(group_counts)
  } else {
    group_by_var <- "sample"
    message("No sample groups provided. Using individual samples for analysis.")
  }
  
  # Dimensionality reduction and clustering
  pbmc <- perform_dimensionality_reduction(
    pbmc, 
    group_by = group_by_var,
    nfeatures = analysis_params$nfeatures,
    npcs = analysis_params$npcs,
    clustering_dims = analysis_params$clustering_dims,
    pca_analysis_dims = analysis_params$pca_analysis_dims
  )
  
  # Hierarchical clustering
  pbmc <- perform_hierarchical_clustering(
    pbmc,
    cell_threshold = analysis_params$cell_threshold,
    resolution = analysis_params$resolution_setting,
    parent_marker = analysis_params$parent_marker,
    pca_analysis_dims = analysis_params$pca_analysis_dims
  )
  
  # Find marker genes
  markers <- find_and_export_markers(
    pbmc,
    ident_col = "hierarchical_cluster",
    min.pct = analysis_params$min_pct,
    logfc.threshold = analysis_params$logfc_threshold,
    output_file = file.path(output_dir, "marker_genes.csv")
  )
  
  # Cell type annotation
  pbmc <- annotate_cell_types(pbmc)
  
  # Save final Seurat object
  saveRDS(pbmc, file.path(output_dir, "final_processed_data.rds"))
  
  # Generate visualizations with custom genes
  generate_visualizations(
    pbmc,
    output_dir = output_dir,
    group_by = c("sample_group", "hierarchical_cluster", "cell_type"),
    dims = analysis_params$clustering_dims,
    custom_genes = custom_genes,
    dotplot_group_by = "cell_type"
  )
  
  return(pbmc)
}