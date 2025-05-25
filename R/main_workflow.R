# R/main_workflow.R
#' Main workflow for single-cell RNA-seq analysis
#'
#' @param data_dir Directory containing 10X Genomics data
#' @param sample_groups_file Path to a file defining sample groups (optional)
#' @param analysis_params List of analysis parameters
#' @param output_dir Directory to save results
#' @param custom_genes Custom gene list for visualization
#' @param annotation_file Path to the cell type annotation file (optional)
#' @param wait_for_annotation Whether to wait for user to create annotation file
#' @param dims_setting Dimensions to use for PCA, clustering, and UMAP (list with pca_dims, clustering_dims, umap_dims)
#' @return A processed Seurat object
#' @export
main_workflow <- function(
    data_dir = "~/workspace",
    sample_groups_file = NULL,
    analysis_params = set_analysis_params(),
    output_dir = "./results",
    custom_genes = NULL,
    annotation_file = NULL,
    wait_for_annotation = TRUE,
    dims_setting = list(pca_dims = 1:30, clustering_dims = 1:30, umap_dims = 1:30)
) {
  # 创建输出目录
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # ----------------------
  # 1. 数据加载与预处理
  # ----------------------
  pbmc <- load_and_preprocess_data(
    data_dir = data_dir,
    min_features = analysis_params$min_features,
    max_features = analysis_params$max_features,
    max_mito_percent = analysis_params$max_mito_percent
  )
  
  # ----------------------
  # 2. 样本分组与质量控制
  # ----------------------
  pbmc <- perform_quality_control(pbmc)
  sample_groups <- handle_sample_groups(sample_groups_file, pbmc)
  
  # ----------------------
  # 3. 降维与批次校正
  # ----------------------
  pbmc <- run_dim_reduction(
    pbmc,
    group_by = sample_groups$group_by_var,
    nfeatures = analysis_params$nfeatures,
    npcs = analysis_params$npcs,
    dims_setting = dims_setting
  )
  
  # ----------------------
  # 4. 层次聚类与标记基因分析
  # ----------------------
  pbmc <- perform_hierarchical_clustering(
    pbmc,
    cell_threshold = analysis_params$cell_threshold,
    resolution = analysis_params$resolution_setting,
    pca_dims = dims_setting$pca_dims
  )
  
  markers_file <- export_markers(pbmc, output_dir)
  
  # ----------------------
  # 5. 细胞类型注释（手动流程）
  # ----------------------
  if (wait_for_annotation) {
    annotation_file <- prompt_for_annotation(markers_file, output_dir)
  }
  
  if (!is.null(annotation_file)) {
    pbmc <- annotate_cell_types(pbmc, annotation_file)
  } else {
    pbmc$cell_type_OK <- "Unannotated"
  }
  
  # ----------------------
  # 6. 结果保存与可视化
  # ----------------------
  save_final_objects(pbmc, output_dir)
  generate_visualizations(pbmc, output_dir, custom_genes, dims_setting$umap_dims)
  
  message("\n=== 分析完成 ===")
  message(paste("结果已保存至:", output_dir))
  return(pbmc)
}


### ----------------------
### 辅助函数：数据加载与预处理
### ----------------------
load_and_preprocess_data <- function(data_dir, min_features, max_features, max_mito_percent) {
  # 设置工作目录（用户自定义路径）
  original_wd <- setwd(data_dir)
  on.exit(setwd(original_wd))  # 恢复原始工作目录
  
  # 遍历子目录读取10X数据
  sub_dirs <- list.dirs(full.names = TRUE, recursive = FALSE)
  seurat_objects <- list()
  
  for (dir in sub_dirs) {
    if (check_10x_files(dir)) {
      sample_name <- basename(dir)
      seurat_obj <- create_seurat_object(dir, sample_name, min_features)
      seurat_objects[[sample_name]] <- seurat_obj
      message(paste("已处理样本:", sample_name))
    }
  }
  
  if (length(seurat_objects) == 0) {
    stop("未找到有效样本数据")
  }
  
  # 合并样本
  pbmc <- merge_seurat_objects(seurat_objects)
  
  # 质量控制：线粒体比例计算与过滤
  pbmc <- calculate_mito_ratio(pbmc)
  pbmc <- filter_cells(pbmc, min_features, max_features, max_mito_percent)
  
  # 标准化与特征选择
  pbmc <- normalize_data(pbmc)
  pbmc <- select_variable_features(pbmc)
  pbmc <- scale_data(pbmc)
  
  return(pbmc)
}

check_10x_files <- function(dir) {
  required_files <- c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
  all(file.exists(file.path(dir, required_files)))
}

create_seurat_object <- function(dir, sample_name, min_features) {
  counts <- Read10X(data.dir = dir)
  CreateSeuratObject(
    counts = counts,
    project = paste0("Project_", sample_name),
    min.cells = 3,
    min.features = min_features
  ) %>% 
    AddMetaData(sample = sample_name)
}

merge_seurat_objects <- function(seurat_objects) {
  Reduce(
    function(x, y) merge(x, y, add.cell.ids = TRUE),
    seurat_objects
  ) %>% 
    JoinLayers()  # 合并多组学数据层（如有）
}

calculate_mito_ratio <- function(pbmc) {
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc
}

filter_cells <- function(pbmc, min_features, max_features, max_mito_percent) {
  subset(pbmc, 
         subset = nFeature_RNA > min_features & 
           nFeature_RNA < max_features & 
           percent.mt < max_mito_percent)
}

normalize_data <- function(pbmc) {
  NormalizeData(pbmc)
}

select_variable_features <- function(pbmc, nfeatures = 2000) {
  FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nfeatures)
}

scale_data <- function(pbmc) {
  ScaleData(pbmc)
}


### ----------------------
### 辅助函数：样本分组处理
### ----------------------
handle_sample_groups <- function(sample_groups_file, pbmc) {
  available_samples <- unique(pbmc$sample)
  
  if (!is.null(sample_groups_file)) {
    sample_groups <- read_sample_groups(sample_groups_file)
    sample_groups <- validate_sample_groups(sample_groups, available_samples)
    
    if (length(sample_groups) == 0) {
      stop("验证后未找到有效样本分组")
    }
    
    pbmc <- assign_sample_groups(pbmc, sample_groups)
    group_by_var <- "sample_group"
    
    save_sample_groups(pbmc, sample_groups_file, output_dir)
    message("样本分组统计:")
    print(table(pbmc$sample_group))
  } else {
    group_by_var <- "sample"
    message("未提供样本分组，使用单个样本分析")
  }
  
  return(list(group_by_var = group_by_var))
}

read_sample_groups <- function(file) {
  read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
    deframe()  # 转换为命名向量（样本名: 分组名）
}

validate_sample_groups <- function(groups, available_samples) {
  groups[names(groups) %in% available_samples]  # 过滤无效样本
}

assign_sample_groups <- function(pbmc, groups) {
  pbmc$sample_group <- groups[pbmc$sample]
  pbmc
}

save_sample_groups <- function(pbmc, file, output_dir) {
  data.frame(
    cell_id = colnames(pbmc),
    sample = pbmc$sample,
    group = pbmc$sample_group
  ) %>% 
    write.csv(file.path(output_dir, "sample_group_assignments.csv"), row.names = FALSE)
}


### ----------------------
### 辅助函数：降维与聚类
### ----------------------
run_dim_reduction <- function(pbmc, group_by, nfeatures, npcs, dims_setting) {
  # PCA 分析
  pbmc <- RunPCA(
    pbmc,
    features = VariableFeatures(pbmc),
    npcs = npcs
  )
  
  # Harmony 批次校正（按样本分组）
  pbmc <- RunHarmony(pbmc, group.by.vars = group_by)
  
  # 邻居查找与聚类
  pbmc <- FindNeighbors(pbmc, dims = dims_setting$pca_dims)
  pbmc <- FindClusters(pbmc, resolution = analysis_params$resolution_setting)
  
  # UMAP 可视化
  pbmc <- RunUMAP(pbmc, dims = dims_setting$umap_dims)
  
  return(pbmc)
}


### ----------------------
### 辅助函数：标记基因与注释
### ----------------------
export_markers <- function(pbmc, output_dir) {
  markers_file <- file.path(output_dir, "filtered_markers_all_0.25.csv")
  find_and_export_markers(
    pbmc,
    ident_col = "hierarchical_cluster",
    min.pct = analysis_params$min_pct,
    logfc.threshold = analysis_params$logfc_threshold,
    output_file = markers_file
  )
  message(paste("标记基因已导出至:", markers_file))
  markers_file
}

prompt_for_annotation <- function(markers_file, output_dir) {
  template_file <- file.path(output_dir, "annotation_template.txt")
  prepare_annotation_template(markers_file, template_file)
  
  message("\n=== 请准备细胞类型注释文件 ===")
  message("1. 已生成标记基因文件:", markers_file)
  message("2. 模板文件已生成:", template_file)
  message("3. 编辑后另存为 annotation.txt 至输出目录")
  
  input <- readline("输入 'Y' 继续或 'Q' 退出: ")
  if (tolower(input) == "q") return(NULL)
  
  default_annotation_file <- file.path(output_dir, "annotation.txt")
  if (!file.exists(default_annotation_file)) {
    stop("未找到 annotation.txt")
  }
  default_annotation_file
}


### ----------------------
### 辅助函数：结果保存与可视化
### ----------------------
save_final_objects <- function(pbmc, output_dir) {
  saveRDS(pbmc, file.path(output_dir, "final_processed_data.rds"))
  message("最终数据已保存为 final_processed_data.rds")
}

generate_visualizations <- function(pbmc, output_dir, custom_genes, umap_dims) {
  # UMAP 可视化
  DimPlot(pbmc, reduction = "umap", group.by = "cell_type_OK") %>% 
    ggsave(file.path(output_dir, "umap_cell_types.png"), width = 10, height = 8)
  
  # 自定义基因 DotPlot
  if (!is.null(custom_genes)) {
    DotPlot(pbmc, features = custom_genes, group.by = "cell_type_OK") %>% 
      ggsave(file.path(output_dir, "dotplot_custom_genes.png"), width = 12, height = 10)
  }
}
