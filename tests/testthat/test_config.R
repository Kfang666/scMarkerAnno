context("Test parameter configuration")

test_that("set_analysis_params returns valid parameters", {
  params <- set_analysis_params()
  expect_is(params, "list")
  expect_named(params, c("cell_threshold", "pca_dims", "resolution_setting", 
                         "parent_marker", "min_pct", "logfc_threshold"))
})

test_that("invalid parameters are caught", {
  expect_error(set_analysis_params(cell_threshold = -100))
  expect_error(set_analysis_params(min_pct = 1.5))
})

# tests/testthat/test_grouping.R
context("Test sample grouping functions")

test_that("read_sample_groups works correctly", {
  # 创建临时文件
  temp_file <- tempfile(fileext = ".txt")
  writeLines("sample\tgroup\nSample1\tGroup1\nSample2\tGroup1\nSample3\tGroup2", temp_file)
  
  groups <- read_sample_groups(temp_file)
  
  expect_is(groups, "list")
  expect_named(groups, c("Group1", "Group2"))
  expect_equal(groups$Group1, c("Sample1", "Sample2"))
  expect_equal(groups$Group2, "Sample3")
})

test_that("read_sample_groups handles missing file", {
  expect_error(read_sample_groups("nonexistent_file.txt"))
})






library(scMarkerAnno)
# 设置分析参数
params <- set_analysis_params(
  cell_threshold = 1500,
  resolution_setting = 0.8,
  min_pct = 0.3
)
# 运行主分析流程
results <- main_workflow(
  data_dir = "~/workspace/NSCLC/OK/In-house-MW/paired-sample",
  sample_groups_file = "~/sample_groups.txt",
  analysis_params = params,
  output_dir = "~/scMarkerAnno_results"
)

# 查看结果
print(results)
DimPlot(results, group_by = "cell_type")