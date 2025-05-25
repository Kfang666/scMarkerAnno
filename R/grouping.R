#' Read sample groups from a text file
#'
#' @param file_path Path to the group definition file
#' @param sep Separator used in the file (default: tab)
#' @param header Does the file contain a header line?
#' @return A named list of sample groups
#' @export
#' @examples
#' groups <- read_sample_groups("path/to/group.txt")
read_sample_groups <- function(file_path, sep = "\t", header = TRUE) {
  # 验证文件存在
  if (!file.exists(file_path)) {
    stop(paste("Group file not found:", file_path))
  }
  # 读取分组信息
  group_data <- read.table(
    file = file_path,
    sep = sep,
    header = header,
    stringsAsFactors = FALSE
  )
  # 验证必需的列存在
  if (!all(c("sample", "group") %in% colnames(group_data))) {
    stop("Group file must contain columns 'sample' and 'group'")
  }
  
  # 转换为命名列表
  group_names <- unique(group_data$group)
  sample_groups <- lapply(group_names, function(grp) {
    group_data$sample[group_data$group == grp]
  })
  
  names(sample_groups) <- group_names
  
  return(sample_groups)
}
