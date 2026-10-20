source('R/load_environment.R')
# 在 teaching.Rproj 打开的项目根目录中运行各课脚本。
options(stringsAsFactors = FALSE, scipen = 6)
Sys.setenv(OMP_NUM_THREADS = "2", OPENBLAS_NUM_THREADS = "2", MKL_NUM_THREADS = "2")
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(2)
read_tsv <- function(path) read.delim(if (grepl("\\.gz$", path)) gzfile(path) else path,
  check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
  quote = FALSE, na = "NA", fileEncoding = "UTF-8")
lesson_output <- function(id) {
  path <- file.path(Sys.getenv("TEACHING_OUTPUT", "results"), id)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
theme_course <- function() ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(colour = "#566573"),
    axis.line = ggplot2::element_line(linewidth = .35),
    legend.position = "bottom", plot.margin = ggplot2::margin(12, 16, 12, 12))
save_plot <- function(p, output, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(output, paste0(name, ".pdf")), p, width = width, height = height)
  ggplot2::ggsave(file.path(output, paste0(name, ".png")), p, width = width, height = height,
                 dpi = 320, bg = "white")
  if(isTRUE(getOption("teaching.inline",FALSE))) {
    options(repr.plot.width=width,repr.plot.height=height)
    print(p)
  }
}
source('R/tool_support.R')
run_tool <- function(program, args, log) {
  code <- system2(program, vapply(args, shQuote, ""), stdout = log, stderr = log)
  if (code != 0) stop("外部程序运行失败；请查看 ", log)
  invisible(code)
}
