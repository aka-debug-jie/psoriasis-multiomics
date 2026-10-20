# 交互教学库按 R 版本与平台隔离，不更新用户其他项目。
course_root <- function(root=NULL) {
  candidates <- if(is.null(root)) c(getwd(),dirname(getwd())) else root
  roots <- candidates[file.exists(file.path(candidates,'teaching.Rproj'))]
  if(!length(roots))stop('请从教学项目或 notebooks 目录打开课程。')
  normalizePath(roots[1],winslash='/',mustWork=TRUE)
}
course_library <- function(root=course_root()) {
  minor <- strsplit(R.version$minor,'.',fixed=TRUE)[[1]][1]
  file.path(root,'.library',paste0('R-',R.version$major,'.',minor),R.version$platform)
}
course_dependencies <- function(lessons=sprintf('%02d',1:12),frontend='both') {
  groups <- list(
    core=c('data.table','ggplot2','scales','jsonlite','digest'),
    `01`=c('edgeR','limma'),`02`=character(),`03`=character(),`04`='plyr',
    `05`=character(),`06`='coloc',`07`=c('Seurat','Matrix','edgeR','limma'),
    `08`=character(),`09`='survival',`10`='igraph',`11`=character(),`12`=character())
  lessons <- sprintf('%02d',as.integer(lessons))
  front <- switch(frontend,jupyter='IRkernel',rstudio=c('knitr','rmarkdown'),both=c('IRkernel','knitr','rmarkdown'))
  unique(c(groups$core,unlist(groups[intersect(lessons,names(groups))]),front))
}
course_use_library <- function(root=course_root()) {
  if(Sys.getenv('TEACHING_BUNDLED')=='1') {
    lib <- Sys.getenv('TEACHING_LIBRARY')
    if(!dir.exists(lib))stop('内置包库缺失，请重新解压完整课程包。')
    .libPaths(c(lib,.Library));return(invisible(lib))
  }
  lib <- course_library(root)
  if(dir.exists(lib)).libPaths(c(lib,.libPaths()))
  invisible(lib)
}
course_package_status <- function(lessons=sprintf('%02d',1:12),frontend='both') {
  p <- course_dependencies(lessons,frontend)
  available <- vapply(p,requireNamespace,logical(1),quietly=TRUE)
  version <- vapply(p,function(k)if(available[[k]])as.character(utils::packageVersion(k)) else NA_character_,'')
  data.frame(package=p,available=unname(available),version=unname(version))
}
course_check_R <- function() {
  v <- getRversion()
  if(v<'4.4.0'||v>='4.7.0')stop('交互课程当前适配 R 4.4–4.6。请另外安装支持的 R，保留现有 R 即可。')
  invisible(TRUE)
}
