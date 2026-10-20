# Jupyter 和 RStudio 共用初始化与显示；分析步骤保留在每课单元格内。
course_init <- function(lesson='00',root=NULL,frontend=Sys.getenv('TEACHING_FRONTEND','rstudio')) {
  if(is.null(root)) {
    candidates<-c(getwd(),dirname(getwd()))
    root<-candidates[file.exists(file.path(candidates,'teaching.Rproj'))][1]
  }
  if(is.na(root))stop('请从教学项目或 notebooks 目录打开课程。')
  if(Sys.getenv('TEACHING_PROFILE','interactive')=='reference' && any(c('Seurat','edgeR','limma','Matrix') %in% loadedNamespaces()))
    stop('从对照环境切换到交互教学前，请重启R会话。')
  root<-normalizePath(root,winslash='/',mustWork=TRUE)
  source(file.path(root,'R/environment.R'),local=.GlobalEnv)
  course_check_R()
  Sys.setenv(TEACHING_PROFILE='interactive')
  course_use_library(root)
  # knitr 每个 chunk 恢复工作目录，必须在 setup 中设置 root.dir。
  if(requireNamespace('knitr',quietly=TRUE)) {
    knitr::opts_knit$set(root.dir=root)
    knitr::opts_chunk$set(cache=FALSE,error=FALSE,warning=TRUE,message=FALSE,fig.width=10,fig.height=6)
  }
  setwd(root)
  options(teaching.root=root,teaching.lesson=lesson,teaching.inline=TRUE,teaching.reference=FALSE,
          repr.plot.width=10,repr.plot.height=6,repr.plot.res=120,jupyter.plot_mimetypes="image/png",repr.matrix.max.rows=10,repr.matrix.max.cols=12)
  outroot<-Sys.getenv('TEACHING_NOTEBOOK_OUTPUT',file.path(root,'results','interactive',frontend))
  Sys.setenv(TEACHING_OUTPUT=outroot,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
  if(lesson!='00') {
    packages<-course_package_status(lesson,frontend)
    if(any(!packages$available) && Sys.getenv('TEACHING_BUNDLED')=='1')stop('课程依赖不完整，请重新解压完整课程包。')
    if(any(!packages$available))stop('本课缺少 ',paste(packages$package[!packages$available],collapse=', '),
      '。在起步页运行 setup_interactive(lessons="',lesson,'", frontend="',frontend,'") 后重启会话。')
  }
  source(file.path(root,'R/common.R'),local=.GlobalEnv)
  if(lesson=='05')course_require_tool('smr')
  if(lesson=='08')course_require_tool('plink')
  invisible(data.frame(lesson=lesson,R=as.character(getRversion()),frontend=frontend,
      profile='interactive',threads=2L))
}
course_table <- function(x,n=8,columns=NULL) {
  if(!is.null(columns))x<-x[,columns,drop=FALSE]
  utils::head(as.data.frame(x),n)
}
course_require_tool <- function(name) {
  tryCatch(tool_path(name),error=function(e)stop('本课需要匹配系统与架构的 ',toupper(name),
    '。请查看起步页“外部工具”；其他课程可以继续。\n',conditionMessage(e),call.=FALSE))
}
course_status <- function() {
  p<-course_package_status(frontend=Sys.getenv('TEACHING_FRONTEND','rstudio'))
  p
}

course_fresh_notebook <- function(lesson) {
  root<-getOption('teaching.root')
  if(is.null(root) || length(root)!=1L || !nzchar(root))stop('缺少 teaching.root 课程根目录。',call.=FALSE)
  id<-as.character(lesson)
  if(length(id)!=1L || is.na(id) || !grepl('^[0-9]{2}$',id))stop('lesson 必须是两位课程编号。',call.=FALSE)
  templates<-list.files(file.path(root,'notebook_templates'),pattern=paste0('^',id,'_.*\\.ipynb$'),full.names=TRUE)
  if(length(templates)!=1L)stop('课程 ',id,' 必须且只能匹配一个 notebook_templates 模板；找到 ',length(templates),' 个。',call.=FALSE)
  stem<-tools::file_path_sans_ext(basename(templates))
  destdir<-file.path(root,'notebooks')
  dir.create(destdir,recursive=TRUE,showWarnings=FALSE)
  stamp<-format(Sys.time(),'%Y%m%d%H%M%S')
  i<-1L
  repeat {
    dest<-file.path(destdir,sprintf('%s_练习_%s_%03d.ipynb',stem,stamp,i))
    if(!file.exists(dest))break
    i<-i+1L
  }
  if(!file.copy(templates,dest,overwrite=FALSE))stop('复制 notebook 模板失败。',call.=FALSE)
  relative<-substring(normalizePath(dest,winslash='/',mustWork=TRUE),nchar(normalizePath(root,winslash='/',mustWork=TRUE))+2L)
  message('已创建练习 notebook：',relative)
  invisible(relative)
}

course_snapshot <- function(label,settings=list()) {
  if(length(label)!=1L || is.na(label) || nchar(label,type='chars')>24L || !grepl('^[\\p{Han}A-Za-z0-9_-]+$',label,perl=TRUE) || grepl('[/\\\\]|\\.\\.',label))
    stop('label 最长 24 个字符，仅可包含中英文、数字、下划线和横线，且不能包含路径分隔符或 ..。',call.=FALSE)
  root<-getOption('teaching.root')
  id<-as.character(getOption('teaching.lesson'))
  if(is.null(root) || length(id)!=1L || is.na(id) || !grepl('^[0-9]{2}$',id))stop('缺少有效的 teaching.root 或 teaching.lesson。',call.=FALSE)
  scripts<-list.files(file.path(root,'lessons'),pattern=paste0('^',id,'_.*\\.R$'))
  if(length(scripts)!=1L)stop('课程 ',id,' 必须且只能匹配一个 lessons R 脚本；找到 ',length(scripts),' 个。',call.=FALSE)
  stem<-sub('\\.R$','',scripts)
  output_root<-Sys.getenv('TEACHING_OUTPUT')
  if(!nzchar(output_root))stop('TEACHING_OUTPUT 未设置。',call.=FALSE)
  source_dir<-file.path(output_root,stem)
  if(!dir.exists(source_dir))stop('当前课程没有结果目录：',source_dir,call.=FALSE)
  entries<-list.files(source_dir,all.files=TRUE,no..=TRUE,full.names=TRUE)
  if(!length(entries))stop('当前课程结果目录为空，未创建快照：',source_dir,call.=FALSE)
  snapshot_parent<-file.path(dirname(output_root),'snapshots',stem)
  dir.create(snapshot_parent,recursive=TRUE,showWarnings=FALSE)
  stamp<-format(Sys.time(),'%Y%m%d%H%M%S')
  i<-1L
  repeat {
    dest<-file.path(snapshot_parent,sprintf('%s_%s_%03d',stamp,label,i))
    if(!file.exists(dest))break
    i<-i+1L
  }
  if(!dir.create(dest,recursive=FALSE,showWarnings=FALSE))stop('创建快照目录失败：',dest,call.=FALSE)
  copied<-file.copy(entries,dest,recursive=TRUE,copy.mode=TRUE,copy.date=TRUE)
  if(!all(copied))stop('复制课程结果到快照时失败。',call.=FALSE)
  saveRDS(settings,file.path(dest,'parameters.rds'))
  dput(settings,file=file.path(dest,'parameters.txt'))
  created_at<-format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')
  utils::write.table(data.frame(lesson=id,label=label,created_at=created_at),
    file=file.path(dest,'snapshot_info.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  relative<-substring(normalizePath(dest,winslash='/',mustWork=TRUE),nchar(normalizePath(dirname(output_root),winslash='/',mustWork=TRUE))+2L)
  message('已保存课程快照：',relative)
  invisible(dest)
}

course_output_path <- function() {
  root<-getOption('teaching.root');id<-getOption('teaching.lesson')
  f<-list.files(file.path(root,'lessons'),pattern=paste0('^',id,'_.*\\.R$'))
  slug<-if(length(f))sub('\\.R$','',f[1])else ''
  p<-file.path(Sys.getenv('TEACHING_OUTPUT'),slug)
  if(Sys.getenv('TEACHING_BUNDLED')=='1')return(paste('当前输出目录：课程资源/work/results/',slug,sep=''))
  if(startsWith(p,paste0(root,'/')))p<-substring(p,nchar(root)+2)
  paste('当前输出目录：',p)
}

course_self_check <- function() {
  packages <- course_package_status(frontend='both')
  output <- Sys.getenv('TEACHING_OUTPUT')
  dir.create(output, recursive=TRUE, showWarnings=FALSE)
  probe <- tempfile('write-check-', tmpdir=output)
  writable <- file.create(probe)
  if (writable) unlink(probe)
  tools <- vapply(c('smr','plink'), function(name) {
    tryCatch({tool_path(name); TRUE}, error=function(e) FALSE)
  }, logical(1))
  data.frame(check=c('R version','Required R packages','SMR','PLINK','Save results'),
    value=c(as.character(getRversion()), paste(sum(packages$available),'/',nrow(packages)),
      ifelse(tools,'available','missing'),ifelse(writable,'writable','not writable')),
    passed=c(getRversion()=='4.4.2',all(packages$available),tools,writable), row.names=NULL)
}
