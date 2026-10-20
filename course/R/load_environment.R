# Notebook默认使用兼容教学库；对照入口显式选择reference。
profile <- Sys.getenv('TEACHING_PROFILE','interactive')
if(profile=='reference') {
  Sys.setenv(RENV_PATHS_CACHE=file.path(getwd(),'.renv-cache'),RENV_CONFIG_CACHE_SYMLINKS='FALSE')
  if(dir.exists('.bootstrap')) .libPaths(c(normalizePath('.bootstrap'),.libPaths()))
  if(dir.exists('renv/library') && requireNamespace('renv',quietly=TRUE)) renv::load(project=getwd())
} else {
  source('R/environment.R')
  course_use_library()
}
