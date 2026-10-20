# Build-time only. Nothing in this script is executed by a student launcher.
args<-commandArgs(trailingOnly=TRUE)
runtime<-normalizePath(args[1],winslash='/',mustWork=TRUE)
lib<-file.path(runtime,'course-library');dir.create(lib,showWarnings=FALSE)
.libPaths(c(lib,.Library))
Sys.setenv(RENV_CONFIG_CACHE_ENABLED='FALSE',RENV_CONFIG_AUTOLOADER_ENABLED='FALSE',
           RENV_CONFIG_SANDBOX_ENABLED='FALSE',RENV_CONFIG_PAK_ENABLED='FALSE')
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=600,Ncpus=2)
# R 4.4.2's macOS Makeconf may leave CC17 empty although Clang supports C17.
# Select that standard explicitly for packages that request it; package sources stay unchanged.
rcommand<-file.path(R.home('bin'),if(.Platform$OS.type=='windows')'R.exe'else'R')
cc<-system2(rcommand,c('CMD','config','CC'),stdout=TRUE)
cflags<-system2(rcommand,c('CMD','config','CFLAGS'),stdout=TRUE)
makevars<-file.path(runtime,'Makevars.build')
settings<-c(paste('CC17 =',paste(cc,collapse=' '),'-std=gnu17'),
            paste('C17FLAGS =',paste(cflags,collapse=' ')))
if(Sys.info()[['sysname']]=='Darwin') {
  # cctools' assembler wrapper misdetects gfortran's stdin input with -pipe.
  # Use a temporary assembly file; optimization and numerical flags are unchanged.
  for(name in c('FFLAGS','FCFLAGS')) {
    value<-paste(system2(rcommand,c('CMD','config',name),stdout=TRUE),collapse=' ')
    value<-gsub('(^|[[:space:]])-pipe(?=[[:space:]]|$)','\\1',value,perl=TRUE)
    settings<-c(settings,paste(name,'=',value))
  }
}
writeLines(settings,makevars)
Sys.setenv(R_MAKEVARS_USER=makevars)
probe_dir<-file.path(dirname(runtime),'compiler-probe');dir.create(probe_dir,showWarnings=FALSE)
probe<-file.path(probe_dir,'course_probe.f')
writeLines(c('      subroutine course_probe(x)','      double precision x','      x=2d0*x','      end'),probe)
object<-file.path(probe_dir,paste0('course_probe',.Platform$dynlib.ext))
status<-system2(rcommand,c('CMD','SHLIB',shQuote(probe),'-o',shQuote(object)))
stopifnot(status==0)
dyn.load(object);stopifnot(.Fortran('course_probe',x=as.double(2))$x==4);dyn.unload(object)


if(!requireNamespace('renv',quietly=TRUE))install.packages('renv',lib=lib)
renv::restore(lockfile='course/renv.lock',library=lib,prompt=FALSE,transactional=FALSE)
renv::install('IRkernel@1.3.2',library=lib,prompt=FALSE)
source('course/R/environment.R')
status<-course_package_status(frontend='both');print(status)
stopifnot(all(status$available),getRversion()=='4.4.2')
write.table(installed.packages(lib.loc=lib)[,c('Package','Version'),drop=FALSE],file.path(runtime,'packages.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
