course_platform <- function(system=Sys.info()[['sysname']],arch=R.version$arch) {
  os<-switch(tolower(system),windows='windows',darwin='macos',linux='linux',tolower(system))
  cpu<-if(tolower(arch)%in%c('x86_64','amd64','x64'))'x86_64' else if(tolower(arch)%in%c('aarch64','arm64'))'arm64' else tolower(arch)
  list(os=os,arch=cpu)
}
course_binary_platform <- function(path) {
  con<-file(path,'rb');on.exit(close(con));b<-readBin(con,'raw',n=64)
  hex<-paste(sprintf('%02x',as.integer(b[1:4])),collapse='')
  integer_at<-function(at,size=4,endian='little') {
    seek(con,where=at,origin='start');readBin(con,integer(),n=1,size=size,signed=(size!=2),endian=endian)
  }
  arch_name<-function(n)if(n%in%c(62L,34404L,16777223L))'x86_64' else if(n%in%c(183L,43620L,16777228L))'arm64' else 'unknown'
  if(hex=='7f454c46')return(list(os='linux',arch=arch_name(integer_at(18,2))))
  if(hex=='cffaedfe')return(list(os='macos',arch=arch_name(integer_at(4))))
  if(hex=='feedfacf')return(list(os='macos',arch=arch_name(integer_at(4,endian='big'))))
  if(hex%in%c('cafebabe','cafebabf')) {
    n<-integer_at(4,endian='big');stride<-if(hex=='cafebabf')32L else 20L
    types<-vapply(seq_len(n),function(i)arch_name(integer_at(8+(i-1)*stride,endian='big')),'')
    return(list(os='macos',arch=types))
  }
  if(paste(sprintf('%02x',as.integer(b[1:2])),collapse='')=='4d5a') {
    offset<-integer_at(60);return(list(os='windows',arch=arch_name(integer_at(offset+4,2))))
  }
  list(os='unknown',arch='unknown')
}
tool_path <- function(name) {
  platform<-course_platform();override<-Sys.getenv(paste0('TEACHING_',toupper(name)),'')
  manifest<-read.delim('tools/manifest.tsv',stringsAsFactors=FALSE)
  row<-manifest[manifest$tool==name & manifest$platform==platform$os & manifest$arch==platform$arch,]
  path<-if(nzchar(override))override else if(nrow(row)==1)row$path else ''
  if(!nzchar(path)||!file.exists(path))stop('未找到 ',name,' 的 ',platform$os,'/',platform$arch,
    ' 程序；请重新解压完整课程包。',call.=FALSE)
  p<-normalizePath(path,winslash='/',mustWork=TRUE);actual<-course_binary_platform(p)
  if(actual$os!=platform$os || !platform$arch%in%actual$arch)stop('程序与当前系统/架构不匹配：',name,'。',call.=FALSE)
  if(!nzchar(override)&&nrow(row)==1&&!is.na(row$sha256)&&nzchar(row$sha256)) {
    if(digest::digest(file=p,algo='sha256')!=row$sha256)stop('程序校验值不符：',name,call.=FALSE)
  }
  p
}
