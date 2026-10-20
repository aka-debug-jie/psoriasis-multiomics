test_arg<-grep('^--file=',commandArgs(),value=TRUE)
test_file<-normalizePath(sub('^--file=','',test_arg[1]),mustWork=TRUE)
support_file<-file.path(dirname(dirname(test_file)),'R','notebook_support.R')
source(support_file,local=TRUE)

root<-tempfile('student-progress-')
dir.create(root)
dir.create(file.path(root,'notebook_templates'))
dir.create(file.path(root,'lessons'))
writeLines('print("starter")',file.path(root,'notebook_templates','01_intro.ipynb'))
writeLines('lesson',file.path(root,'lessons','01_intro.R'))
options(teaching.root=root,teaching.lesson='01')

first<-course_fresh_notebook('01')
first_path<-file.path(root,first)
stopifnot(file.exists(first_path),startsWith(first,'notebooks/01_intro_练习_'))
writeLines('student edits',first_path)
second<-course_fresh_notebook('01')
stopifnot(first!=second,file.exists(file.path(root,second)),
          identical(readLines(first_path), 'student edits'))

output_root<-file.path(root,'results')
source_dir<-file.path(output_root,'01_intro')
dir.create(file.path(source_dir,'nested'),recursive=TRUE)
writeLines('result',file.path(source_dir,'result.txt'))
writeLines('nested result',file.path(source_dir,'nested','detail.txt'))
Sys.setenv(TEACHING_OUTPUT=output_root)
settings<-list(alpha=0.25,cohort='测试')
snap1<-course_snapshot('复核_run-1',settings)
snap2<-course_snapshot('复核_run-1',settings)
stopifnot(snap1!=snap2)
snap_path<-snap1
stopifnot(identical(readLines(file.path(snap_path,'result.txt')),'result'),
          identical(readLines(file.path(snap_path,'nested','detail.txt')),'nested result'),
          identical(readRDS(file.path(snap_path,'parameters.rds')),settings),
          file.exists(file.path(snap_path,'parameters.txt')),
          file.exists(file.path(snap_path,'snapshot_info.tsv')))
before<-tools::md5sum(c(file.path(source_dir,'result.txt'),file.path(source_dir,'nested','detail.txt')))
for(bad in c('../bad','a/b','a\\b','bad..label')) {
  rejected<-inherits(try(course_snapshot(bad),silent=TRUE),'try-error')
  stopifnot(rejected)
}
after<-tools::md5sum(c(file.path(source_dir,'result.txt'),file.path(source_dir,'nested','detail.txt')))
stopifnot(identical(before,after))

unlink(source_dir,recursive=TRUE)
snapshot_root<-file.path(dirname(output_root),'snapshots','01_intro')
existing_snapshots<-if(dir.exists(snapshot_root))list.files(snapshot_root,full.names=TRUE)else character()
missing_error<-inherits(try(course_snapshot('missing'),silent=TRUE),'try-error')
current_snapshots<-if(dir.exists(snapshot_root))list.files(snapshot_root,full.names=TRUE)else character()
stopifnot(missing_error,identical(existing_snapshots,current_snapshots))
unlink(root,recursive=TRUE)
cat('test_student_progress.R: PASS\n')
