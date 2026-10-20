source('R/common.R')
source('R/genetics.R')
out<-lesson_output('03_proteome_mr')
x<-read_tsv('data/genetics/exposures.tsv')
ld<-load_ld()
# 2709 是事先固定的受检蛋白家族；不按显著性或结局可得性缩小分母。
stopifnot(length(unique(x$uniprot))==2709L)
plan<-freeze_instruments(x,ld,max_abs_r=.95)
write_tsv(plan,file.path(out,'selected_instruments.tsv'))
results<-list();harmonized<-list()
for(label in c('mvp','stuart')) {
  y<-read_tsv(paste0('data/genetics/',label,'.tsv'))
  fit<-analyse_outcome(plan,y,ld)
  fit$results$outcome<-label;fit$harmonized$outcome<-label
  results[[label]]<-fit$results;harmonized[[label]]<-fit$harmonized
}
all<-do.call(rbind,results)
write_tsv(all,file.path(out,'all_outcomes.tsv'))
write_tsv(do.call(rbind,harmonized),file.path(out,'harmonization.tsv'))
main<-results$mvp
for(label in c('stuart')) {
  v<-results[[label]][match(main$uniprot,results[[label]]$uniprot),]
  main[[paste0(label,'_beta')]]<-v$beta;main[[paste0(label,'_p')]]<-v$p
}
main$discovery<-main$p<.05/2709
main$external_nominal<-main$stuart_p<.05
main$candidate<-main$discovery & main$external_nominal
main$direction_agrees<-  ((main$stuart_p<.05)&main$beta*main$stuart_beta>0)
write_tsv(main,file.path(out,'all_proteins.tsv'))
selected<-main[which(main$candidate),]
write_tsv(selected,file.path(out,'candidates.tsv'))
plotdata<-all[all$uniprot %in% selected$uniprot,]
plotdata$gene<-factor(plotdata$gene,levels=rev(selected$gene[order(selected$beta)]))
p<-ggplot2::ggplot(plotdata,ggplot2::aes(x=beta,y=gene,colour=outcome))+
  ggplot2::geom_vline(xintercept=0,colour='#A8B2B8',linetype=2)+
  ggplot2::geom_pointrange(ggplot2::aes(xmin=CI_log_low,xmax=CI_log_high),position=ggplot2::position_dodge(.55))+
  ggplot2::scale_colour_manual(values=c(mvp='#184E77',stuart='#D97948'))+
  ggplot2::labs(x='Log odds ratio per unit protein',y=NULL,colour='Outcome',title='Discovery and external outcomes',
    subtitle='One observation per variant; signed-LD generalized IVW')+theme_course()
save_plot(p,out,'discovery_and_validation',8,5.5)
print(selected[,c('gene','beta','se','p','stuart_p','direction_agrees')])
