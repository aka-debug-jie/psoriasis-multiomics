source('R/common.R')
source('R/genetics.R')
out<-lesson_output('06_colocalization')
rows<-list();points<-list();k<-0L
for(gene in c('IL7R','ERAP2')) for(window in c(500,1000)) {
  d<-read_tsv(paste0('data/coloc/',gene,'_',window,'kb.tsv.gz'))
  stopifnot(!anyDuplicated(d$SNP),nrow(d)>1000L)
  d$by<-d$by*mapply(allele_sign,d$A1,d$A2,d$EA_y,d$OA_y)
  stopifnot(all(is.finite(d$by)))
  for(p12 in c(1e-6,1e-5,1e-4)) {
    fit<-coloc::coloc.abf(
      list(snp=d$SNP,position=d$position,type='quant',beta=d$bx,varbeta=d$sx^2,sdY=1),
      list(snp=d$SNP,position=d$position,type='cc',beta=d$by,varbeta=d$sy^2,s=15635/444887),
      p1=1e-4,p2=1e-4,p12=p12)
    pp<-fit$summary[paste0('PP.H',0:4,'.abf')]
    stopifnot(abs(sum(pp)-1)<1e-10,all(pp>=0 & pp<=1))
    k<-k+1;rows[[k]]<-data.frame(gene=gene,window_kb=window,n_variants=nrow(d),p12=p12,
      as.list(pp),strong_H4=unname(pp['PP.H4.abf'])>.8,check.names=FALSE)
    if(p12==1e-5)write_tsv(fit$results,file.path(out,paste0(gene,'_',window,'kb_posteriors.tsv')))
  }
  if(window==500)points[[gene]]<-rbind(data.frame(gene,position=d$position,signal='Molecular trait',neglogp=-log10(pmax(2*pnorm(-abs(d$bx/d$sx)),1e-300))),
    data.frame(gene,position=d$position,signal='Psoriasis',neglogp=-log10(pmax(2*pnorm(-abs(d$by/d$sy)),1e-300))))
}
res<-do.call(rbind,rows);write_tsv(res,file.path(out,'coloc_summary.tsv'))
p<-ggplot2::ggplot(do.call(rbind,points),ggplot2::aes(position/1e6,neglogp,colour=signal))+
 ggplot2::geom_point(size=.6,alpha=.6)+ggplot2::facet_grid(signal~gene,scales='free')+
 ggplot2::scale_colour_manual(values=c('Molecular trait'='#148F8A',Psoriasis='#394D72'))+
 ggplot2::labs(x='Position (GRCh37, Mb)',y=expression(-log[10](P)),title='Regional association patterns',subtitle='All shared, harmonizable variants within index ±500 kb')+theme_course()+ggplot2::theme(legend.position='none')
save_plot(p,out,'regional_associations',9,6)
print(res)
