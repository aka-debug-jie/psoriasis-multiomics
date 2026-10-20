script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'../..'))
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
out<-file.path(root,'research/results/covariance_resolution');fig<-file.path(out,'figures');dir.create(fig,showWarnings=FALSE)
read<-function(p)read.delim(file.path(root,p),stringsAsFactors=FALSE,check.names=FALSE)
sources<-c('research/results/covariance_resolution/all_proteins_comparison.tsv','research/results/covariance_resolution/clinical_original_17.tsv',
           'research/results/clinical_mr_substitute/all_meta_results.tsv','research/results/mvp_precision_audit/mvp_precision_audit.tsv')
d<-read(sources[1]);cl<-read(sources[2]);oldcl<-read(sources[3]);precision<-read(sources[4])
logp<-function(p)-log10(pmax(p,1e-300));cut<--log10(.05/2709)
theme_set(theme_minimal(base_family='Liberation Sans',base_size=8)+theme(panel.grid.minor=element_blank(),
 panel.grid.major=element_line(color='#E8EEF1',linewidth=.25),plot.title=element_text(face='bold',size=10,color='#243746'),
 plot.subtitle=element_text(size=7,color='#52616D'),axis.title=element_text(size=8),axis.text=element_text(color='#465564',size=7),
 legend.position='bottom',legend.title=element_text(size=7),legend.text=element_text(size=7),plot.margin=margin(8,8,8,8)))
base<-d[d$branch=='unique_variant_givw'&d$scenario=='CI_Wald',]
a<-base[is.finite(base$p)&is.finite(base$original_p),];a$x<-pmin(30,logp(a$original_p));a$y<-pmin(30,logp(a$p))
a$clipped<-logp(a$original_p)>30|logp(a$p)>30;a$group<-ifelse(a$gene %in% base$gene[base$followup=='True'],'Follow-up','Other')
p1<-ggplot(a,aes(x,y))+geom_abline(slope=1,intercept=0,color='#B8C3CB',linewidth=.35)+
 geom_hline(yintercept=cut,linetype=2,color='#AEB8C0',linewidth=.3)+geom_vline(xintercept=cut,linetype=2,color='#AEB8C0',linewidth=.3)+
 geom_point(aes(color=group,shape=clipped),size=1.1,alpha=.65)+scale_color_manual(values=c('Follow-up'='#D07945','Other'='#7897A3'),name=NULL)+
 scale_shape_manual(values=c('FALSE'=16,'TRUE'=17),guide='none')+coord_equal(xlim=c(0,30),ylim=c(0,30))+
 labs(x=expression('Original pooled '*-log[10]*italic(P)),y=expression('Unique-variant '*-log[10]*italic(P)),
 title='Shared outcomes counted once',subtitle=paste(nrow(a),'matched proteins; axes capped at 30'))
targets<-base[base$original_followup=='True',];genes<-targets$gene[order(targets$original_p)]
b<-rbind(data.frame(gene=targets$gene,method='Original',p=targets$original_p,followup=TRUE),
 data.frame(gene=targets$gene,method='LD-aware',p=targets$p,followup=targets$followup=='True'))
for(s in c('P_equivalent','CI_compatible_only')){
 z<-d[d$branch=='unique_variant_givw'&d$scenario==s&d$original_followup=='True',]
 b<-rbind(b,data.frame(gene=z$gene,method=if(s=='P_equivalent')'P-equivalent' else 'Compatible only',p=z$p,followup=z$followup=='True'))}
b$method<-factor(b$method,levels=c('Original','LD-aware','P-equivalent','Compatible only'));b$display<-pmin(15,logp(b$p))
p2<-ggplot(b,aes(method,gene,fill=display))+geom_tile(width=.93,height=.9,color='white',linewidth=.2)+
 geom_text(data=b[!is.na(b$followup)&b$followup,],label='+',size=2.6,color='white',fontface='bold')+
 geom_text(data=b[is.na(b$p),],label='\u00d7',size=2.4,color='#85929C')+
 scale_fill_gradient(low='#E6EEF1',high='#245A70',limits=c(0,15),na.value='#F1F3F4',name='-log10 P (cap 15)')+
 scale_y_discrete(limits=rev(genes))+labs(x=NULL,y=NULL,title='The original 17-protein set',subtitle='+ discovery and external nominal support; \u00d7 unavailable')+
 theme(axis.text.x=element_text(angle=35,hjust=1),panel.grid=element_blank())
c<-precision[tolower(precision$currently_used_rsid)=='true',];c$source_P<-as.numeric(c$original_p_value);c$ci_P<-as.numeric(c$p_ci_point_estimate)
c$x<-pmin(20,logp(c$source_P));c$y<-pmin(20,logp(c$ci_P));shown<-c[is.finite(c$x)&is.finite(c$y),]
p3<-ggplot(shown,aes(x,y))+geom_abline(slope=1,intercept=0,color='#AAB6BF',linewidth=.35)+
 geom_point(aes(color=status),size=1,alpha=.65)+scale_color_manual(values=c(consistent='#94A9B3',inconsistent='#D07945',unresolved='#7B6D96'),name=NULL)+
 coord_equal(xlim=c(0,20),ylim=c(0,20))+labs(x=expression('Source '*-log[10]*italic(P)),y=expression('CI-derived '*-log[10]*italic(P)),
 title='MVP uncertainty remains',subtitle='218 incompatible, 19 unresolved; axes capped at 20')
e<-data.frame(method=c('Original','LD-aware','LD-pruned'),BH05=c(sum(oldcl$BH_all_requested_pairs<.05),
 sum(cl$BH_fixed_374[cl$branch=='unique_variant_givw']<.05,na.rm=TRUE),sum(cl$BH_fixed_374[cl$branch=='ld_pruned_givw']<.05,na.rm=TRUE)),
 estimable=c(sum(is.finite(oldcl$p)),sum(is.finite(cl$p[cl$branch=='unique_variant_givw'])),sum(is.finite(cl$p[cl$branch=='ld_pruned_givw']))))
e$method<-factor(e$method,levels=e$method)
p4<-ggplot(e,aes(method,BH05,fill=method))+geom_col(width=.55)+geom_text(aes(label=paste0(BH05,' / ',estimable)),vjust=-.5,size=2.8,color='#243746')+
 scale_fill_manual(values=c('#9EAFB9','#2E7486','#7B6D96'),guide='none')+scale_y_continuous(expand=expansion(mult=c(0,.15)))+
 labs(x=NULL,y='Associations with BH FDR < 0.05',title='Clinical associations: original 17 proteins',subtitle='374 tests; labels show significant / estimable pairs')
plot<-(p1+p2)/(p3+p4)+plot_layout(heights=c(1.35,1))+plot_annotation(title='Resolving shared-outcome weighting',
 subtitle='Method extension; original outputs retained. P-equivalent and compatible-only panels are diagnostics.',tag_levels='a',
 theme=theme(plot.title=element_text(family='Liberation Sans',face='bold',size=15,color='#243746'),plot.subtitle=element_text(family='Liberation Sans',size=8),plot.tag=element_text(size=10,face='bold')))
for(ext in c('pdf','svg','png')){
 path<-file.path(fig,paste0('METHOD_COMPARISON.',ext))
 if(ext=='pdf')cairo_pdf(path,width=7.2,height=8.7,family='Liberation Sans')
 if(ext=='svg')svg(path,width=7.2,height=8.7,family='Liberation Sans')
 if(ext=='png')png(path,width=7.2,height=8.7,units='in',res=600,type='cairo',bg='white')
 print(plot);dev.off()
}
tables<-list(discovery=a,followup=b,mvp_ci=c,clinical=e)
for(n in names(tables))write.table(tables[[n]],file.path(fig,paste0(n,'.tsv')),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
files<-list.files(fig,full.names=TRUE);files<-files[basename(files)!='manifest.json'];record<-list(status='generated',width_in=7.2,height_in=8.7,png_dpi=600,
 inputs=setNames(as.list(vapply(file.path(root,c(sources,'research/scripts/plot_method_resolution.R')),sha,character(1))),c(sources,'research/scripts/plot_method_resolution.R')),
 outputs=setNames(as.list(vapply(files,sha,character(1))),basename(files)))
jsonlite::write_json(record,file.path(fig,'manifest.json'),pretty=TRUE,auto_unbox=TRUE)
