args<-commandArgs(trailingOnly=FALSE)
script<-sub('^--file=','',grep('^--file=',args,value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),'../..'))
source(file.path(root,'research/scripts/mr_core.R'))
input<-file.path(root,'research/data/processed/substitute_mr')
out<-file.path(root,'research/results/mr_fenland_substitute');dir.create(out,recursive=TRUE,showWarnings=FALSE)
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
write_tsv<-function(d,n)write.table(d,file.path(out,n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
selection<-jsonlite::fromJSON(file.path(input,'selection.json'))
x<-read.delim(file.path(input,'frozen_instruments.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
x$instrument_id<-paste(x$cohort,x$uniprot,sep='::')
stopifnot(!anyDuplicated(x$instrument_id))
registry<-read.delim(file.path(input,'protein_registry.tsv'),stringsAsFactors=FALSE)
threshold<-selection$bonferroni_threshold
rejections<-list();effects<-list();metas<-list()
for(outcome in c('mvp','finngen','stuart')){
  y<-read.delim(file.path(input,paste0(outcome,'_outcome.tsv')),stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(y$rsid))
  e<-x[x$eligibility=='eligible',]
  if(outcome!='mvp')e<-e[e$instrument_id %in% effects$mvp$instrument_id,]
  e<-do.call(rbind,lapply(seq_len(nrow(e)),function(i){ids<-strsplit(e$rsid[i],',',fixed=TRUE)[[1]];r<-e[rep(i,length(ids)),];r$rsid<-ids;r}))
  d<-merge(e,y,by='rsid',all.x=TRUE,suffixes=c('_exposure','_outcome'))
  orient<-lapply(seq_len(nrow(d)),function(i)orient_alleles(d$EA_exposure[i],d$OA_exposure[i],d$EA_outcome[i],d$OA_outcome[i],d$beta_outcome[i]))
  d$harmonization<-vapply(orient,`[[`,character(1),'status')
  d$beta_y_aligned<-vapply(orient,`[[`,numeric(1),'beta')
  d$status<-ifelse(is.na(d$EA_outcome),'frozen_index_missing_from_outcome',d$harmonization)
  valid<-is.finite(d$beta_y_aligned)&is.finite(d$se)&d$se>0
  # Multiple explicit aliases for the same source index are allowed only when
  # the outcome effect is identical. No secondary exposure index is substituted.
  groups<-split(which(valid),d$instrument_id[valid])
  conflicts<-unlist(lapply(groups,function(i)if(nrow(unique(d[i,c('beta_y_aligned','se')]))>1)i else integer()))
  if(length(conflicts)){valid[conflicts]<-FALSE;d$status[conflicts]<-'conflicting_outcome_aliases'}
  write_tsv(d,paste0(outcome,'_harmonization.tsv'))
  d<-d[valid,];d<-d[!duplicated(d$instrument_id),]
  w<-wald_ratio(d$beta_exposure,d$se_exposure,d$beta_y_aligned,d$se)
  names(w)[1:2]<-c('beta_mr','se_mr')
  d<-cbind(d,w);effects[[outcome]]<-d
  write_tsv(d,paste0(outcome,'_cohort_wald.tsv'))
  combine<-function(g)cbind(uniprot=g$uniprot[1],gene=g$gene[1],cohorts=paste(sort(g$cohort),collapse=';'),fixed_effect(g$beta_mr,g$se_mr))
  m<-do.call(rbind,lapply(split(d,d$uniprot),combine))
  m$passes_frozen_bonferroni<-m$p<threshold
  m$passes_paper_threshold<-m$p<1.89e-5
  m$scope<-'authorized_Fenland_substitute_estimated_scale_no_Steiger'
  metas[[outcome]]<-m
  write_tsv(m,paste0(outcome,'_meta.tsv'))
  two<-d[d$cohort!='Fenland_Olink_substitute_485',]
  write_tsv(do.call(rbind,lapply(split(two,two$uniprot),combine)),paste0(outcome,'_two_source_control.tsv'))
  f<-d[d$cohort=='Fenland_Olink_substitute_485',]
  raw<-wald_ratio(f$raw_beta,f$raw_se,f$beta_y_aligned,f$se)
  write_tsv(cbind(f[c('uniprot','gene','assay','rsid','raw_beta','raw_se','sdY_estimated')],raw),paste0(outcome,'_Fenland_raw_unit_wald.tsv'))
}
m<-metas$mvp
for(o in c('finngen','stuart')){
  r<-metas[[o]][c('uniprot','beta','p')];names(r)[2:3]<-paste0(c('beta_','p_'),o)
  m<-merge(m,r,by='uniprot',all.x=TRUE)
}
m$externally_supported<-(is.finite(m$p_finngen)&m$p_finngen<.05)|(is.finite(m$p_stuart)&m$p_stuart<.05)
m$replication_direction_concordant<-(is.finite(m$p_finngen)&m$p_finngen<.05&m$beta*m$beta_finngen>0)|(is.finite(m$p_stuart)&m$p_stuart<.05&m$beta*m$beta_stuart>0)
m<-merge(m,registry[c('uniprot','gene_id','gene_chr','gene_start','gene_end')],by='uniprot',all.x=TRUE)
write_tsv(m,'all_proteins_comparison.tsv')
chosen<-m[m$passes_frozen_bonferroni&m$externally_supported,]
write_tsv(chosen,'followup_proteins.tsv')
pdf(file.path(out,'substitute_MR.pdf'),width=8,height=6)
plot(m$beta,m$neg_log10_p,pch=16,cex=.5,col=ifelse(m$uniprot %in% chosen$uniprot,'#B2182B','#666666'),xlab='Combined log odds ratio',ylab='-log10(P)',main='Public pQTL reconstruction with Fenland Olink substitute')
abline(h=-log10(threshold),lty=2);dev.off()
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
inputs<-c(script,file.path(root,'research/scripts/mr_core.R'),file.path(input,c('selection.json','frozen_instruments.tsv','protein_registry.tsv','outcome_extraction.json',paste0(c('mvp','finngen','stuart'),'_outcome.tsv'))))
record<-list(status='substituted_three_source_MR_computed',original_data_reproduction=FALSE,
             frozen_denominator=selection$fixed_testing_denominator,bonferroni_threshold=threshold,
             proteins_by_outcome=lapply(metas,nrow),cohort_estimates=lapply(effects,nrow),
             discovery=sum(m$passes_frozen_bonferroni),externally_supported=nrow(chosen),
             steiger_applied=FALSE,missing_individual_validation='UKB PGS/Cox omitted by user',
             unit_boundary='Fenland per estimated SD, not measured SD; raw-unit single-source and two-source controls supplied; SD estimation uncertainty omitted',
             inputs=setNames(as.list(vapply(inputs,sha,character(1))),inputs),
             outputs=setNames(as.list(vapply(list.files(out,full.names=TRUE,pattern='[.](tsv|pdf|txt)$'),sha,character(1))),list.files(out,pattern='[.](tsv|pdf|txt)$')))
jsonlite::write_json(record,file.path(out,'run.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
print(record[c('status','proteins_by_outcome','discovery','externally_supported')])
