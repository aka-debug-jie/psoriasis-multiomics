script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'../..'))
source(file.path(root,'research/scripts/mr_core.R'))
substituted<-'--substitute' %in% commandArgs(trailingOnly=TRUE)
input<-file.path(root,if(substituted)'research/data/processed/coloc_protein_substitute' else 'research/data/processed/coloc_protein_fenland');out<-file.path(root,if(substituted)'research/results/coloc_protein_substitute' else 'research/results/coloc_protein_fenland')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
s<-jsonlite::fromJSON(file.path(root,'research/config/reconstruction_defaults.json'))$coloc
stopifnot(as.character(packageVersion('coloc'))=='5.2.3')
x<-read.delim(gzfile(file.path(input,'full_cis.tsv.gz')),stringsAsFactors=FALSE)
y<-read.delim(gzfile(file.path(input,'mvp_matched.ma.gz')),stringsAsFactors=FALSE)
names(y)<-c('SNP','EA_gwas','OA_gwas','freq_gwas','beta_gwas','se_gwas','p_gwas','n_gwas')
stopifnot(!anyDuplicated(y$SNP))
registry<-read.delim(file.path(input,'index_registry.tsv'),stringsAsFactors=FALSE)
indices<-read.delim(file.path(root,'research/data/processed/Fenland_Olink/source_indices.tsv'),stringsAsFactors=FALSE)
result<-list();counter<-0
for(i in seq_len(nrow(registry))){r<-registry[i,];for(window in c(s$index_window_kb,s$sensitivity_window_kb)){
  counter<-counter+1
  info<-list(gene_id=r$gene_id,symbol=r$symbol,uniprot=r$uniprot,assay_id=r$assay_id,window_kb=window,
             source='Fenland Olink 485 authorized substitute',sdY_estimated=r$sdY_estimated,analysis_sdY=1)
  if(r$status!='selected'){info$status<-r$status;result[[counter]]<-info;next}
  index<-indices[indices$id==r$assay_id,];stopifnot(nrow(index)==1)
  q<-x[x$assay_id==r$assay_id,];stopifnot(!anyDuplicated(q$SNP))
  info$index_snp<-index$rsid;info$index_pos<-index$pos;info$index_chr<-index$chr
  info$requested_start<-max(1,index$pos-window*1000);info$requested_end<-index$pos+window*1000
  info$source_min_pos<-min(q$pos);info$source_max_pos<-max(q$pos)
  info$source_registered_start<-index$region_start;info$source_registered_end<-index$region_end
  info$extends_beyond_registered_source<-info$requested_start<index$region_start | info$requested_end>index$region_end
  info$extends_beyond_observed_source<-info$requested_start<min(q$pos) | info$requested_end>max(q$pos)
  q<-q[q$pos>=info$requested_start & q$pos<=info$requested_end & q$chr==index$chr,];info$source_pairs_in_window<-nrow(q)
  d<-merge(q,y,by='SNP',sort=FALSE);info$shared_before_harmonization<-nrow(d)
  orientation<-lapply(seq_len(nrow(d)),function(j)orient_alleles(d$A1[j],d$A2[j],d$EA_gwas[j],d$OA_gwas[j],d$beta_gwas[j]))
  d$beta_gwas_aligned<-vapply(orientation,`[[`,numeric(1),'beta');d$harmonization<-vapply(orientation,`[[`,character(1),'status')
  d<-d[is.finite(d$beta_gwas_aligned)&is.finite(d$beta)&is.finite(d$se)&d$se>0 & is.finite(d$se_gwas)&d$se_gwas>0,]
  info$nsnps<-nrow(d)
  if(nrow(d)<2){info$status<-'insufficient_variants';result[[counter]]<-info;next}
  warnings<-character()
  fit<-withCallingHandlers(coloc::coloc.abf(
    list(snp=d$SNP,position=d$pos,type='quant',beta=d$beta,varbeta=d$se^2,sdY=1),
    list(snp=d$SNP,position=d$pos,type='cc',beta=d$beta_gwas_aligned,varbeta=d$se_gwas^2,s=15635/(15635+429252)),
    p1=s$p1,p2=s$p2,p12=s$p12),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  for(k in paste0('PP.H',0:4,'.abf'))info[[k]]<-unname(fit$summary[k])
  for(k in c('p1','p2','p12')){info[[paste0('requested_',k)]]<-s[[k]];info[[paste0('effective_',k)]]<-unname(fit$priors[k])}
  info$strong_coloc<-info[['PP.H4.abf']]>.8;info$status<-'computed';info$warnings<-paste(unique(warnings),collapse=' | ')
  d<-merge(d,fit$results[c('snp','SNP.PP.H4')],by.x='SNP',by.y='snp',sort=FALSE)
  name<-paste0(r$gene_id,'__',window,'kb.tsv.gz');info$variant_table<-name
  data.table::fwrite(d,file.path(out,name),sep='\t',compress='gzip');result[[counter]]<-info
}}
res<-data.table::rbindlist(result,fill=TRUE);data.table::fwrite(res,file.path(out,'coloc_summary.tsv'),sep='\t',na='NA')
pdf(file.path(out,'protein_H4.pdf'),width=8,height=6)
print(ggplot2::ggplot(res[status=='computed'],ggplot2::aes(x=factor(window_kb),y=symbol,fill=PP.H4.abf))+ggplot2::geom_tile()+ggplot2::scale_fill_viridis_c(limits=c(0,1))+ggplot2::theme_bw()+ggplot2::labs(x='Index window (kb)',y=NULL,title='Protein colocalization: Fenland Olink substitute'));dev.off()
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
inputs<-c(script,file.path(root,c('research/scripts/mr_core.R','research/config/reconstruction_defaults.json','research/data/processed/Fenland_Olink/source_indices.tsv')),file.path(input,c('run.json','index_registry.tsv','full_cis.tsv.gz','mvp_matched.ma.gz')))
outputs<-list.files(out,full.names=TRUE);outputs<-outputs[!grepl('/run.json$',outputs)]
record<-list(status='protein_coloc_public_substitute_computed',regions=nrow(res),computed=sum(res$status=='computed'),strong_H4_regions=sum(res$strong_coloc,na.rm=TRUE),
             original_data=FALSE,phenotype_scale='estimated sdY; raw beta/SE retained; estimator uncertainty not propagated',
             scope='Available follow-up plus author-reference genes; not deCODE; original gene-centred source window limits index-centred coverage',
             inputs=setNames(as.list(vapply(inputs,sha,character(1))),inputs),outputs=setNames(as.list(vapply(outputs,sha,character(1))),basename(outputs)))
jsonlite::write_json(record,file.path(out,'run.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA);print(record[c('status','regions','computed','strong_H4_regions')])
