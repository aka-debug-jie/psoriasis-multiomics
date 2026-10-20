script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'../..'))
source(file.path(root,'research/scripts/mr_core.R'))
methoddir<-file.path(root,'research/tools/TwoSampleMR-methods/R')
source(file.path(methoddir,'mr.R'));source(file.path(methoddir,'mr_mode.R'))
args<-commandArgs(trailingOnly=TRUE);substituted<-'--substitute' %in% args
input<-file.path(root,if(substituted)'research/data/processed/sensitivity_substitute' else 'research/data/processed/sensitivity_inputs')
out<-file.path(root,if(substituted)'research/results/multi_instrument_substitute' else 'research/results/multi_instrument_available')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(file.path(out,'clumps'),showWarnings=FALSE)
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
write_tsv<-function(d,n)write.table(d,file.path(out,n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
x<-read.delim(file.path(input,'all_loci.tsv'),stringsAsFactors=FALSE);y<-read.delim(file.path(input,'mvp_subset.tsv'),stringsAsFactors=FALSE)
x$se_for_analysis<-if('se_used' %in% names(x))x$se_used else x$se_reconstructed
x<-x[x$neg_log10_p > -log10(5e-6),]
x<-do.call(rbind,lapply(seq_len(nrow(x)),function(i){ids<-strsplit(x$variant[i],',',fixed=TRUE)[[1]];ids<-ids[grepl('^rs[0-9]+$',ids)];if(!length(ids))return(NULL);r<-x[rep(i,length(ids)),];r$rsid<-ids;r}))
plink<-file.path(root,'research/tools/plink-1.9.0/plink');ref<-file.path(root,'research/data/processed/LD/EUR')
stopifnot(file.exists(plink),!anyDuplicated(y$SNP))
writeLines(sort(unique(x$rsid)),file.path(out,'requested_LD_rsids.txt'))
small<-file.path(out,'LD_subset')
status<-system2(plink,c('--bfile',shQuote(ref),'--extract',shQuote(file.path(out,'requested_LD_rsids.txt')),'--make-bed','--out',shQuote(small),'--threads','4'),stdout=file.path(out,'LD_extract.log'),stderr=file.path(out,'plink_stderr.log'))
stopifnot(status==0)
bim<-read.table(paste0(small,'.bim'),stringsAsFactors=FALSE)
x$in_LD_reference<-x$rsid %in% bim$V2
write_tsv(x,'source_loci_with_LD_coverage.tsv')
result<-list();retained<-list();coverage<-list();counter<-0L
groups<-split(x,interaction(x$cohort,x$uniprot,drop=TRUE))
parameters<-default_parameters();parameters$nboot<-1000L
for(g in groups){for(mode in c('cis','cis_trans')){
  counter<-counter+1L
  z<-if(mode=='cis')g[g$locus_type=='cis',] else g
  id<-paste(g$cohort[1],g$uniprot[1],mode,sep='__');prefix<-file.path(out,'clumps',id)
  cov<-data.frame(cohort=g$cohort[1],uniprot=g$uniprot[1],gene=g$gene[1],mode=mode,source_rows=nrow(z),LD_present=sum(z$in_LD_reference),clumped=0,harmonized=0)
  if(g$cohort[1]=='Fenland_Olink_substitute_485' && mode=='cis_trans'){cov$status<-'trans_source_unavailable';coverage[[id]]<-cov;next}
  z<-z[z$in_LD_reference,];z<-z[order(-z$neg_log10_p,z$rsid),];z<-z[!duplicated(z$rsid),]
  if(!nrow(z)){cov$status<-'no_LD_variants';coverage[[id]]<-cov;next}
  # Real P threshold was applied above. Strict monotonic rank values preserve
  # exposure significance order for PLINK without underflow at log10 P > 300.
  clumpinput<-data.frame(SNP=z$rsid,P=seq_len(nrow(z))/(nrow(z)+1))
  write.table(clumpinput,paste0(prefix,'.input.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  status<-system2(plink,c('--bfile',shQuote(small),'--clump',shQuote(paste0(prefix,'.input.tsv')),'--clump-p1','1','--clump-p2','1','--clump-r2','0.1','--clump-kb','10000','--out',shQuote(prefix),'--threads','4'),stdout=paste0(prefix,'.stdout.log'),stderr=file.path(out,'plink_stderr.log'))
  stopifnot(status==0)
  clumped<-read.table(paste0(prefix,'.clumped'),header=TRUE,stringsAsFactors=FALSE)
  z<-z[z$rsid %in% clumped$SNP,];cov$clumped<-nrow(z)
  d<-merge(z,y,by.x='rsid',by.y='SNP',suffixes=c('_exposure','_outcome'))
  orientation<-lapply(seq_len(nrow(d)),function(i)orient_alleles(d$EA[i],d$OA[i],d$A1[i],d$A2[i],d$b[i]))
  d$beta_outcome_aligned<-vapply(orientation,`[[`,numeric(1),'beta')
  d$harmonization<-vapply(orientation,`[[`,character(1),'status')
  d<-d[is.finite(d$beta_outcome_aligned)&is.finite(d$se)&d$se>0,]
  cov$harmonized<-nrow(d);cov$status<-if(nrow(d))'computed' else 'no_harmonized_variants';coverage[[id]]<-cov
  if(!nrow(d))next
  d$mode<-mode;retained[[id]]<-d
  for(method in c('IVW','weighted_median','weighted_mode','MR_Egger')){
    r<-data.frame(cohort=g$cohort[1],uniprot=g$uniprot[1],gene=g$gene[1],mode=mode,method=method,nsnp=nrow(d),seed=20261017L+counter)
    if(nrow(d)<3 && method!='IVW'){r$status<-'requires_at_least_three_IVs';result[[paste(id,method)]]<-r;next}
    set.seed(r$seed)
    if(nrow(d)==1){w<-wald_ratio(d$beta,d$se_for_analysis,d$beta_outcome_aligned,d$se);fit<-list(b=w$beta,se=w$se,pval=w$p);r$method<-'single_IV_Wald'}
    else{
      fn<-switch(method,IVW=mr_ivw,weighted_median=mr_weighted_median,weighted_mode=mr_weighted_mode,MR_Egger=mr_egger_regression)
      fit<-tryCatch(fn(d$beta,d$beta_outcome_aligned,d$se_for_analysis,d$se,parameters=parameters),error=function(e)e)
    }
    if(inherits(fit,'error')){r$status<-'estimator_failed';r$error<-conditionMessage(fit)}
    else{for(k in intersect(c('b','se','pval','Q','Q_df','Q_pval','b_i','se_i','pval_i'),names(fit)))r[[k]]<-fit[[k]];r$status<-if(is.finite(fit$b))'computed' else 'estimator_unavailable'}
    result[[paste(id,method)]]<-r
  }
}}
res<-data.table::rbindlist(result,fill=TRUE);cov<-data.table::rbindlist(coverage,fill=TRUE)
write_tsv(res,'method_results.tsv');write_tsv(cov,'coverage.tsv');write_tsv(data.table::rbindlist(retained,fill=TRUE),'harmonized_instruments.tsv')
pdf(file.path(out,'sensitivity_estimates.pdf'),width=10,height=8)
for(selected_mode in c('cis','cis_trans')){plotdata<-res[res$status=='computed' & res$mode==selected_mode,];if(nrow(plotdata))print(ggplot2::ggplot(plotdata,ggplot2::aes(x=gene,y=b,color=method))+ggplot2::geom_pointrange(ggplot2::aes(ymin=b-1.96*se,ymax=b+1.96*se),position=ggplot2::position_dodge(width=.7))+ggplot2::facet_wrap(~cohort,scales='free')+ggplot2::coord_flip()+ggplot2::theme_bw()+ggplot2::labs(title=paste('Published-locus sensitivity:',selected_mode),y='Log odds ratio',x=NULL))}
dev.off()
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
inputs<-c(script,if(substituted)file.path(input,'run.json') else file.path(root,'research/records/advance_20261017/sensitivity_inputs_run.json'),file.path(input,c('all_loci.tsv','mvp_subset.tsv')),file.path(methoddir,c('mr.R','mr_mode.R')),plink,paste0(ref,c('.bed','.bim','.fam')))
outputs<-list.files(out,recursive=TRUE,full.names=TRUE);outputs<-outputs[!grepl('/run.json$',outputs)]
record<-list(status='multi_instrument_available_sources_computed',scope=if(substituted)'17 substituted-source follow-up proteins; Eldjarn published loci and Fenland full cis regions; Fenland trans absent' else '12 two-source follow-up proteins; published conditional/sentinel loci, not complete genomewide P<5e-6 coverage',groups=nrow(cov),computed_method_rows=sum(res$status=='computed'),IVW='multiplicative random effects with underdispersion correction (TwoSampleMR mr_ivw)',LD=list(reference='1000G Phase3 EUR 503 individuals MAF>0.01',r2=.1,kb=10000,exposure_p=5e-6),parameters=parameters,rank_surrogate='Actual exposure P filter first; PLINK order from strict exposure significance ranks; p1=p2=1',
             inputs=setNames(as.list(vapply(inputs,sha,character(1))),inputs),outputs=setNames(as.list(vapply(outputs,sha,character(1))),substring(outputs,nchar(out)+2)))
jsonlite::write_json(record,file.path(out,'run.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA);print(cov)
