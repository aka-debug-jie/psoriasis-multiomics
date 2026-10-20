script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'../..'))
source(file.path(root,'research/scripts/mr_core.R'))
source(file.path(root,'research/tools/TwoSampleMR-methods/R/add_rsq.R'))
input<-file.path(root,'research/results/mr_fenland_substitute');out<-file.path(root,'research/results/steiger_substitute')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
d<-read.delim(file.path(input,'mvp_cohort_wald.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
replications<-lapply(c('finngen','stuart'),function(o)read.delim(file.path(input,paste0(o,'_cohort_wald.tsv')),stringsAsFactors=FALSE));names(replications)<-c('finngen','stuart')
n<-c(Iceland=35892,'UKB-PPP'=46218,Fenland_Olink_substitute_485=485)
d$n_exposure_cohort_approximation<-unname(n[d$cohort])
d$R2_exposure_standardized_effect<-2*d$eaf_exposure*(1-d$eaf_exposure)*d$beta_exposure^2
d$R2_exposure_bse_N<-(d$beta_exposure/d$se_exposure)^2/((d$beta_exposure/d$se_exposure)^2+d$n_exposure_cohort_approximation-2)
all<-list();summaries<-list();followups<-list();counter<-0
threshold<-jsonlite::fromJSON(file.path(root,'research/data/processed/substitute_mr/selection.json'))$bonferroni_threshold
combine<-function(x)do.call(rbind,lapply(split(x,x$uniprot),function(g)cbind(uniprot=g$uniprot[1],gene=g$gene[1],fixed_effect(g$beta_mr,g$se_mr))))
for(K in c(.01,.02,.03,.05)){
  # MVP EA frequency and beta remain in their original, mutually matched orientation.
  ry<-get_r_from_lor(d$beta_outcome,d$eaf_outcome,15635,429252,K,model='logit',correction=FALSE)
  for(estimator in c('standardized_effect','bse_N')){
    counter<-counter+1
    r<-d[c('instrument_id','cohort','uniprot','gene','rsid','n_exposure_cohort_approximation')]
    r$prevalence<-K;r$exposure_R2_estimator<-estimator
    r$R2_exposure<-d[[paste0('R2_exposure_',estimator)]];r$R2_outcome<-ry^2
    valid<-is.finite(r$R2_exposure)&r$R2_exposure>=0&r$R2_exposure<1&is.finite(ry)&abs(ry)<1
    r$forward_direction<-ifelse(valid,r$R2_exposure>r$R2_outcome,NA)
    z<-(atanh(sqrt(r$R2_exposure[valid]))-atanh(abs(ry[valid])))/sqrt(1/(r$n_exposure_cohort_approximation[valid]-3)+1/(15635+429252-3))
    r$approximate_direction_p<-NA_real_;r$approximate_direction_p[valid]<-2*pnorm(-abs(z))
    r$status<-ifelse(valid,'computed_under_declared_approximations','invalid_R2_not_filtered_as_negative')
    all[[counter]]<-r
    keep<-r$instrument_id[!is.na(r$forward_direction)&r$forward_direction]
    m<-combine(d[d$instrument_id %in% keep,])
    for(o in names(replications)){
      v<-replications[[o]];v<-combine(v[v$instrument_id %in% keep,]);v<-v[c('uniprot','p','beta')];names(v)[2:3]<-paste0(c('p_','beta_'),o)
      m<-merge(m,v,by='uniprot',all.x=TRUE)
    }
    m$external_support<-(is.finite(m$p_finngen)&m$p_finngen<.05)|(is.finite(m$p_stuart)&m$p_stuart<.05)
    m$external_direction_agreement<-(is.finite(m$p_finngen)&m$p_finngen<.05&m$beta*m$beta_finngen>0)|(is.finite(m$p_stuart)&m$p_stuart<.05&m$beta*m$beta_stuart>0)
    m$selected<-m$p<threshold&m$external_support;m$prevalence<-K;m$exposure_R2_estimator<-estimator
    followups[[counter]]<-m
    summaries[[counter]]<-data.frame(prevalence=K,exposure_R2_estimator=estimator,cohort_rows=nrow(r),forward=sum(r$forward_direction,na.rm=TRUE),reverse=sum(!r$forward_direction,na.rm=TRUE),invalid=sum(!valid),discovery=sum(m$p<threshold),externally_supported=sum(m$selected))
  }
}
write_tsv<-function(x,name)data.table::fwrite(data.table::rbindlist(x,fill=TRUE),file.path(out,name),sep='\t',na='NA')
write_tsv(all,'instrument_direction_grid.tsv');write_tsv(summaries,'sensitivity_summary.tsv');write_tsv(followups,'all_meta_grid.tsv')
sel<-data.table::rbindlist(followups);registry<-read.delim(file.path(root,'research/data/processed/substitute_mr/protein_registry.tsv'))
robust<-sel[selected==TRUE,.(scenarios_selected=.N,scenarios_with_direction_agreement=sum(external_direction_agreement)),by=uniprot];robust<-merge(robust,registry,by='uniprot');robust$selected_all_eight<-robust$scenarios_selected==8
robust$direction_conflict_in_any_selected_scenario<-robust$scenarios_with_direction_agreement<robust$scenarios_selected
data.table::fwrite(robust,file.path(out,'candidate_robustness.tsv'),sep='\t',na='NA')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
inputs<-c(script,file.path(root,c('research/scripts/mr_core.R','research/tools/TwoSampleMR-methods/R/add_rsq.R','research/data/processed/substitute_mr/selection.json','research/data/processed/substitute_mr/protein_registry.tsv')),file.path(input,paste0(c('mvp','finngen','stuart'),'_cohort_wald.tsv')))
outputs<-list.files(out,full.names=TRUE);outputs<-outputs[!grepl('/run.json$',outputs)]
record<-list(status='Steiger_parameter_sensitivity_computed',primary_prevalence_chosen=FALSE,prevalence_grid=c(.01,.02,.03,.05),
             exposure_sample_size_approximation=as.list(n),exposure_N_note='Eldjarn Table1 normalized Iceland / UKB-BI totals; per-assay N unavailable; Fenland reported N=485',
             direction_filter='R2_exposure>R2_outcome; approximate direction P reported separately, not a significance gate',
             binary_method='Pinned TwoSampleMR get_r_from_lor logit correction=FALSE; case fraction never substituted for population prevalence',
             sources=list(Eldjarn='https://doi.org/10.1038/s41586-023-06563-x',method='https://raw.githubusercontent.com/MRCIEU/TwoSampleMR/d8b07182da23fb984157d9335fefcb83e894a944/R/add_rsq.r'),
             inputs=setNames(as.list(vapply(inputs,sha,character(1))),inputs),outputs=setNames(as.list(vapply(outputs,sha,character(1))),basename(outputs)))
jsonlite::write_json(record,file.path(out,'run.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA);print(data.table::rbindlist(summaries))
