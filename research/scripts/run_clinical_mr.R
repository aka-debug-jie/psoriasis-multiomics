script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'../..'))
source(file.path(root,'research/scripts/mr_core.R'))
input<-file.path(root,'research/data/processed/clinical_finngen');out<-file.path(root,'research/results/clinical_mr_substitute');dir.create(out,recursive=TRUE,showWarnings=FALSE)
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
retrieval<-jsonlite::fromJSON(file.path(input,'run.json'));stopifnot(retrieval$status=='complete')
panel<-read.delim(file.path(input,'phenotype_manifest.tsv'),stringsAsFactors=FALSE)
follow<-read.delim(file.path(root,'research/results/mr_fenland_substitute/followup_proteins.tsv'),stringsAsFactors=FALSE)
x<-read.delim(file.path(root,'research/results/mr_fenland_substitute/mvp_cohort_wald.tsv'),stringsAsFactors=FALSE)
x<-x[x$uniprot %in% follow$uniprot,c('instrument_id','cohort','assay','uniprot','gene','rsid','EA_exposure','OA_exposure','beta_exposure','se_exposure')]
inputs<-c(script,file.path(root,'research/scripts/mr_core.R'),file.path(root,paste0('research/results/mr_fenland_substitute/',c('followup_proteins.tsv','mvp_cohort_wald.tsv'))),file.path(input,c('run.json','phenotype_manifest.tsv')))
cohorts<-list();estimates<-list();harm<-list()
for(i in seq_len(nrow(panel))){
  code<-panel$phenocode[i];path<-file.path(input,code,'regional_rows.tsv');inputs<-c(inputs,path,file.path(input,code,'retrieval.json'))
  raw<-read.delim(path,stringsAsFactors=FALSE,check.names=FALSE)
  y<-do.call(rbind,lapply(seq_len(nrow(raw)),function(j){ids<-strsplit(gsub(';',',',raw$rsids[j],fixed=TRUE),',',fixed=TRUE)[[1]];r<-raw[rep(j,length(ids)),];r$rsid<-ids;r}))
  y<-unique(y[c('rsid','alt','ref','beta','sebeta','pval','af_alt')]);bad<-y$rsid[duplicated(y$rsid)|duplicated(y$rsid,fromLast=TRUE)];y<-y[!y$rsid %in% bad,]
  d<-merge(x,y,by='rsid',all.x=TRUE)
  orientation<-lapply(seq_len(nrow(d)),function(j)orient_alleles(d$EA_exposure[j],d$OA_exposure[j],d$alt[j],d$ref[j],d$beta[j]))
  d$beta_aligned<-vapply(orientation,`[[`,numeric(1),'beta');d$harmonization<-vapply(orientation,`[[`,character(1),'status');d$phenocode<-code
  d$status<-ifelse(d$rsid %in% bad,'conflicting_outcome_rows',ifelse(is.na(d$alt),'exact_instrument_unavailable',d$harmonization));harm[[code]]<-d
  valid<-d[is.finite(d$beta_aligned)&is.finite(d$sebeta)&d$sebeta>0,]
  w<-wald_ratio(valid$beta_exposure,valid$se_exposure,valid$beta_aligned,valid$sebeta);names(w)[1:2]<-c('beta_mr','se_mr')
  valid<-cbind(valid,w);cohorts[[code]]<-valid
  m<-do.call(rbind,lapply(split(valid,valid$uniprot),function(g)cbind(uniprot=g$uniprot[1],fixed_effect(g$beta_mr,g$se_mr))))
  all<-merge(follow[c('uniprot','gene')],m,by='uniprot',all.x=TRUE);all$phenocode<-code;all$phenotype<-panel$phenotype[i]
  all$num_cases<-panel$num_cases[i];all$num_controls<-panel$num_controls[i]
  all$status<-ifelse(is.finite(all$beta),'computed_public_R12_reconstruction','no_harmonized_instrument')
  estimates[[code]]<-all
}
res<-data.table::rbindlist(estimates,fill=TRUE)
res$BH_all_requested_pairs<-p.adjust(res$p,method='BH',n=nrow(res))
res$paper_nominal_marker<-ifelse(is.na(res$p),NA,ifelse(res$p<.001,'***',ifelse(res$p<.01,'**',ifelse(res$p<.05,'*',''))))
res$direction_with_psoriasis<-ifelse(is.na(res$beta),NA,sign(res$beta)==sign(follow$beta[match(res$uniprot,follow$uniprot)]))
for(n in c('all_meta_results.tsv','cohort_wald.tsv','harmonization.tsv')){
  d<-switch(n,all_meta_results.tsv=res,cohort_wald.tsv=data.table::rbindlist(cohorts,fill=TRUE),harmonization.tsv=data.table::rbindlist(harm,fill=TRUE));data.table::fwrite(d,file.path(out,n),sep='\t',na='NA')
}
pdf(file.path(out,'clinical_associations.pdf'),width=14,height=9)
labels<-c(K11_CD_STRICT2='Crohn disease (strict)',K11_IBD_STRICT='IBD (strict)',K11_UC_STRICT2='Ulcerative colitis (strict)',E4_HYPERLIPNAS='Hyperlipidaemia\n(other/unspecified)',NAFLD='NAFLD',M13_RHEUMA='Rheumatoid arthritis',I9_CHD='Major CHD events',I9_HEARTFAIL='Heart failure (strict)',C_STROKE='Stroke',T2D='Type 2 diabetes',C3_OTHER_SKIN_EXALLC='Non-melanoma\nskin cancer',C3_MELANOMA_SKIN_EXALLC='Melanoma',I9_CORATHER='Coronary atherosclerosis',I9_HYPTENS='Hypertension',E4_THYROIDITAUTOIM='Autoimmune thyroiditis',E4_OBESITY='Obesity',M13_SJOGREN='Sjogren syndrome',SLE_FG='SLE',M13_GOUT='Gout',E4_GRAVES_STRICT='Graves disease (strict)',L12_PSORIASIS='Psoriasis',M13_PSORIARTH='Psoriatic arthropathies')
stopifnot(all(res$phenocode %in% names(labels)))
print(ggplot2::ggplot(res,ggplot2::aes(x=phenocode,y=gene,fill=beta))+ggplot2::geom_tile()+ggplot2::geom_text(ggplot2::aes(label=paper_nominal_marker),size=3,na.rm=TRUE)+ggplot2::scale_x_discrete(limits=names(labels),labels=labels)+ggplot2::scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',na.value='grey80')+ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1))+ggplot2::labs(x=NULL,y=NULL,fill='Log OR',title='FinnGen R12 clinical associations (substituted-source instruments)',subtitle='Stars: nominal P; all-pair BH and full endpoint definitions are saved in the accompanying tables'));dev.off()
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
outputs<-list.files(out,full.names=TRUE);outputs<-outputs[!grepl('/run.json$',outputs)]
record<-list(status='public_clinical_MR_computed',outcomes=nrow(panel),proteins=nrow(follow),requested_pairs=nrow(res),computed_pairs=sum(res$status=='computed_public_R12_reconstruction'),nominal_P05=sum(res$p<.05,na.rm=TRUE),BH05=sum(res$BH_all_requested_pairs<.05,na.rm=TRUE),
             original_panel_reproduction=FALSE,definition='Exact FinnGen R12 endpoint descriptions and counts retained; original paper did not report panel accessions; no neighboring-variant substitution',
             UKB_PGS_Cox_replaced=FALSE,inputs=setNames(as.list(vapply(inputs,sha,character(1))),inputs),outputs=setNames(as.list(vapply(outputs,sha,character(1))),basename(outputs)))
jsonlite::write_json(record,file.path(out,'run.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA);print(record[c('outcomes','proteins','requested_pairs','computed_pairs','nominal_P05','BH05')])
