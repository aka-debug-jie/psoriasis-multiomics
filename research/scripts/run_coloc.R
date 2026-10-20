args<-commandArgs(trailingOnly=TRUE)
script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),"../.."))
value<-function(key,default){i<-match(key,args);if(is.na(i))default else args[i+1]}
config_path<-value("--config",file.path(root,"research/config/coloc_mrna.json"))
s<-jsonlite::fromJSON(config_path)
input<-file.path(root,s$input_directory);output<-value("--output",file.path(root,s$output_directory))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
sha<-function(path)strsplit(system2("sha256sum",shQuote(path),stdout=TRUE)," +")[[1]][1]
source_manifest<-jsonlite::fromJSON(file.path(input,"run.json"))
stopifnot(as.character(packageVersion("coloc"))=="5.2.3")
approved<-jsonlite::fromJSON(file.path(root,s$reconstruction_config))$coloc
stopifnot(s$p1==approved$p1,s$p2==approved$p2,s$p12==approved$p12,
          identical(as.numeric(s$windows_kb),as.numeric(c(approved$index_window_kb,approved$sensitivity_window_kb))))
# Derive the observed sample fraction from actual outcome counts, not prevalence.
s$case_fraction<-s$MVP_cases/(s$MVP_cases+s$MVP_controls)
started<-format(Sys.time(),tz="UTC",usetz=TRUE)
jsonlite::write_json(list(status="running",started_utc=started),file.path(output,"run.json"),auto_unbox=TRUE,pretty=TRUE)
y<-data.table::fread(file.path(input,"mvp_regional.tsv"),data.table=FALSE)
names(y)<-c("SNP","A1_gwas","A2_gwas","freq_gwas","beta_gwas","se_gwas","p_gwas","n_gwas")
stopifnot(!anyDuplicated(y$SNP))
inputs<-c(config_path,script,file.path(root,s$reconstruction_config),file.path(input,c("run.json","gene_registry.tsv","mvp_regional.tsv")))
results<-list();generated<-character();counter<-0L
comp<-c(A="T",T="A",C="G",G="C")
for(tissue in source_manifest$tissues$tissue){
  message("Coloc mRNA: ",tissue)
  raw_path<-file.path(input,tissue,"full_cis.txt");idx_path<-file.path(input,tissue,"index_registry.tsv")
  inputs<-c(inputs,raw_path,idx_path)
  qtls<-data.table::fread(raw_path,data.table=FALSE)
  indices<-read.delim(idx_path,stringsAsFactors=FALSE)
  for(i in seq_len(nrow(indices))){
    idx<-indices[i,]
    for(window in s$windows_kb){
      counter<-counter+1L
      record<-list(tissue=tissue,gene_id=idx$gene_id,symbol=idx$symbol,ProbeID=idx$ProbeID,
                   available_followup=idx$available_followup,author_reference=idx$author_reference,
                   window_kb=window,index_snp=idx$topSNP,index_chr=idx$topSNP_chr,index_bp=idx$topSNP_bp,
                   requested_p1=s$p1,requested_p2=s$p2,requested_p12=s$p12,
                   sdY=s$quant_sdY,case_fraction=s$case_fraction)
      if(is.na(idx$ProbeID) || is.na(idx$topSNP)){
        record$status<-idx$index_status;results[[counter]]<-record;next
      }
      q<-qtls[qtls$Probe==idx$ProbeID & qtls$Chr==idx$topSNP_chr,]
      record$source_cis_pairs<-nrow(q)
      if(!nrow(q)){record$status<-"no_source_pairs";results[[counter]]<-record;next}
      start<-max(1,idx$topSNP_bp-window*1000);end<-idx$topSNP_bp+window*1000
      record$requested_start<-start;record$requested_end<-end
      record$source_min_bp<-min(q$BP);record$source_max_bp<-max(q$BP)
      record$extends_beyond_observed_source_span<-start<min(q$BP) || end>max(q$BP)
      q<-q[q$BP>=start & q$BP<=end,];record$source_pairs_in_window<-nrow(q)
      q<-unique(q)
      bad<-q$SNP[duplicated(q$SNP) | duplicated(q$SNP,fromLast=TRUE)]
      record$conflicting_qtl_rsids<-length(unique(bad));q<-q[!q$SNP %in% bad,]
      d<-merge(q,y,by="SNP",sort=FALSE);record$shared_before_harmonization<-nrow(d)
      single<-d$A1 %in% names(comp) & d$A2 %in% names(comp)
      palindromic<-paste0(d$A1,d$A2) %in% c("AT","TA","CG","GC")
      direct<-d$A1==d$A1_gwas & d$A2==d$A2_gwas
      swap<-d$A1==d$A2_gwas & d$A2==d$A1_gwas
      complement<-unname(comp[d$A1])==d$A1_gwas & unname(comp[d$A2])==d$A2_gwas
      flipcomp<-unname(comp[d$A1])==d$A2_gwas & unname(comp[d$A2])==d$A1_gwas
      ok<-single & !palindromic & (direct|swap|complement|flipcomp) & is.finite(d$b) & is.finite(d$SE) & d$SE>0 & is.finite(d$beta_gwas) & is.finite(d$se_gwas) & d$se_gwas>0
      ok[is.na(ok)]<-FALSE
      record$palindromic_excluded<-sum(palindromic)
      d$beta_gwas_aligned<-ifelse(swap|flipcomp,-d$beta_gwas,d$beta_gwas)
      d<-d[ok,];d<-d[order(d$BP,d$SNP),]
      record$nsnps<-nrow(d)
      if(nrow(d)<2){record$status<-"insufficient_shared_variants";results[[counter]]<-record;next}
      d1<-list(snp=d$SNP,position=d$BP,type="quant",beta=d$b,varbeta=d$SE^2,sdY=s$quant_sdY)
      d2<-list(snp=d$SNP,position=d$BP,type="cc",beta=d$beta_gwas_aligned,varbeta=d$se_gwas^2,s=s$case_fraction)
      fit<-tryCatch(coloc::coloc.abf(d1,d2,p1=s$p1,p2=s$p2,p12=s$p12),error=function(e)e)
      if(inherits(fit,"error")){record$status<-"model_failed";record$error<-conditionMessage(fit);results[[counter]]<-record;next}
      for(key in paste0("PP.H",0:4,".abf"))record[[key]]<-unname(fit$summary[key])
      record$effective_p1<-unname(fit$priors['p1']);record$effective_p2<-unname(fit$priors['p2']);record$effective_p12<-unname(fit$priors['p12'])
      record$strong_coloc<-record[['PP.H4.abf']]>.8;record$moderate_coloc<-record[['PP.H4.abf']]>.5
      record$status<-"computed_full_available_cis_region"
      aligned<-merge(d,fit$results[,c("snp","SNP.PP.H4")],by.x="SNP",by.y="snp",sort=FALSE)
      name<-paste0(tissue,"__",idx$ProbeID,"__",window,"kb.tsv.gz")
      data.table::fwrite(aligned,file.path(output,name),sep="\t",compress="gzip",na="NA")
      record$variant_table<-name;generated<-c(generated,name);results[[counter]]<-record
    }
  }
}
summary<-data.table::rbindlist(results,fill=TRUE)
data.table::fwrite(summary,file.path(output,"coloc_summary.tsv"),sep="\t",na="NA")
generated<-c(generated,"coloc_summary.tsv")
completed<-summary[status=="computed_full_available_cis_region"]
if(nrow(completed)){
  plotdata<-as.data.frame(completed)
  pdf(file.path(output,"mRNA_coloc_H4.pdf"),width=11,height=10)
  print(ggplot2::ggplot(plotdata,ggplot2::aes(x=tissue,y=symbol,fill=PP.H4.abf))+ggplot2::geom_tile()+
          ggplot2::facet_wrap(~window_kb)+ggplot2::scale_fill_viridis_c(limits=c(0,1),na.value="grey80")+
          ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))+
          ggplot2::labs(x=NULL,y=NULL,fill="PP.H4",title="mRNA colocalization on available cis regions"))
  dev.off();generated<-c(generated,"mRNA_coloc_H4.pdf")
}
capture.output(sessionInfo(),file=file.path(output,"sessionInfo.txt"));generated<-c(generated,"sessionInfo.txt")
record<-list(status=if(any(summary$status=="model_failed"))"completed_with_model_failures" else "mRNA_colocalization_computed",
             started_utc=started,finished_utc=format(Sys.time(),tz="UTC",usetz=TRUE),settings=s,
             coloc=as.character(packageVersion("coloc")),regions=nrow(summary),computed=nrow(completed),
             strong_H4_regions=sum(completed$strong_coloc),formal_candidate_list_used=FALSE,
             inputs=setNames(as.list(vapply(inputs,sha,character(1))),inputs),
             outputs=setNames(as.list(vapply(file.path(output,generated),sha,character(1))),generated))
jsonlite::write_json(record,file.path(output,"run.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
print(record[c("status","regions","computed","strong_H4_regions")])
