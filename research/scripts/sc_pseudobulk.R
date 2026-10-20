# Independent donor-level extension, leaving the paper's cell-level method intact.
sc_pseudobulk <- function(root,config,output) {
  s<-config$pseudobulk
  tab<-data.table::fread(cmd=paste("gzip -dc",shQuote(file.path(root,s$counts))),data.table=FALSE,check.names=FALSE)
  counts<-as.matrix(tab[,-1]);rownames(counts)<-tab$gene_id
  samples<-read.delim(file.path(root,s$samples),stringsAsFactors=FALSE)
  features<-read.delim(file.path(root,s$features),stringsAsFactors=FALSE)
  stopifnot(identical(colnames(counts),samples$sample_id),identical(rownames(counts),features$gene_id),
            all(counts>=0),all(counts==floor(counts)),!anyDuplicated(samples[c("donor_id","cell_type")]))
  samples$included<-samples$n_cells>=s$min_cells_per_donor_type & samples$cell_type!="Unassigned"
  write_tsv<-function(d,name)data.table::fwrite(d,file.path(output,name),sep="\t",na="NA",quote=FALSE)
  write_tsv(samples,"sample_inclusion.tsv")
  results<-list();families<-list();norms<-list()
  for(type in sort(unique(samples$cell_type[samples$included]))){
    select<-samples$included & samples$cell_type==type
    meta<-samples[select,]
    group<-factor(meta$condition,levels=c("normal","psoriasis"))
    n<-table(group)
    if(any(n<s$min_donors_per_condition)){
      families[[type]]<-data.frame(cell_type=type,normal=unname(n[1]),psoriasis=unname(n[2]),tested_genes=0,significant_genes=NA_integer_,status="insufficient_donors")
      next
    }
    y<-edgeR::DGEList(counts[,select,drop=FALSE])
    design<-model.matrix(~group)
    keep<-edgeR::filterByExpr(y,design=design)
    y<-y[keep,,keep.lib.sizes=FALSE]
    y<-edgeR::calcNormFactors(y,method=s$normalization)
    y<-edgeR::estimateDisp(y,design,robust=s$robust)
    fit<-edgeR::glmQLFit(y,design,robust=s$robust)
    test<-edgeR::glmQLFTest(fit,coef="grouppsoriasis")
    d<-edgeR::topTags(test,n=Inf,sort.by="none")$table
    d$gene_id<-rownames(d);d$symbol<-features$symbol[match(d$gene_id,features$gene_id)]
    d$cell_type<-type;d$comparison<-"psoriasis_vs_normal";d$normal_donors<-unname(n[1]);d$psoriasis_donors<-unname(n[2])
    d$significant<-d$FDR<s$fdr & abs(d$logFC)>=s$abs_log2fc
    d$analysis_role<-s$analysis_role
    results[[type]]<-d
    families[[type]]<-data.frame(cell_type=type,normal=unname(n[1]),psoriasis=unname(n[2]),tested_genes=nrow(d),significant_genes=sum(d$significant),status="tested")
    norms[[type]]<-data.frame(sample_id=meta$sample_id,donor_id=meta$donor_id,cell_type=type,condition=meta$condition,n_cells=meta$n_cells,y$samples)
  }
  combined<-do.call(rbind,results);family<-do.call(rbind,families)
  write_tsv(combined,"donor_level_differential_expression.tsv")
  write_tsv(family,"comparison_summary.tsv")
  write_tsv(do.call(rbind,norms),"normalization.tsv")
  cell<-data.table::fread(file.path(root,s$cell_DE),data.table=FALSE)
  joint<-merge(cell[,c("gene_id","cell_type","log2FC","p_value","fdr_bh","significant")],
               combined[,c("gene_id","cell_type","logFC","PValue","FDR","significant")],
               by=c("gene_id","cell_type"),all.x=TRUE,suffixes=c("_cell","_donor"),sort=FALSE)
  joint$direction_agrees<-ifelse(is.na(joint$logFC),NA,sign(joint$log2FC)==sign(joint$logFC))
  joint$donor_result_missing_reason<-ifelse(is.na(joint$logFC),"low_expression_filter_or_insufficient_donors","")
  write_tsv(joint,"cell_vs_donor_comparison.tsv")
  pdf(file.path(output,"donor_level_volcanoes.pdf"),width=7,height=6)
  for(type in unique(combined$cell_type)){
    d<-combined[combined$cell_type==type,]
    plot(d$logFC,-log10(pmax(d$FDR,.Machine$double.xmin)),pch=16,cex=.35,col=ifelse(d$significant,"#216B85","#888888"),
         xlab="log2FC: psoriasis - normal",ylab="-log10 BH FDR",main=paste(type,"(donors)"))
    abline(v=c(-s$abs_log2fc,s$abs_log2fc),h=-log10(s$fdr),lty=2)
  }
  dev.off()
  list(status="donor_level_extension_completed",tested_cell_types=family$cell_type[family$status=="tested"],
       tested_gene_contrasts=nrow(combined),significant_gene_contrasts=sum(combined$significant),settings=s,
       versions=list(edgeR=as.character(packageVersion("edgeR")),limma=as.character(packageVersion("limma"))),
       generated_files=c("sample_inclusion.tsv","donor_level_differential_expression.tsv","comparison_summary.tsv",
                         "normalization.tsv","cell_vs_donor_comparison.tsv","donor_level_volcanoes.pdf"))
}
