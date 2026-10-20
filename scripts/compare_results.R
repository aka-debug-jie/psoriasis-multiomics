# Same inputs and parameters; joins use biological identifiers rather than row order.
new<-Sys.getenv('TEACHING_NOTEBOOK_OUTPUT')
specs<-list(
 list('01_bulk','all_gene_differential_expression.tsv','GeneID',c('logFC','AveExpr','t','P.Value','adj.P.Val','B')),
 list('03_proteome_mr','all_outcomes.tsv',c('uniprot','outcome'),c('beta','se','p','Q')),
 list('04_sensitivity','multi_instrument_methods.tsv',c('cohort','uniprot','mode','method'),c('b','se','pval')),
 list('05_smr','two_probe.smr','probeID',c('b_SMR','se_SMR','p_SMR','p_HEIDI')),
 list('06_colocalization','coloc_summary.tsv',c('gene','window_kb','p12'),paste0('PP.H',0:4,'.abf')),
 list('07_single_cell','donor_pseudobulk.tsv',c('gene_id','cell_type'),c('logFC','logCPM','F','PValue','FDR')),
 list('10_functional','enrichment_all_families.tsv',c('family','term'),c('overlap','query_annotated','term_background','background_annotated','p_value','fdr_bh')),
 list('11_clinical_mr','clinical_candidates_154.tsv',c('uniprot','phenocode'),c('beta','se','p','BH_original_family_374')))
checks<-list();k<-0L
for(s in specs) {
 a<-read.delim(file.path('results',s[[1]],s[[2]]),check.names=FALSE)
 b<-read.delim(file.path(new,s[[1]],s[[2]]),check.names=FALSE)
 key<-function(x)do.call(paste,c(x[s[[3]]],sep='\r'))
 stopifnot(!anyDuplicated(key(a)),!anyDuplicated(key(b)),setequal(key(a),key(b)))
 b<-b[match(key(a),key(b)),,drop=FALSE]
 grouping<-if(s[[1]]=='10_functional')'family' else if(s[[1]]=='06_colocalization')'gene' else NULL
 groups<-if(is.null(grouping))'all' else unique(a[[grouping]])
 for(g in groups)for(v in s[[4]]) {
  use<-if(is.null(grouping))rep(TRUE,nrow(a)) else a[[grouping]]==g
  av<-a[[v]][use];bv<-b[[v]][use];stopifnot(identical(is.na(av),is.na(bv)))
  ok<-is.na(av)|abs(av-bv)<=1e-10+1e-8*abs(av)
  k<-k+1L;checks[[k]]<-data.frame(module=s[[1]],group=g,column=v,passed=all(ok),max_abs_error=max(abs(av-bv),na.rm=TRUE))
 }

}
checks<-do.call(rbind,checks);write.table(checks,file.path(new,'numeric_comparison.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(checks);stopifnot(all(checks$passed))

# Stochastic embeddings may move; retained cells, annotations and test decisions are checked separately.
sc<-function(base,name)read.delim(file.path(base,'07_single_cell',name),check.names=FALSE)
a<-sc('results','cell_qc.tsv');b<-sc(new,'cell_qc.tsv');stopifnot(setequal(a$cell_id,b$cell_id));b<-b[match(a$cell_id,b$cell_id),]
qc_same<-identical(a$keep,b$keep)&&identical(a$donor_id,b$donor_id)&&identical(a$condition,b$condition)
a<-sc('results','umap.tsv');b<-sc(new,'umap.tsv');stopifnot(setequal(a$cell_id,b$cell_id));b<-b[match(a$cell_id,b$cell_id),]
pairs<-unique(data.frame(reference=a$cluster,current=b$cluster))
partition_same<-!anyDuplicated(pairs$reference)&&!anyDuplicated(pairs$current)
types_same<-identical(a$cell_type,b$cell_type)
sc_checks<-list(data.frame(check=c('retained_cells_and_donors','cluster_partition','cell_types'),passed=c(qc_same,partition_same,types_same),max_abs_error=NA_real_))
for(name in c('candidate_cell_wilcoxon.tsv','candidate_type_vs_rest_wilcoxon.tsv')) {
 a<-sc('results',name);b<-sc(new,name);key<-function(x)paste(x$gene,x$cell_type,sep='\r')
 stopifnot(!anyDuplicated(key(a)),!anyDuplicated(key(b)),setequal(key(a),key(b)));b<-b[match(key(a),key(b)),]
 for(column in c('p','BH')) {
  av<-a[[column]];bv<-b[[column]];same_na<-identical(is.na(av),is.na(bv));finite<-is.finite(av)&is.finite(bv)
  sc_checks[[length(sc_checks)+1L]]<-data.frame(check=paste(name,column,sep=':'),passed=same_na&&identical(av<.05,bv<.05),max_abs_error=if(any(finite))max(abs(av[finite]-bv[finite]))else NA_real_)
 }
}
sc_checks<-do.call(rbind,sc_checks)
write.table(sc_checks,file.path(new,'single_cell_comparison.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(sc_checks);stopifnot(all(sc_checks$passed))
