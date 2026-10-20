# Complete cell-level comparisons on the frozen, broadly annotated reconstruction.
sc_de <- function(root, config, output) {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(presto))
  source(file.path(root,"research/scripts/sc_de_helpers.R"))
  s<-config$sc_de
  stopifnot(s$unknown_in_type_vs_rest, s$p_adjust=="BH_within_each_contrast_over_all_feature_rows",as.character(packageVersion("presto"))=="1.1.0")
  object<-readRDS(file.path(root,s$checkpoint))
  features<-read.delim(file.path(root,s$features),stringsAsFactors=FALSE)
  labels<-read.delim(file.path(root,s$annotation),stringsAsFactors=FALSE)
  cells<-object[[]]; cells$cell_id<-colnames(object); cells$cluster<-as.integer(as.character(Idents(object)))
  stopifnot(!anyDuplicated(labels$cluster),all(cells$cluster %in% labels$cluster))
  cells$cell_type<-labels$cell_type[match(cells$cluster,labels$cluster)]
  cells$condition<-ifelse(cells$tissue_type=="Normal skin","normal","psoriasis")
  stopifnot(all(cells$tissue_type %in% c("Normal skin","Psoriatic skin")),
            all(vapply(split(cells$condition,cells$donor_id),function(g)length(unique(g)),integer(1))==1))
  tested_types<-sort(unique(labels$cell_type[as.logical(labels$test_contrast)]))
  write_tsv<-function(d,name) data.table::fwrite(d,file.path(output,name),sep="\t",na="NA",quote=FALSE)
  write_tsv(cells,"annotated_cells.tsv")
  write_tsv(labels,"cluster_annotations.tsv")
  normalized<-LayerData(object,layer="data")
  stopifnot(identical(rownames(normalized),features$gene_id),ncol(normalized)==nrow(cells))
  message("Testing all ",nrow(normalized)," genes: type versus all other cells, including Unassigned")
  all_wilcox<-presto::wilcoxauc(normalized,cells$cell_type)
  enrichment<-list();disease<-list();family<-list()
  for(type in tested_types){
    message("Cell-level contrasts: ",type)
    enrichment[[type]]<-sc_contrast_table(normalized,cells$cell_type,type,type,"type_vs_rest",features,all_wilcox)
    idx<-cells$cell_type==type
    disease[[type]]<-sc_contrast_table(normalized[,idx,drop=FALSE],cells$condition[idx],"psoriasis",type,"psoriasis_vs_normal",features)
    family[[type]]<-data.frame(cell_type=type,cells=sum(idx),normal=sum(idx & cells$condition=="normal"),
                              psoriasis=sum(idx & cells$condition=="psoriasis"),
                              normal_donors=length(unique(cells$donor_id[idx & cells$condition=="normal"])),
                              psoriasis_donors=length(unique(cells$donor_id[idx & cells$condition=="psoriasis"])),
                              tested_gene_family=nrow(features))
    gc()
  }
  enrichment<-do.call(rbind,enrichment);disease<-do.call(rbind,disease);family<-do.call(rbind,family)
  for(name in c("enrichment","disease")){
    d<-get(name)
    d$significant<-ifelse(is.na(d$fdr_bh),NA,d$fdr_bh<s$fdr & abs(d$log2FC)>=s$abs_log2fc)
    assign(name,d)
  }
  write_tsv(enrichment,"all_gene_type_vs_rest.tsv")
  write_tsv(disease,"all_gene_psoriasis_vs_normal.tsv")
  write_tsv(family,"contrast_sample_sizes.tsv")
  # Author targets are a separate comparison panel, never a reconstructed candidate list.
  author<-read.delim(file.path(root,s$author_targets),stringsAsFactors=FALSE)
  target_mapping<-do.call(rbind,lapply(author$target,function(target){
    symbol<-switch(target,IL23="IL23A",IL12="IL12A",DDX58="RIGI",target)
    hits<-features$gene_id[!is.na(features$symbol) & features$symbol==symbol]
    data.frame(author_target=target,requested_symbol=symbol,gene_id=if(length(hits)==1)hits else NA_character_,
               mapping_status=if(length(hits)==1)"unique" else if(length(hits)==0)"missing" else "ambiguous")
  }))
  write_tsv(target_mapping,"author_target_mapping.tsv")
  full<-rbind(enrichment,disease)
  panel<-merge(expand.grid(author_target=target_mapping$author_target,cell_type=tested_types,
                           comparison=c("type_vs_rest","psoriasis_vs_normal"),stringsAsFactors=FALSE),target_mapping,by="author_target",sort=FALSE)
  panel<-merge(panel,full,by=c("gene_id","cell_type","comparison"),all.x=TRUE,sort=FALSE)
  panel$author_panel_fdr_bh<-NA_real_
  for(comparison in unique(panel$comparison)){
    idx<-panel$comparison==comparison
    panel$author_panel_fdr_bh[idx]<-p.adjust(panel$p_value[idx],"BH",n=length(author$target)*length(tested_types))
  }
  panel$evidence_origin<-"recomputed_author_reference_panel_not_formal_discovery"
  write_tsv(panel,"author_reference_all_contrasts.tsv")
  # Pseudobulk retains raw counts and biological donors for a separate extension.
  keys<-paste(cells$donor_id,cells$cell_type,sep="::")
  groups<-sort(unique(keys));membership<-match(keys,groups)
  design<-Matrix::sparseMatrix(i=seq_len(nrow(cells)),j=membership,x=1,dims=c(nrow(cells),length(groups)))
  counts<-LayerData(object,layer="counts") %*% design
  colnames(counts)<-sprintf("PB%03d",seq_along(groups));rownames(counts)<-features$gene_id
  meta<-do.call(rbind,lapply(seq_along(groups),function(i){
    members<-cells[membership==i,]
    data.frame(sample_id=colnames(counts)[i],donor_id=members$donor_id[1],cell_type=members$cell_type[1],
               condition=members$condition[1],n_cells=nrow(members),gsm=paste(unique(members$orig.ident),collapse=";"),
               eligible_cells=nrow(members)>=s$pseudobulk_min_cells)
  }))
  # GSM sample names are kept explicitly; orig.ident is not donor identification.
  meta$gsm<-vapply(seq_along(groups),function(i)paste(unique(cells$gsm[membership==i]),collapse=";"),character(1))
  write_tsv(meta,"pseudobulk_samples.tsv")
  pb<-data.frame(gene_id=rownames(counts),as.matrix(counts),check.names=FALSE)
  data.table::fwrite(pb,file.path(output,"pseudobulk_raw_counts.tsv.gz"),sep="\t",compress="gzip",quote=FALSE)
  stopifnot(sum(counts)==sum(LayerData(object,layer="counts")))
  object$cell_type<-cells$cell_type
  pdf(file.path(output,"umap_broad_cell_types.pdf"),width=10,height=7)
  print(DimPlot(object,group.by="cell_type",label=TRUE,raster=TRUE))
  dev.off()
  pdf(file.path(output,"cell_level_disease_volcanoes.pdf"),width=7,height=6)
  for(type in tested_types){
    d<-disease[disease$cell_type==type & is.finite(disease$p_value),]
    plot(d$log2FC,-log10(pmax(d$fdr_bh,.Machine$double.xmin)),pch=16,cex=.3,
         col=ifelse(d$significant,"#B23A48","#888888"),xlab="log2FC: psoriasis - normal",ylab="-log10 BH FDR",main=type)
    abline(v=c(-1,1)*s$abs_log2fc,h=-log10(s$fdr),lty=2,col="grey40")
  }
  dev.off()
  summary<-do.call(rbind,lapply(tested_types,function(type){
    en<-enrichment[enrichment$cell_type==type,];de<-disease[disease$cell_type==type,]
    data.frame(cell_type=type,type_enriched=sum(en$significant & en$log2FC>0,na.rm=TRUE),
               disease_significant=sum(de$significant,na.rm=TRUE),
               noninformative_genes=sum(is.na(de$p_value)))
  }))
  write_tsv(summary,"comparison_summary.tsv")
  list(status="reconstructed_single_cell_DE_completed",cells=nrow(cells),genes=nrow(features),
       tested_cell_types=tested_types,unassigned_cells=sum(cells$cell_type=="Unassigned"),
       type_vs_rest_rows=nrow(enrichment),disease_comparison_rows=nrow(disease),
       pseudobulk_samples=nrow(meta),formal_candidates_used=FALSE,
       test="Presto 1.1.0 asymptotic Wilcoxon; ties and continuity corrected",settings=s,
       limitations=c("Broad reconstructed labels; paper Keratinocytes 1/2 not recovered",
                     "Mural label includes pericyte/smooth muscle-like cells",
                     "Cell-level tests follow the paper but do not model donor dependence",
                     "Full-gene BH differs from unspecified original candidate testing family"),
       generated_files=c("annotated_cells.tsv","cluster_annotations.tsv","all_gene_type_vs_rest.tsv",
                         "all_gene_psoriasis_vs_normal.tsv","contrast_sample_sizes.tsv","author_target_mapping.tsv",
                         "author_reference_all_contrasts.tsv","pseudobulk_samples.tsv","pseudobulk_raw_counts.tsv.gz",
                         "umap_broad_cell_types.pdf","cell_level_disease_volcanoes.pdf","comparison_summary.tsv"))
}
