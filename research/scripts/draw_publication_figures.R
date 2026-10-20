script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),'../..'))
script<-normalizePath(script)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(cowplot);library(ggrepel)})
cfg_path<-file.path(root,'research/config/publication_figures.json');cfg<-jsonlite::fromJSON(cfg_path)
C<-cfg$colors;FONT<-cfg$font;SEED<-cfg$seed;set.seed(SEED)
out<-file.path(root,'research/release');figdir<-file.path(out,'figures');datadir<-file.path(out,'figure_data')
dir.create(figdir,recursive=TRUE,showWarnings=FALSE);dir.create(datadir,recursive=TRUE,showWarnings=FALSE)
inputs<-c(script,cfg_path);panels<-list();figures<-list()
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
read_tab<-function(rel){
 p<-file.path(root,rel);inputs<<-unique(c(inputs,p));d<-read.delim(if(grepl('[.]gz$',p))gzfile(p) else p,stringsAsFactors=FALSE,check.names=FALSE)
 # pandas emits True/False while R emits TRUE/FALSE. Preserve missingness and
 # convert only columns whose observed values are exclusively boolean tokens.
 for(k in names(d))if(is.character(d[[k]])){
  values<-unique(na.omit(d[[k]]))
  if(length(values)&&all(tolower(values)%in%c('true','false')))d[[k]]<-as.logical(toupper(d[[k]]))
 }
 d
}
panel_data<-function(name,d){p<-file.path(datadir,paste0(name,'.tsv'));write.table(d,p,sep='\t',row.names=FALSE,quote=TRUE,na='NA');panels[[name]]<<-list(file=paste0('figure_data/',name,'.tsv'),rows=nrow(d),columns=names(d),sha256=sha(p));invisible(d)}
theme_pub<-function(){theme_minimal(base_size=cfg$base_font_pt,base_family=FONT)+theme(
 text=element_text(color=C$ink),panel.grid.minor=element_blank(),panel.grid.major=element_line(color=C$grid,linewidth=.22),
 axis.text=element_text(color='#465564',size=7),axis.title=element_text(size=8),axis.line.x=element_line(color='#BBC4CA',linewidth=.25),
 plot.title=element_text(size=10,face='bold',hjust=0),plot.subtitle=element_text(size=7.5,color='#52616D',margin=margin(b=5)),
 legend.title=element_text(size=7.5),legend.text=element_text(size=7),legend.key.height=grid::unit(3,'mm'),legend.key.width=grid::unit(4,'mm'),
 strip.text=element_text(size=8,face='bold',color=C$ink),strip.background=element_rect(fill='#F2F5F7',color=NA),
 plot.margin=margin(7,7,5,7),panel.spacing=grid::unit(6,'pt'))}
theme_set(theme_pub())
title_page<-function(p,title,subtitle=NULL){p+plot_annotation(title=title,subtitle=subtitle,tag_levels='a',theme=theme(plot.title=element_text(family=FONT,size=14,face='bold',color=C$ink),plot.subtitle=element_text(family=FONT,size=8,color='#52616D'),plot.tag=element_text(family=FONT,size=10,face='bold'),plot.margin=margin(10,9,8,9)))}
save_fig<-function(id,p,height,caption,table_names,width=cfg$width_inches){
 message('Drawing ',id)
 for(format in cfg$formats){path<-file.path(figdir,paste0(id,'.',format))
  if(format=='pdf')grDevices::cairo_pdf(path,width=width,height=height,family=FONT,onefile=TRUE)
  if(format=='svg')grDevices::svg(path,width=width,height=height,family=FONT,onefile=TRUE)
  if(format=='png')grDevices::png(path,width=width,height=height,units='in',res=cfg$png_dpi,type='cairo',bg='white')
  print(p);grDevices::dev.off()
 }
 files<-file.path(figdir,paste0(id,'.',cfg$formats))
 figures[[id]]<<-list(width_in=width,height_in=height,png_dpi=cfg$png_dpi,caption=caption,panel_tables=table_names,
    outputs=setNames(as.list(vapply(files,sha,character(1))),paste0('figures/',basename(files))))
}
logp<-function(x)-log10(pmax(x,1e-300))
targets<-read_tab('research/results/mr_fenland_substitute/followup_proteins.tsv')
targets<-targets[order(-targets$beta,targets$gene),]
genes<-targets$gene;uids<-targets$uniprot;gene_labels<-setNames(ifelse(genes=='TGFBI','TGFBI\u2020',genes),genes)
y_genes<-function()scale_y_discrete(limits=rev(genes),labels=gene_labels,drop=FALSE)
threshold<-.05/2709

# 01: preserve the reconstructed screen and show all three disease datasets.
allmr<-read_tab('research/results/mr_fenland_substitute/all_proteins_comparison.tsv')
allmr$followup<-allmr$uniprot %in% uids;allmr$display_p<-logp(allmr$p)
allmr$group<-ifelse(allmr$followup,'Follow-up set','Other tested proteins')
allmr$display_label<-allmr$gene %in% c('IL12B','ICAM5','DDX58','CTF1','CD8A','RARRES2','IL7R','TGFBI','NFKB1')
panel_data('F01a_discovery',allmr)
p1a<-ggplot(allmr,aes(beta,display_p))+geom_point(data=allmr[!allmr$followup,],color='#BAC4CB',size=.6,alpha=.65)+
 geom_point(data=allmr[allmr$followup,],aes(color=gene=='TGFBI'),size=1.6)+
 geom_hline(yintercept=-log10(threshold),linetype=2,linewidth=.35,color=C$orange)+
 geom_text_repel(data=allmr[allmr$followup & allmr$display_label,],aes(label=gene),size=2.45,family=FONT,seed=SEED,min.segment.length=0,max.overlaps=Inf,box.padding=.45,point.padding=.2,segment.size=.2)+
 scale_color_manual(values=c('FALSE'=C$teal,'TRUE'=C$orange),guide='none')+
 labs(x='Reconstructed pooled log odds ratio',y=expression(-log[10](italic(P))),title='Proteome-wide reconstruction',subtitle='2,213 estimable proteins; fixed testing family of 2,709')
forest<-do.call(rbind,lapply(c('mvp','finngen','stuart'),function(name){d<-read_tab(paste0('research/results/mr_fenland_substitute/',name,'_meta.tsv'));d<-d[d$uniprot %in% uids,];d$dataset<-name;d}))
forest$dataset<-factor(forest$dataset,levels=c('mvp','finngen','stuart'),labels=c('MVP discovery','FinnGen replication','Stuart replication'))
forest$conflict<-forest$gene=='TGFBI';panel_data('F01b_forest',forest)
forest_panel<-function(d,ys,title){ggplot(d,aes(beta,gene))+geom_vline(xintercept=0,color='#A9B4BC',linewidth=.35)+
 geom_errorbar(aes(xmin=beta-1.96*se,xmax=beta+1.96*se,color=dataset),orientation='y',width=.17,linewidth=.45)+
 geom_point(aes(fill=dataset,shape=conflict),size=1.8,color='white',stroke=.25)+
 scale_shape_manual(values=c('FALSE'=21,'TRUE'=23),guide='none')+
 scale_color_manual(values=c(C$ink,C$teal,C$purple),guide='none')+scale_fill_manual(values=c(C$ink,C$teal,C$purple),guide='none')+
 facet_grid(.~dataset)+scale_y_discrete(limits=rev(ys),labels=gene_labels)+labs(x='Log odds ratio (95% CI)',y=NULL,title=title)+theme(panel.grid.major.y=element_line(color='#EEF1F3'),panel.grid.major.x=element_blank())}
p1b<-forest_panel(forest[forest$gene!='CTF1',],genes[genes!='CTF1'],'Outcome-specific estimates')
p1c<-forest_panel(forest[forest$gene=='CTF1',],'CTF1','CTF1 on a separate axis range')+theme(strip.text=element_blank(),strip.background=element_blank())
save_fig('F01_genetic_reconstruction',title_page(p1a/p1b/p1c+plot_layout(heights=c(1,1.55,.34)),'Genetic associations across outcome datasets','Paper-style pooling shown for comparison; source dependence is examined in Figure 2.'),8.9,
 'a, All estimable proteins under the reconstructed independent-weight fixed-effect model. The dashed line is P=0.05/2709; labels identify selected diagnostic examples, while all follow-ups are shown below. b, Three outcome datasets for 16 reconstructed follow-up proteins; intervals are 95% normal-approximation confidence intervals on the standardized protein scale (estimated SD for Fenland). c, CTF1 is displayed separately with its full confidence intervals and a different x-axis range, avoiding compression of the other estimates. TGFBI\u2020 has a direction conflict in Stuart. Pooled effects reuse disease GWAS estimates across exposure sources and do not account for covariance; these are reconstruction outputs, not independence-validated causal findings. MVP SEs are inferred from reported confidence intervals.',c('F01a_discovery','F01b_forest'))

# 02: a separate fixed-single-IV sensitivity, not a silent replacement of pooling.
audit<-read_tab('research/results/shared_outcome_audit/all_proteins_comparison.tsv')
audit$pooled_logp<-logp(audit$p);audit$single_logp<-logp(audit$mvp_single_p)
audit$highlight<-ifelse(audit$gene %in% c('CD8A','RARRES2'),'Repeated outcome SNP',ifelse(audit$pooled_followup,'Other follow-up','Other protein'))
panel_data('F02a_paired_significance',audit)
p2a<-ggplot(audit,aes(pooled_logp,single_logp))+geom_abline(slope=1,intercept=0,color='#9EABB5',linewidth=.4)+
 geom_hline(yintercept=-log10(threshold),color=C$orange,linetype=2,linewidth=.3)+geom_vline(xintercept=-log10(threshold),color=C$orange,linetype=2,linewidth=.3)+
 geom_point(aes(color=highlight),size=.8,alpha=.65,na.rm=TRUE)+
 geom_text_repel(data=audit[audit$gene %in% c('CD8A','RARRES2','IL7R','IL12B'),],aes(label=gene),family=FONT,size=2.5,seed=SEED,box.padding=.4,min.segment.length=0,segment.size=.2)+
 scale_color_manual(values=c('Other protein'='#C7CED3','Other follow-up'=C$teal,'Repeated outcome SNP'=C$orange),guide='none')+
 labs(x=expression('Pooled reconstruction: '*-log[10](italic(P))),y=expression('Fixed single IV: '*-log[10](italic(P))),title='Dependence on pooling',subtitle='One exposure-only choice per protein; no outcome-based replacement')
f<-audit[audit$pooled_followup,];f$available<-is.finite(f$mvp_single_p);f$state<-ifelse(!f$available,'Not estimable',ifelse(f$single_IV_followup,'Retained under single IV','Not retained under single IV'))
stopifnot(nrow(f)==17,sum(f$available)==14)
panel_data('F02b_followup_sensitivity',f)
p2b<-ggplot(f,aes(y=gene))+geom_vline(xintercept=-log10(threshold),color=C$orange,linetype=2,linewidth=.35)+
 geom_segment(aes(x=pooled_logp,xend=single_logp,yend=gene),color='#B9C4CB',linewidth=.6,na.rm=TRUE)+
 geom_point(aes(x=pooled_logp),shape=1,size=1.7,color='#83929E',stroke=.5)+
 geom_point(aes(x=single_logp,color=state,shape=shared_same_SNP),size=2,na.rm=TRUE)+
 geom_text(data=f[!f$available,],aes(x=0,label='not estimable'),hjust=0,size=2.1,color='#7A8994',family=FONT)+
 scale_color_manual(values=c('Retained under single IV'=C$teal,'Not retained under single IV'=C$orange),name='Single-IV outcome')+
 scale_shape_manual(values=c('FALSE'=16,'TRUE'=18),labels=c('FALSE'='No repeated SNP','TRUE'='Repeated SNP'),name='Source indices')+y_genes()+
 labs(x=expression(-log[10](italic(P))),y=NULL,title='The 17-protein follow-up set',subtitle='Open points: pooled; filled points: fixed single IV')+
 theme(legend.position='bottom',legend.box='vertical',panel.grid.major.x=element_blank())
save_fig('F02_shared_outcome_sensitivity',title_page(p2a/p2b+plot_layout(heights=c(1,1.65)),'Sensitivity to shared outcome statistics'),8.1,
 'This extension freezes one instrument per protein using the largest exposure F statistic, then matches outcomes; missing choices are not replaced. a, Pooled and single-IV P values among estimable proteins. b, The original 17-protein follow-up set: 6 retain discovery plus nominal external support, 8 no longer pass, and 3 cannot be evaluated using the frozen choice. Diamonds identify repeated identical outcome SNPs. This sensitivity removes cross-source repeated-outcome weighting but is not a covariance-aware pooled estimator or an independently validated target list.',c('F02a_paired_significance','F02b_followup_sensitivity'))

# 03: every retained cell is shown; composition is summarized per donor.
emb<-read_tab('research/results/sc_embedding/cell_embeddings.tsv');ann<-read_tab('research/results/sc_de/annotated_cells.tsv')
stopifnot(!anyDuplicated(emb$cell_id),!anyDuplicated(ann$cell_id),setequal(emb$cell_id,ann$cell_id))
cells<-merge(emb,ann[c('cell_id','cell_type','condition')],by='cell_id',sort=FALSE)
cell_order<-c('Keratinocytes','Fibroblasts','Immune cells','Vascular endothelial cells','Mural cells','Melanocytes','Unassigned')
cell_colors<-setNames(c('#56B4E9','#E69F00','#CC79A7','#009E73','#D55E00','#0072B2','#CBD1D6'),cell_order)
cells$cell_type<-factor(cells$cell_type,levels=cell_order)
cells$condition<-factor(cells$condition,levels=c('normal','psoriasis'),labels=c('Normal skin','Psoriatic lesion'))
panel_data('F03a_cells',cells[c('cell_id','donor_id','condition','cell_type','UMAP_1','UMAP_2')])
centres<-aggregate(cbind(UMAP_1,UMAP_2)~cell_type,data=cells,FUN=median);centres<-centres[centres$cell_type!='Unassigned',]
centres$label<-c('Fibroblasts'='Fibroblasts','Immune cells'='Immune','Keratinocytes'='Keratinocytes','Melanocytes'='Melanocytes','Mural cells'='Mural','Vascular endothelial cells'='Endothelial')[as.character(centres$cell_type)]
atlas_theme<-theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.title=element_text(size=7),legend.position='bottom')
cells<-cells[sample.int(nrow(cells)),]
p3a<-ggplot(cells,aes(UMAP_1,UMAP_2,color=cell_type))+geom_point(size=.13,alpha=.7,stroke=0)+
 geom_label_repel(data=centres,aes(label=label),color=C$ink,fill=scales::alpha('white',.8),family=FONT,size=2.4,linewidth=0,seed=SEED,box.padding=.3,max.overlaps=Inf,min.segment.length=0,segment.size=.2)+
 scale_color_manual(values=cell_colors,drop=FALSE)+coord_equal()+labs(title='Skin cell atlas',subtitle='46,225 cells from 22 samples',x='UMAP 1',y='UMAP 2',color=NULL)+atlas_theme+guides(color=guide_legend(override.aes=list(size=2,alpha=1),nrow=3))
p3b<-ggplot(cells,aes(UMAP_1,UMAP_2))+geom_point(data=cells[c('UMAP_1','UMAP_2')],color='#E3E7EA',size=.1)+geom_point(aes(color=condition),size=.12,alpha=.65,stroke=0)+facet_wrap(~condition,ncol=1)+scale_color_manual(values=c(C$teal,C$orange),guide='none')+coord_equal()+labs(title='Shared embedding',x='UMAP 1',y='UMAP 2')+atlas_theme
counts<-as.data.frame(table(cells$donor_id,cells$condition,cells$cell_type));names(counts)<-c('donor','condition','cell_type','cells');counts<-counts[counts$cells>0,]
donor_meta<-unique(data.frame(donor=cells$donor_id,condition=cells$condition));donor_grid<-merge(donor_meta,data.frame(cell_type=cell_order),by=NULL)
counts<-merge(donor_grid,counts,by=c('donor','condition','cell_type'),all.x=TRUE);counts$cells[is.na(counts$cells)]<-0
counts$fraction<-counts$cells/ave(counts$cells,counts$donor,FUN=sum)
counts$donor<-factor(counts$donor,levels=unique(counts$donor[order(counts$condition,counts$donor)]));panel_data('F03c_donor_composition',counts)
p3c<-ggplot(counts,aes(donor,cell_type,fill=fraction))+geom_tile(width=.92,height=.86,color='white',linewidth=.15)+facet_grid(.~condition,scales='free_x',space='free_x')+
 scale_fill_gradient(low='#F2F6F7',high='#1E596E',trans='sqrt',limits=c(0,1),breaks=c(0,.1,.25,.5,1),labels=scales::percent,name='Cell fraction')+
 scale_y_discrete(limits=rev(cell_order),labels=c('Vascular endothelial cells'='Endothelial','Keratinocytes'='Keratinocytes','Fibroblasts'='Fibroblasts','Immune cells'='Immune','Mural cells'='Mural','Melanocytes'='Melanocytes','Unassigned'='Unassigned'))+
 labs(x='GEO subject label',y=NULL,title='Cell composition varies across samples')+theme(axis.text.x=element_text(angle=60,hjust=1,size=6.3),panel.grid=element_blank(),axis.line=element_blank(),legend.position='bottom')+
 guides(fill=guide_colorbar(barwidth=grid::unit(42,'mm'),barheight=grid::unit(2.5,'mm'),title.position='top'))
save_fig('F03_single_cell_atlas',title_page((p3a+p3b+plot_layout(widths=c(1.7,1)))/p3c+plot_layout(heights=c(2.2,1.2)),'Cellular context in psoriatic skin'),7.6,
 'a, UMAP of all 46,225 retained cells with broad reconstructed annotations; cells are shown in a fixed random drawing order, without subsampling. b, The same coordinates separated by condition; no independent re-embedding. c, Cell fractions for each of the 8 normal and 14 lesion samples, displayed on a square-root color scale with percentage labels. Columns use GEO subject labels. The source study reports 14 psoriasis patients; relationships between the labels 30696/30696V3 and 6215A/6215B are not separately explained in its sample table. Zero counts are observed zeros for known sample/type combinations; fractions describe sampled cells, not inferred population proportions. Unassigned cells are gray in the atlas and retained as a separate composition row. The six broad types do not recover the original paper\u2019s KC1/KC2 split.',c('F03a_cells','F03c_donor_composition'))

# 04: align tissue, cell-level and donor-level results without manufacturing zeros.
ev<-read_tab('research/results/evidence_substitute/followup_evidence.tsv')
cell_de<-read_tab('research/results/evidence_substitute/cell_Wilcoxon_followup_contrasts.tsv')
donor_de<-read_tab('research/results/evidence_substitute/donor_pseudobulk_followup_contrasts.tsv')
types<-cell_order[cell_order!='Unassigned'];short_types<-c('Keratinocytes'='Keratinocyte','Fibroblasts'='Fibroblast','Immune cells'='Immune','Vascular endothelial cells'='Endothelial','Mural cells'='Mural','Melanocytes'='Melanocyte')
expand_expression<-function(d,fold,fdr,layer){
 d$gene<-targets$gene[match(d$gene_id,targets$gene_id)];d<-d[c('gene','cell_type',fold,fdr,'significant','source_gene_id')]
 names(d)[3:4]<-c('log2FC','FDR');g<-expand.grid(gene=genes,cell_type=types,stringsAsFactors=FALSE)
 d<-merge(g,d,by=c('gene','cell_type'),all.x=TRUE,sort=FALSE);d$layer<-layer
 d$tested<-is.finite(d$FDR);d$display_FC<-ifelse(d$tested,pmax(-3,pmin(3,d$log2FC)),NA_real_)
 d$clipped<-d$tested & abs(d$log2FC)>3;d
}
cell_mat<-expand_expression(cell_de,'log2FC','fdr_bh','Cell-level Wilcoxon')
donor_mat<-expand_expression(donor_de,'logFC','FDR','Sample-level edgeR')
panel_data('F04a_bulk',ev[c('gene','bulk_log2FC','bulk_FDR','bulk_significant','bulk_tested')]);panel_data('F04b_cell',cell_mat);panel_data('F04c_donor',donor_mat)
p4a<-ggplot(ev,aes(bulk_log2FC,gene))+geom_vline(xintercept=0,color='#BBC4CA',linewidth=.3)+
 geom_segment(aes(x=0,xend=bulk_log2FC,yend=gene),color='#C4CDD3',linewidth=.6,na.rm=TRUE)+geom_point(aes(fill=bulk_significant),shape=21,color='white',stroke=.2,size=2,na.rm=TRUE)+
 geom_text(data=ev[!ev$bulk_tested,],aes(x=0,label='\u00d7'),size=2.6,color='#7E8992')+
 scale_fill_manual(values=c('FALSE'='#BCC6CE','TRUE'=C$ink),guide='none')+y_genes()+labs(x=expression(log[2]*' fold change'),y=NULL,title='Bulk RNA-seq',subtitle='90 lesions / 81 controls')+theme(panel.grid.major.x=element_blank())
expression_heat<-function(d,title,subtitle){d$sign_color<-ifelse(abs(d$display_FC)>1.8,'white',C$ink);ggplot(d,aes(cell_type,gene,fill=display_FC))+
 geom_tile(width=.92,height=.88,color='white',linewidth=.2)+geom_text(data=d[d$tested & !is.na(d$significant) & d$significant,],aes(label=ifelse(log2FC>0,'+','\u2212'),color=sign_color),size=2.1,fontface='bold')+scale_color_identity(guide='none')+
 geom_text(data=d[!d$tested,],label='\u00d7',color='#89959E',size=2)+
 geom_point(data=d[d$clipped,],shape=17,size=.55,color=C$ink,position=position_nudge(x=.27,y=.25))+
 scale_fill_gradient2(low=C$negative,mid='#F9FAFA',high=C$positive,midpoint=0,limits=cfg$expression_color_limits,na.value=C$missing,name=expression(log[2]*' FC'))+
 scale_x_discrete(limits=types,labels=short_types)+y_genes()+labs(x=NULL,y=NULL,title=title,subtitle=subtitle)+
 theme(panel.grid=element_blank(),axis.line=element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(angle=65,hjust=1,size=7),legend.position='bottom')+guides(fill=guide_colorbar(barwidth=grid::unit(35,'mm'),barheight=grid::unit(2.5,'mm'),title.position='top'))}
p4b<-expression_heat(cell_mat,'Cell-level comparison','All expressed genes; Wilcoxon')
p4c<-expression_heat(donor_mat,'Sample-level comparison','Pseudobulk; edgeR')
save_fig('F04_expression_evidence',title_page((p4a+p4b+p4c+plot_layout(widths=c(1.35,1.75,1.75),guides='collect'))&theme(legend.position='bottom'),'Expression evidence across biological scales','The same gene order is retained across panels; TGFBI\u2020 marks external genetic direction conflict.'),6.6,
 'a, Bulk fold changes from the 171-sample NCBI count reconstruction; dark points satisfy the recorded bulk FDR and fold-change criteria. b\u2013c, Disease contrasts in six broad cell types from cell-level Wilcoxon and sample-level edgeR, respectively; pseudobulk models include 6\u20138 normal and 12\u201314 lesion samples depending on cell type. GEO subject labels define the modeled units, following the source study\u2019s reported patient count; independence has not been separately verified at the person-identity level. Plus/minus signs indicate effect direction for comparisons satisfying FDR<0.05 and absolute log2FC>0.25; gray crosses denote untested/undefined comparisons, never zero effects. Colors are clipped to \u00b13, with corner triangles identifying clipped cells. Cell-level tests do not account for donor dependence. MMP12\u2019s version-specific gene identity mapping is retained in the source data.',c('F04a_bulk','F04b_cell','F04c_donor'))

# 05: paired-window posterior evidence, plus actual regional profiles.
mrna<-read_tab('research/results/coloc_mrna_substitute/coloc_summary.tsv');prot<-read_tab('research/results/coloc_protein_substitute/coloc_summary.tsv')
mrna$source_label<-c('eQTLGen_Blood'='eQTLGen\nblood','GTEx_Whole_Blood'='GTEx\nblood','GTEx_Skin_Sun_Exposed_Lower_leg'='GTEx skin\nexposed','GTEx_Skin_Not_Sun_Exposed_Suprapubic'='GTEx skin\nunexposed')[mrna$tissue]
prot$source_label<-'Plasma\nprotein'
posterior<-rbind(mrna[c('gene_id','symbol','window_kb','source_label','PP.H3.abf','PP.H4.abf','status','nsnps')],prot[c('gene_id','symbol','window_kb','source_label','PP.H3.abf','PP.H4.abf','status','nsnps')])
source_order<-c('Plasma\nprotein','eQTLGen\nblood','GTEx\nblood','GTEx skin\nexposed','GTEx skin\nunexposed')
posterior$source_label<-factor(posterior$source_label,levels=source_order);posterior$window<-factor(posterior$window_kb,levels=c(500,1000),labels=c('0.5','1.0'))
posterior$gene<-targets$gene[match(posterior$gene_id,targets$gene_id)]
selected_post<-posterior[!is.na(posterior$gene),];panel_data('F05a_coloc_windows',selected_post)
coloc_heat<-function(d,allgenes=genes){ggplot(d,aes(window,gene,fill=PP.H4.abf))+geom_tile(width=.94,height=.9,color='white',linewidth=.15)+
 geom_point(data=d[!is.na(d$PP.H4.abf)&d$PP.H4.abf>.8,],shape=21,size=1.1,fill='white',stroke=.3,color=C$ink)+
 geom_text(data=d[is.na(d$PP.H4.abf),],label='\u00d7',size=2,color='#89959E')+
 facet_grid(.~source_label)+scale_fill_gradientn(colors=c('#F3F6F7','#A6D0D3','#4A929C','#164E63'),values=c(0,.5,.8,1),limits=c(0,1),na.value=C$missing,name='PP.H4')+
 scale_y_discrete(limits=rev(allgenes),labels=function(x)ifelse(x=='TGFBI','TGFBI\u2020',x))+
 labs(x='Half-window around the source index (Mb)',y=NULL)+theme(panel.grid=element_blank(),axis.line=element_blank(),legend.position='right',strip.text=element_text(size=7.2),panel.spacing=grid::unit(4,'pt'))}
p5a<-coloc_heat(selected_post)+labs(title='Shared-variant posterior probability',subtitle='Open dots: PP.H4 > 0.8; crosses: no eligible source region/probe')
regional_plots<-list();regional_names<-character()
for(symbol in c('IL7R','HAVCR1')){
 row<-prot[prot$symbol==symbol&prot$window_kb==500,];stopifnot(nrow(row)==1)
 d<-read_tab(paste0('research/results/coloc_protein_substitute/',row$variant_table));best<-d$SNP[which.max(d$SNP.PP.H4)]
 z<-rbind(data.frame(SNP=d$SNP,position_Mb=d$pos/1e6,z=abs(d$beta/d$se),trait='Protein',SNP_PP_H4=d$SNP.PP.H4),data.frame(SNP=d$SNP,position_Mb=d$pos/1e6,z=abs(d$beta_gwas_aligned/d$se_gwas),trait='Psoriasis',SNP_PP_H4=d$SNP.PP.H4))
 z$trait<-factor(z$trait,levels=c('Protein','Psoriasis'));z$highlight<-z$SNP==best;z$genome_build<-'GRCh37'
 name<-paste0('F05_region_',symbol);panel_data(name,z);regional_names<-c(regional_names,name)
 regional_plots[[symbol]]<-ggplot(z,aes(position_Mb,z))+geom_vline(xintercept=d$pos[match(best,d$SNP)]/1e6,linetype=2,color='#A9B7C1',linewidth=.3)+
 geom_point(aes(color=trait),size=.45,alpha=.65)+geom_point(data=z[z$highlight,],shape=21,fill=C$orange,color='white',stroke=.2,size=1.8)+
 facet_wrap(~trait,ncol=1,scales='free_y')+scale_color_manual(values=c(C$teal,C$purple),guide='none')+
 labs(x=paste0('Chr ',row$index_chr,' (GRCh37), position (Mb)'),y='|z|',title=paste0(symbol,'  |  PP.H4 = ',formatC(row$PP.H4.abf,digits=2,format='f')),subtitle=paste0(nrow(d),' shared variants; highlighted ',best))+
 theme(strip.text=element_text(size=7),panel.grid.minor=element_blank(),plot.subtitle=element_text(size=6.5))
}
save_fig('F05_colocalization',title_page(p5a/(regional_plots[[1]]+regional_plots[[2]])+plot_layout(heights=c(1.55,1)),'Regional genetic evidence'),8.5,
 'a, PP.H4 for the 17 reconstructed follow-up proteins, retaining paired \u00b10.5 Mb and \u00b11 Mb windows and missing regions. The protein layer uses 485-person Fenland Olink substitute data; mRNA layers use eQTLGen and GTEx. b\u2013c, Available \u00b10.5 Mb protein/GWAS regions for IL7R and HAVCR1, plotted as absolute beta/SE to avoid conflating source P values with normal-approximation tests. Coordinates are GRCh37. Orange points mark the maximum conditional SNP PP.H4, not a proven causal variant. Coloc priors are p1=p2=1e-4, p12=1e-5; source cis coverage may truncate requested windows. These are single-causal-variant model results.',c('F05a_coloc_windows',regional_names))
fullpost<-posterior;fullpost$gene<-fullpost$symbol;all_genes<-sort(unique(fullpost$gene));panel_data('S01_all_coloc',fullpost)
save_fig('S01_all_registered_colocalization',title_page(coloc_heat(fullpost,all_genes),'Colocalization across all registered genes','The 33-gene registry includes the 17 reconstructed follow-ups and author-reference genes.'),9.2,
 'Full registered-gene counterpart to Figure 5a. Missing regions remain gray crosses; author-reference genes are not automatically nominated by the reconstructed discovery screen. Window labels are half-widths in Mb.',c('S01_all_coloc'))

# 06: effect magnitudes from the reconstruction, with transparent single-IV support.
clinical<-read_tab('research/results/shared_outcome_audit/clinical_comparison.tsv')
definition<-data.frame(phenocode=c('K11_CD_STRICT2','K11_IBD_STRICT','K11_UC_STRICT2','E4_HYPERLIPNAS','NAFLD','T2D','E4_OBESITY','I9_CHD','I9_HEARTFAIL','C_STROKE','I9_CORATHER','I9_HYPTENS','M13_RHEUMA','E4_THYROIDITAUTOIM','M13_SJOGREN','SLE_FG','M13_GOUT','E4_GRAVES_STRICT','C3_OTHER_SKIN_EXALLC','C3_MELANOMA_SKIN_EXALLC','L12_PSORIASIS','M13_PSORIARTH'),
 label=c('Crohn disease','IBD','Ulcerative colitis','Other hyperlipidaemia','NAFLD','Type 2 diabetes','Obesity','Major CHD events','Heart failure','Stroke','Coronary atherosclerosis','Hypertension','Rheumatoid arthritis','Autoimmune thyroiditis','Sjogren syndrome','SLE','Gout','Graves disease','Non-melanoma skin cancer','Melanoma','Psoriasis','Psoriatic arthropathies'),
 domain=c(rep('Intestinal',3),rep('Metabolic',4),rep('Cardiovascular',5),rep('Autoimmune / inflammatory',6),rep('Skin / joint',4)),stringsAsFactors=FALSE)
clinical<-merge(clinical,definition,by='phenocode',all.x=TRUE);stopifnot(!anyNA(clinical$domain))
clinical$phenocode<-factor(clinical$phenocode,levels=definition$phenocode);clinical$domain<-factor(clinical$domain,levels=unique(definition$domain))
clinical$display_beta<-pmax(-1,pmin(1,clinical$beta));clinical$clipped<-abs(clinical$beta)>1
clinical$support<-ifelse(is.na(clinical$single_p),'Single IV unavailable',ifelse(clinical$single_BH_all_requested_pairs<.05,'Single-IV BH < 0.05',ifelse(clinical$BH_all_requested_pairs<.05,'Pooled-only BH < 0.05','No displayed support')))
clinical$symbol_status<-ifelse(clinical$support=='Single IV unavailable','Single IV unavailable',ifelse(clinical$support=='Single-IV BH < 0.05',ifelse(clinical$single_beta>0,'Single IV: positive','Single IV: negative'),ifelse(clinical$support=='Pooled-only BH < 0.05',ifelse(clinical$beta>0,'Pooled only: positive','Pooled only: negative'),'None')))
panel_data('F06_clinical',clinical);panel_data('F06_endpoint_labels',definition)
p6<-ggplot(clinical,aes(phenocode,gene,fill=display_beta))+geom_tile(width=.96,height=.9,color='white',linewidth=.1)+
 geom_point(data=clinical[clinical$support!='No displayed support',],aes(shape=symbol_status),size=1.3,color=C$ink,fill=C$ink,stroke=.4)+
 geom_point(data=clinical[clinical$clipped,],shape=15,size=.45,color=C$ink,position=position_nudge(x=.28,y=.26))+
 facet_grid(.~domain,scales='free_x',space='free_x')+y_genes()+scale_x_discrete(labels=setNames(definition$label,definition$phenocode))+
 scale_fill_gradient2(low=C$negative,mid='#FBFCFC',high=C$positive,limits=cfg$clinical_color_limits,na.value=C$missing,name='Pooled log OR')+
 scale_shape_manual(values=c('Single IV: positive'=17,'Single IV: negative'=25,'Pooled only: positive'=2,'Pooled only: negative'=6,'Single IV unavailable'=4),breaks=c('Single IV: positive','Single IV: negative','Pooled only: positive','Pooled only: negative','Single IV unavailable'),name='BH < 0.05')+
 labs(x=NULL,y=NULL)+theme(panel.grid=element_blank(),axis.line=element_blank(),axis.text.x=element_text(angle=65,hjust=1,size=6.5),strip.text=element_text(size=6.6),panel.spacing=grid::unit(4,'pt'),legend.position='bottom',legend.box='vertical')+
 guides(fill=guide_colorbar(barwidth=grid::unit(35,'mm'),barheight=grid::unit(2.5,'mm'),title.position='top'),shape=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1.7)))
save_fig('F06_clinical_associations',title_page(p6,'Clinical associations and sensitivity to pooling','17 proteins across 22 FinnGen R12 phenotypes; all-pair multiple-testing family retained.'),7.4,
 'Tiles show reconstructed pooled log odds ratios, clipped to \u00b11; tiny corner squares identify clipped values. Filled up/down triangles denote positive/negative single-IV effects with BH<0.05; open triangles denote pooled-only BH<0.05 and its direction. Crosses denote an unavailable fixed single IV. The BH family remains all 374 pairs; 308 single-IV pairs are calculable and 31 satisfy BH<0.05, compared with 89 in the independent-weight pooled reconstruction. These differences also reflect exposure-source choice and coverage, not only covariance. Endpoint labels are abbreviated; full case/control definitions and sample sizes are in clinical_endpoints.tsv. I9_CHD and E4_HYPERLIPNAS are approved replacement phenotypes.',c('F06_clinical','F06_endpoint_labels'))

# 07: retain only actually significant enrichment terms; use a network matrix instead of a hairball.
enriched<-list()
for(lib in c('GO_biological_process','GO_molecular_function','GO_cellular_component','KEGG')){
 d<-read_tab(paste0('research/results/functional_substitute/',lib,'.tsv'));d<-d[d$fdr_bh<.05,];if(nrow(d)){d$library<-ifelse(lib=='KEGG','KEGG','GO biological process');enriched[[lib]]<-d}}
enriched<-do.call(rbind,enriched);enriched$description<-sub(' - Homo sapiens \\(human\\)$','',enriched$description)
enriched$label<-vapply(enriched$description,function(s)paste(strwrap(s,width=53),collapse='\n'),character(1));enriched$label<-factor(enriched$label,levels=rev(unique(enriched$label)))
panel_data('F07a_enrichment',enriched)
p7a<-ggplot(enriched,aes(-log10(fdr_bh),label))+geom_segment(aes(x=0,xend=-log10(fdr_bh),yend=label),linewidth=.7,color='#D4DDE2')+
 geom_point(aes(size=overlap,color=library),alpha=.95)+scale_color_manual(values=c('GO biological process'=C$teal,'KEGG'=C$purple),name=NULL)+
 scale_size_continuous(range=c(2.3,4.4),breaks=sort(unique(enriched$overlap)),name='Genes')+labs(x=expression(-log[10]*' BH FDR'),y=NULL,title='Enrichment in the reconstructed follow-up set',subtitle='Only terms with BH FDR < 0.05 are displayed')+theme(panel.grid.major.y=element_blank(),legend.position='bottom',legend.box='vertical')+guides(color=guide_legend(nrow=1,order=1),size=guide_legend(nrow=1,order=2))
network<-read_tab('research/results/drug_network_substitute/STRING_network.tsv');mapping<-read_tab('research/results/drug_network_substitute/STRING_mapping.tsv');relations<-read_tab('research/results/drug_network_substitute/all_drug_target_relations.tsv')
target_ids<-unique(relations$uniprot[relations$organism=='Homo sapiens']);mapped_targets<-mapping[mapping$queryItem %in% target_ids,]
mapped_targets<-mapped_targets[order(mapped_targets$preferredName),];target_ids<-mapped_targets$queryItem
grid<-expand.grid(uniprot=uids,target_uniprot=target_ids,stringsAsFactors=FALSE)
grid$gene<-targets$gene[match(grid$uniprot,targets$uniprot)];grid$target_label<-mapping$preferredName[match(grid$target_uniprot,mapping$queryItem)]
string_a<-mapping$stringId[match(grid$uniprot,mapping$queryItem)];string_b<-mapping$stringId[match(grid$target_uniprot,mapping$queryItem)]
edge_key<-paste(pmin(network$stringId_A,network$stringId_B),pmax(network$stringId_A,network$stringId_B),sep='::')
key<-paste(pmin(string_a,string_b),pmax(string_a,string_b),sep='::');grid$score<-network$score[match(key,edge_key)]
grid$component_identity<-grid$uniprot==grid$target_uniprot
grid$source_status<-ifelse(grid$component_identity,'same_protein_component_identity',ifelse(is.na(grid$score),'no_returned_edge_above_0.4','returned_functional_association'))
panel_data('F07b_network_matrix',grid)
p7b<-ggplot(grid,aes(target_label,gene,fill=score))+geom_tile(width=.95,height=.9,color='white',linewidth=.12)+
 geom_point(data=grid[grid$component_identity,],shape=23,size=1.5,color=C$purple,fill='white',stroke=.5)+
 y_genes()+scale_fill_gradientn(colors=c('#D9E9EC','#76A9B4','#244B61'),limits=c(.4,1),breaks=c(.4,.6,.8,1),na.value='#F3F5F7',name='STRING score')+
 labs(x='Proteins mapped to ChEMBL target records',y=NULL,title='Functional links to drug-target proteins',subtitle='No returned edge is left blank; diamonds mark identity with a listed complex component')+
 theme(panel.grid=element_blank(),axis.line=element_blank(),axis.text.x=element_text(angle=65,hjust=1,size=6.5),legend.position='bottom')+guides(fill=guide_colorbar(barwidth=grid::unit(38,'mm'),barheight=grid::unit(2.5,'mm'),title.position='top'))
save_fig('F07_functional_context',title_page(p7a/p7b+plot_layout(heights=c(1,2.6)),'Pathway and drug-target context'),8.5,
 'a, The two GO BP and two KEGG terms passing their separate BH families using the recorded tested-protein backgrounds. Many terms share genes and are not independent mechanisms; GO MF/CC have no FDR<0.05 terms. b, STRING12.5 functional association scores >0.4 between the 17 reconstructed follow-up proteins and mapped components of ChEMBL psoriasis-indication phase-4 target records. Blank cells mean no returned edge above this cutoff, not biological noninteraction. Diamonds mark the same protein listed as a complex component, not a STRING self-edge. There are no direct single-protein candidate drug-binding records in this snapshot; complex membership and functional-network links are distinct. VDR is unmapped and omitted from the matrix. These analyses are conditional on the original reconstructed 17-protein set, not a new causal or treatment validation.',c('F07a_enrichment','F07b_network_matrix'))

capture.output(sessionInfo(),file=file.path(out,'figure_sessionInfo.txt'))
manifest<-list(status='publication_figures_generated',created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),settings=cfg,
  scientific_scope='Paper-style reconstruction, explicitly separated post hoc single-IV sensitivity, observed expression and available regional evidence',
  inputs=setNames(as.list(vapply(inputs,sha,character(1))),substring(inputs,nchar(root)+2)),panel_tables=panels,figures=figures,
  session_info_sha256=sha(file.path(out,'figure_sessionInfo.txt')))
jsonlite::write_json(manifest,file.path(out,'figures_manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
captions<-unlist(lapply(names(figures),function(name)c(paste0('## ',name),'',figures[[name]]$caption,'',paste0('![',name,'](figures/',name,'.png)'),'')))
writeLines(c('# Figure captions','',captions),file.path(out,'FIGURE_CAPTIONS.md'))
