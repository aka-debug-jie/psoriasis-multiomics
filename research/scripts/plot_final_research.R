script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'../..'))
script<-normalizePath(script)
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
args<-commandArgs(trailingOnly=TRUE);out<-if(length(args))file.path(root,args[1]) else file.path(root,'research/final')
fig<-file.path(out,'figures');panels<-file.path(out,'figure_data');dir.create(fig,showWarnings=FALSE);dir.create(panels,showWarnings=FALSE)
inputs<-c(script);records<-list();panel_records<-list()
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
read<-function(name){p<-file.path(out,'data',paste0(name,'.tsv'));inputs<<-unique(c(inputs,p));d<-read.delim(p,stringsAsFactors=FALSE,check.names=FALSE)
 for(k in names(d))if(is.character(d[[k]])){v<-unique(na.omit(d[[k]]));if(length(v)&&all(tolower(v)%in%c('true','false')))d[[k]]<-as.logical(toupper(d[[k]]))};d}
panel<-function(name,d){p<-file.path(panels,paste0(name,'.tsv'));write.table(d,p,sep='\t',row.names=FALSE,quote=TRUE,na='NA');panel_records[[name]]<<-list(file=paste0('figure_data/',name,'.tsv'),rows=nrow(d),sha256=sha(p))}
C<-list(ink='#243746',teal='#28788A',orange='#CE7849',purple='#806B9B',blue='#2B6B89',missing='#EDF1F3')
theme_set(theme_minimal(base_family='Liberation Sans',base_size=8)+theme(text=element_text(color=C$ink),
 panel.grid.minor=element_blank(),panel.grid.major=element_line(color='#E8EEF1',linewidth=.25),
 axis.text=element_text(size=7,color='#435563'),axis.title=element_text(size=8),strip.text=element_text(face='bold',size=8),
 strip.background=element_rect(fill='#F2F5F7',color=NA),plot.title=element_text(face='bold',size=10),
 plot.subtitle=element_text(size=7,color='#52616D'),legend.position='bottom',legend.title=element_text(size=7),legend.text=element_text(size=7),
 panel.spacing=grid::unit(5,'pt'),plot.margin=margin(8,8,6,8)))
page<-function(p,title,subtitle=NULL)p+plot_annotation(title=title,subtitle=subtitle,tag_levels='a',
 theme=theme(plot.title=element_text(family='Liberation Sans',face='bold',size=14,color=C$ink),plot.subtitle=element_text(family='Liberation Sans',size=8),plot.tag=element_text(size=10,face='bold')))
save<-function(id,p,height,caption){caption<-paste(caption,'Gene markers: dagger indicates DDX58 residual-dispersion sensitivity; double dagger indicates IL7R strict-LD sensitivity.');for(ext in c('pdf','svg','png')){
 path<-file.path(fig,paste0(id,'.',ext));if(ext=='pdf')cairo_pdf(path,width=7.2,height=height,family='Liberation Sans')
 if(ext=='svg')svg(path,width=7.2,height=height,family='Liberation Sans')
 if(ext=='png')png(path,width=7.2,height=height,units='in',res=600,type='cairo',bg='white')
 print(p);dev.off()}
 paths<-file.path(fig,paste0(id,c('.pdf','.svg','.png')))
 records[[id]]<<-list(width_in=7.2,height_in=height,png_dpi=600,caption=caption,outputs=setNames(as.list(vapply(paths,sha,character(1))),paste0('figures/',basename(paths))))}
e<-read('integrated_evidence');genes<-e$gene;labels<-setNames(ifelse(genes=='DDX58','DDX58\u2020',ifelse(genes=='IL7R','IL7R\u2021',genes)),genes)
ys<-function()scale_y_discrete(limits=rev(genes),labels=labels,drop=FALSE)
logp<-function(p)-log10(pmax(p,1e-300))

# 1. Primary genetic associations and explicit sensitivity.
f<-read('genetic_effects');f$outcome<-factor(f$outcome,levels=c('mvp','finngen','stuart'),labels=c('MVP: CI-derived SE','FinnGen: reported SE','Stuart: reported SE'));panel('F01a_genetic_effects',f)
p1a<-ggplot(f,aes(beta,gene))+geom_vline(xintercept=0,color='#A8B5BE',linewidth=.35)+geom_segment(aes(x=CI_log_low,xend=CI_log_high,yend=gene),color=C$teal,linewidth=.6)+geom_point(color=C$teal,size=1.8)+
 facet_grid(.~outcome)+ys()+labs(x='Log odds ratio (95% CI)',y=NULL,title='Primary and external genetic associations',subtitle='One outcome observation per variant; signed LD included')+theme(panel.grid.major.x=element_blank())
s<-read('primary_sensitivity');s<-s[(s$branch=='unique_variant_givw'&s$scenario=='CI_Wald')|(s$branch=='ld_pruned_givw'&s$scenario=='CI_Wald'&s$model=='fixed_effect')|(s$branch=='unique_variant_givw'&s$scenario=='P_equivalent'&s$model=='fixed_effect'),]
s$method<-ifelse(s$scenario=='P_equivalent','P-equivalent',ifelse(s$branch=='ld_pruned_givw','Strict LD',ifelse(s$model=='multiplicative_random_effect','Residual dispersion','Primary')))
s<-s[c('gene','method','p','followup')];s<-rbind(s,data.frame(gene=e$gene,method='FinnGen-led control',p=NA_real_,followup=e$MVP_free_FinnGen_discovery&e$MVP_free_Stuart_concordant))
s$method<-factor(s$method,levels=c('Primary','Residual dispersion','Strict LD','P-equivalent','FinnGen-led control'))
s$status<-ifelse(is.na(s$followup),'Unavailable',ifelse(s$followup,'Meets stated rule','Does not meet rule'));panel('F01b_sensitivity',s)
p1b<-ggplot(s,aes(method,gene,fill=status))+geom_tile(width=.92,height=.86,color='white',linewidth=.2)+ys()+
 geom_text(aes(label=ifelse(status=='Meets stated rule','+',ifelse(status=='Unavailable','\u00d7',''))),size=3,color='white')+
 scale_fill_manual(values=c('Meets stated rule'=C$teal,'Does not meet rule'='#DAE4E9','Unavailable'='#9DAAB3'),name=NULL)+
 labs(x=NULL,y=NULL,title='Sensitivity is part of the result',subtitle='\u2020 DDX58: dispersion-sensitive.  \u2021 IL7R: LD-sensitive.')+
 theme(panel.grid=element_blank(),axis.text.x=element_text(angle=25,hjust=1))
save('F01_genetic_evidence',page(p1a/p1b+plot_layout(heights=c(1,1.25)),'Genetic evidence for the primary candidate set','The main analysis uses CI-derived MVP SE; source P/CI convention remains unconfirmed.'),7.1,
 'a, Fixed-effect generalized IVW and both external outcomes, with one GWAS estimate per unique variant and reference LD. b, Each column applies its stated discovery and external-support rule. The FinnGen-led column uses no MVP statistics but reuses existing outcome datasets. Symbols identify DDX58 sensitivity to residual dispersion and IL7R sensitivity to strict LD. P-equivalent analysis adds ICAM3 outside this fixed primary set. These markers are not an evidence score or a claim of validated causal targets.')

# 2. Whole-study expression tests and evidence coverage, with no seven-gene FDR shortcut.
types<-c('Keratinocytes','Fibroblasts','Immune cells','Vascular endothelial cells','Mural cells','Melanocytes')
short<-c('Keratinocytes'='Keratinocyte','Fibroblasts'='Fibroblast','Immune cells'='Immune','Vascular endothelial cells'='Endothelial','Mural cells'='Mural','Melanocytes'='Melanocyte')
expression_data<-function(name,fc,fdr){d<-read(name);d$gene<-e$gene[match(d$gene_id,e$gene_id)];d<-d[c('gene','cell_type',fc,fdr)];names(d)[3:4]<-c('logFC','FDR')
 d<-merge(expand.grid(gene=genes,cell_type=types),d,by=c('gene','cell_type'),all.x=TRUE);d$tested<-is.finite(d$FDR);d$display<-ifelse(d$tested,pmin(3,pmax(-3,d$logFC)),NA_real_);d$clipped<-d$tested&abs(d$logFC)>3;d}
ce<-expression_data('expression_cell','log2FC','fdr_bh');pb<-expression_data('expression_pseudobulk','logFC','FDR');panel('F02b_cell_expression',ce);panel('F02c_pseudobulk',pb)
heat<-function(d,title){ggplot(d,aes(cell_type,gene,fill=display))+geom_tile(width=.92,height=.85,color='white')+
 geom_text(data=d[d$tested&d$FDR<.05&abs(d$logFC)>.25,],aes(label=ifelse(logFC>0,'+','\u2212')),size=2.6,color=C$ink)+
 geom_text(data=d[!d$tested,],label='\u00d7',size=2.4,color='#8C9AA4')+geom_point(data=d[d$clipped,],shape=17,size=.7,position=position_nudge(x=.25,y=.2))+
 scale_fill_gradient2(low=C$blue,mid='white',high=C$orange,limits=c(-3,3),na.value=C$missing,name=expression(log[2]*' FC'))+
 ys()+scale_x_discrete(limits=types,labels=short)+labs(x=NULL,y=NULL,title=title)+theme(panel.grid=element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(angle=60,hjust=1))}
bulk<-e[c('gene','bulk_log2FC','bulk_FDR','bulk_significant')];panel('F02a_bulk',bulk)
p2a<-ggplot(bulk,aes(bulk_log2FC,gene))+geom_vline(xintercept=0,color='#A8B5BE',linewidth=.3)+geom_segment(aes(x=0,xend=bulk_log2FC,yend=gene),color='#C4D0D7')+geom_point(aes(fill=bulk_significant),shape=21,size=2.3,color='white')+
 scale_fill_manual(values=c('FALSE'='#B8C5CD','TRUE'=C$ink),guide='none')+ys()+labs(x=expression(log[2]*' FC'),y=NULL,title='Bulk RNA-seq')+theme(panel.grid.major.x=element_blank())
matrix<-read('evidence_matrix');component_labels<-c('Bulk'='Bulk','Cell Wilcoxon'='Cells','Pseudobulk'='Pseudobulk',
 'GTEx_Skin_Not_Sun_Exposed_Suprapubic'='Skin: unexposed','GTEx_Skin_Sun_Exposed_Lower_leg'='Skin: exposed','GTEx_Whole_Blood'='Blood: GTEx','eQTLGen_Blood'='Blood: eQTLGen','coloc_mRNA'='mRNA','coloc_protein'='Protein')
matrix$component<-factor(matrix$component,levels=names(component_labels));matrix$domain<-factor(matrix$domain,levels=c('Expression','SMR / HEIDI','Colocalization'));panel('F02d_evidence_coverage',matrix)
p2d<-ggplot(matrix,aes(component,gene,fill=status))+geom_tile(width=.92,height=.85,color='white')+ys()+facet_grid(.~domain,scales='free_x',space='free_x')+
 scale_x_discrete(labels=component_labels)+scale_fill_manual(values=c(supported=C$teal,not_supported_in_tested_data='#D8E2E7',partial_no_observed_support='#B8A9CC',not_evaluable='white'),
 labels=c(supported='Observed support',not_supported_in_tested_data='No support in tested data',partial_no_observed_support='Partial coverage, no observed support',not_evaluable='Unavailable'),name=NULL)+
 geom_text(data=matrix[matrix$status=='not_evaluable',],label='\u00d7',color='#87949E',size=3)+
 labs(x=NULL,y=NULL,title='Evidence availability and observed support',subtitle='SMR and colocalization reuse MVP statistics; columns are not independent votes.')+
 theme(panel.grid=element_blank(),axis.text.x=element_text(angle=35,hjust=1),legend.text=element_text(size=6.5))+guides(fill=guide_legend(nrow=2,override.aes=list(color='#C1CDD4')))
top<-p2a+heat(ce,'Cell-level Wilcoxon')+heat(pb,'Sample-level edgeR')+plot_layout(widths=c(1.15,1.8,1.8))
save('F02_expression_and_coverage',page(top/p2d+plot_layout(heights=c(1.1,1)),'Expression and genetic-context evidence','Original whole-study testing families are retained; missing comparisons are never filled with zero.'),8.5,
 'a, Integer-count bulk reconstruction (90 lesions, 81 controls); dark points meet the original full-gene significance rule. DDX58 is mapped by GeneID 23586 to the source symbol RIGI. b-c, Six broad cell types; +/minus marks FDR<0.05 and absolute log2FC>0.25, crosses are undefined or untested, and corner triangles mark clipping at +/-3. Pseudobulk units are GEO subject labels; the original study reports 14 patients. d, Support and coverage under the recorded bulk, cell, SMR/HEIDI and H4>0.8 rules. Partial coverage is separate from complete tested nonsupport. SMR/colocalization share MVP outcome statistics and are not independent confirmation.')

# 3. Newly computed enrichment and an induced graph from the archived STRING query.
enriched<-list();for(lib in c('biological_process','molecular_function','cellular_component','KEGG')){
 d<-read(paste0('enrichment/primary/',lib));d<-d[d$fdr_bh<.05,];if(nrow(d)){d$library<-if(lib=='KEGG')'KEGG' else paste0('GO ',gsub('_',' ',lib));enriched[[lib]]<-d}}
en<-do.call(rbind,enriched);panel('F03a_enrichment',en)
en$label<-vapply(sub(' - Homo sapiens \\(human\\)$','',en$description),function(x)paste(strwrap(x,34),collapse='\n'),character(1));en$label<-factor(en$label,levels=rev(en$label))
p3a<-ggplot(en,aes(-log10(fdr_bh),label))+geom_segment(aes(x=0,xend=-log10(fdr_bh),yend=label),color='#CDD9DF',linewidth=.7)+geom_point(aes(size=overlap,color=library))+
 scale_color_manual(values=c('GO biological process'=C$teal,'KEGG'=C$purple),name=NULL)+scale_size_continuous(range=c(2.5,4.5),breaks=sort(unique(en$overlap)),name='Overlap')+
 labs(x='-log10 BH FDR',y=NULL,title='Enrichment after candidate-set correction',subtitle='Background requires MVP and at least one external outcome to be estimable')+theme(panel.grid.major.y=element_blank())
rel<-read('drugs/all_drug_relations');mapping<-read('drugs/STRING_mapping');net<-read('drugs/STRING_induced_network')
targets<-sort(unique(rel$uniprot[rel$organism=='Homo sapiens']));lab<-mapping$preferredName[match(targets,mapping$queryItem)];lab[is.na(lab)]<-paste0(targets[is.na(lab)],' (unmapped)')
grid<-expand.grid(uniprot=e$uniprot,target_uniprot=targets,stringsAsFactors=FALSE);grid$gene<-e$gene[match(grid$uniprot,e$uniprot)];grid$target<-lab[match(grid$target_uniprot,targets)]
a<-mapping$stringId[match(grid$uniprot,mapping$queryItem)];b<-mapping$stringId[match(grid$target_uniprot,mapping$queryItem)]
key<-paste(pmin(net$stringId_A,net$stringId_B),pmax(net$stringId_A,net$stringId_B),sep='::');grid$score<-net$score[match(paste(pmin(a,b),pmax(a,b),sep='::'),key)]
grid$unmapped<-is.na(a)|is.na(b);grid$identity<-grid$uniprot==grid$target_uniprot;panel('F03b_network',grid)
p3b<-ggplot(grid,aes(target,gene,fill=score))+geom_tile(width=.94,height=.85,color='white')+ys()+
 geom_point(data=grid[grid$identity,],shape=23,size=2,color=C$purple,fill='white')+geom_text(data=grid[grid$unmapped,],label='\u00d7',size=2.5,color='#83939E')+
 scale_fill_gradient(low='#BFD6DE',high=C$ink,limits=c(.4,1),na.value='#F1F4F6',name='STRING score')+
 labs(x=NULL,y=NULL,title='Functional links to drug-target components',subtitle='Diamonds: same listed component; crosses: unmapped target; blank tiles: no returned edge')+
 theme(panel.grid=element_blank(),axis.text.x=element_text(angle=65,hjust=1,size=6.5))
save('F03_functional_context',page(p3a/p3b+plot_layout(heights=c(1,1.25)),'Functional context of the primary candidates','No direct single-protein drug-target record is present in the archived snapshot.'),7.7,
 'a, All terms with BH FDR<0.05 after rerunning enrichment on the primary candidate set; all eligible terms, including zero overlap, enter the appropriate GO or KEGG family. Protein-level opportunity background is 2,045, followed by unique-gene and annotation filtering. b, Exact induced graph from the archived STRING v12.5 query, not a new API response. IL12B has complex-component drug records; these and network links do not establish direct drug binding or efficacy. P11473 is unmapped and explicitly retained. Blank mapped pairs have no returned edge above 0.4.')

# 4. Clinical family preserved and a complete MVP-free role-reversed control.
clinical<-read('clinical_primary');clinical$display<-pmin(1,pmax(-1,clinical$beta));clinical$clipped<-abs(clinical$beta)>1
clinical$label<-ifelse(clinical$BH_original_family_374<.05,'\u25cf',ifelse(clinical$BH_candidate_family_exploratory<.05,'\u25cb',''))
endpoint<-read('clinical_endpoints');endlabels<-c(K11_CD_STRICT2='Crohn disease',K11_IBD_STRICT='IBD',K11_UC_STRICT2='Ulcerative colitis',E4_HYPERLIPNAS='Other hyperlipidaemia',NAFLD='NAFLD',M13_RHEUMA='Rheumatoid arthritis',I9_CHD='Major CHD events',I9_HEARTFAIL='Heart failure',C_STROKE='Stroke',T2D='Type 2 diabetes',C3_OTHER_SKIN_EXALLC='Non-melanoma skin cancer',C3_MELANOMA_SKIN_EXALLC='Melanoma',I9_CORATHER='Coronary atherosclerosis',I9_HYPTENS='Hypertension',E4_THYROIDITAUTOIM='Autoimmune thyroiditis',E4_OBESITY='Obesity',M13_SJOGREN='Sjogren syndrome',SLE_FG='SLE',M13_GOUT='Gout',E4_GRAVES_STRICT='Graves disease',L12_PSORIASIS='Psoriasis',M13_PSORIARTH='Psoriatic arthropathies')
panel('F04a_clinical',clinical)
p4a<-ggplot(clinical,aes(phenocode,gene,fill=display))+geom_tile(width=.94,height=.85,color='white')+ys()+
 geom_text(aes(label=label),size=2,color=C$ink)+geom_point(data=clinical[clinical$clipped,],shape=15,size=.5,position=position_nudge(x=.27,y=.2))+
 scale_x_discrete(limits=endpoint$phenocode,labels=endlabels)+scale_fill_gradient2(low=C$blue,mid='white',high=C$orange,limits=c(-1,1),name='Log OR (cap 1)')+
 labs(x=NULL,y=NULL,title='Exploratory clinical associations',subtitle='Filled: BH<0.05 in the earlier 374-test family. Open: only in the post hoc 154-test family.')+
 theme(panel.grid=element_blank(),axis.text.x=element_text(angle=60,hjust=1,size=6.8))
control<-read('MVP_free_control_discoveries');control<-control[control$branch=='unique_variant_givw'&control$model=='fixed_effect',];control<-control[order(control$FinnGen_p),];control$gene<-factor(control$gene,levels=rev(control$gene));panel('F04b_MVP_free',control)
p4b<-ggplot(control,aes(-log10(FinnGen_p),gene))+geom_vline(xintercept=-log10(.05/2709),linetype=2,color='#AAB8C0',linewidth=.3)+
 geom_segment(aes(x=0,xend=-log10(FinnGen_p),yend=gene),color='#D2DDE2')+geom_point(aes(shape=control_followup,color=in_MVP_primary_set),size=2.4)+
 scale_shape_manual(values=c('TRUE'=16,'FALSE'=1),labels=c('TRUE'='Supported','FALSE'='Not supported'),name='Stuart same-direction P<0.05')+scale_color_manual(values=c('TRUE'=C$teal,'FALSE'='#8B9DA8'),labels=c('TRUE'='Primary candidate','FALSE'='Control only'),name=NULL)+
 labs(x='-log10 FinnGen P',y=NULL,title='Control without MVP statistics',subtitle='FinnGen discovery; Stuart nominal same-direction support; existing datasets with roles reversed')+
 theme(panel.grid.major.y=element_blank(),legend.box='vertical')
save('F04_clinical_and_MVP_free',page(p4a/p4b+plot_layout(heights=c(1.15,1)),'Clinical context and dependence on MVP','The role-reversed control is post hoc; associated proteins may share the same genetic variant.'),9.0,
 'a, Primary candidates across the fixed 22 FinnGen R12 endpoints. Filled circles retain BH from the earlier reconstruction of 374 comparisons; open circles indicate significance only under the new 154-comparison exploratory family. Color saturates at +/-1 with corner squares marking clipping. Two phenotype proxies retain their approved definitions. b, All FinnGen Bonferroni discoveries at the unchanged 2,709 denominator, with Stuart same-direction nominal support and primary-set membership indicated. FinnGen/Stuart were already used elsewhere; this is not a new independent replication. ERP29, PPP1CC and BRAP share a frozen outcome variant, so protein counts are not counts of independent loci.')

capture.output(sessionInfo(),file=file.path(out,'figure_sessionInfo.txt'))
manifest<-list(status='final_figures_generated',inputs=setNames(as.list(vapply(unique(inputs),sha,character(1))),substring(unique(inputs),nchar(root)+2)),figures=records,panel_tables=panel_records)
jsonlite::write_json(manifest,file.path(out,'figures_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
captions<-unlist(lapply(names(records),function(n)c(paste0('## ',n),'',records[[n]]$caption,'',paste0('![',n,'](figures/',n,'.png)'),'')))
writeLines(c('# Figure captions','',captions),file.path(out,'FIGURE_CAPTIONS.md'))
stopifnot(system2('pdfunite',shQuote(c(file.path(fig,paste0(names(records),'.pdf')),file.path(out,'FIGURE_ATLAS.pdf'))))==0)
exports<-c(list.files(fig,pattern='[.](pdf|svg|png)$',full.names=TRUE),list.files(panels,full.names=TRUE),file.path(out,c('figures_manifest.json','FIGURE_ATLAS.pdf','FIGURE_CAPTIONS.md','figure_sessionInfo.txt')))
child<-list(status='final_figures_generated',inputs=manifest$inputs,outputs=setNames(as.list(vapply(exports,sha,character(1))),paste0('../',substring(exports,nchar(out)+2))))
jsonlite::write_json(child,file.path(fig,'run.json'),pretty=TRUE,auto_unbox=TRUE)
