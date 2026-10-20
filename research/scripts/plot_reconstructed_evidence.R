script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'../..'))
out<-file.path(root,'research/results/evidence_substitute');dir.create(out,recursive=TRUE,showWarnings=FALSE)
sha<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1]
inputs<-character();read_tsv<-function(rel){p<-file.path(root,rel);inputs<<-c(inputs,p);read.delim(p,stringsAsFactors=FALSE,check.names=FALSE)}
targets<-read_tsv('research/results/mr_fenland_substitute/followup_proteins.tsv')
cells<-read_tsv('research/results/evidence_substitute/cell_Wilcoxon_followup_contrasts.tsv')
donors<-read_tsv('research/results/evidence_substitute/donor_pseudobulk_followup_contrasts.tsv')
cells<-cells[cells$gene_id %in% targets$gene_id,c('gene_id','cell_type','log2FC','fdr_bh')];cells$analysis<-'Cell-level Wilcoxon'
donors<-donors[donors$gene_id %in% targets$gene_id,c('gene_id','cell_type','logFC','FDR')];names(donors)[3:4]<-c('log2FC','fdr_bh');donors$analysis<-'Donor-level edgeR'
data<-rbind(cells,donors)
grid<-expand.grid(gene_id=targets$gene_id,cell_type=unique(cells$cell_type),analysis=c('Cell-level Wilcoxon','Donor-level edgeR'),stringsAsFactors=FALSE)
data<-merge(grid,data,by=c('gene_id','cell_type','analysis'),all.x=TRUE)
data$gene<-targets$gene[match(data$gene_id,targets$gene_id)];data$evidence<--log10(pmax(data$fdr_bh,1e-10));data$display_log2FC<-pmax(-3,pmin(3,data$log2FC))
pdf(file.path(out,'cell_and_donor_expression.pdf'),width=12,height=8)
print(ggplot2::ggplot(data,ggplot2::aes(x=cell_type,y=gene,fill=display_log2FC,size=evidence))+ggplot2::geom_point(shape=21,color='grey35',na.rm=TRUE)+ggplot2::geom_point(data=data[is.na(data$fdr_bh),],mapping=ggplot2::aes(x=cell_type,y=gene),inherit.aes=FALSE,shape=4,color='grey65',size=1.5)+ggplot2::facet_wrap(~analysis)+ggplot2::scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-3,3))+ggplot2::scale_size_continuous(range=c(.5,5),limits=c(0,10))+ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=50,hjust=1))+ggplot2::labs(x=NULL,y=NULL,fill='log2 FC (clipped)',size='-log10(FDR), cap 10',title='Expression evidence for 17 follow-up proteins',subtitle='Grey crosses: gene/contrast not tested; source gene identities retained in mapping table'));dev.off()
enrichment<-list()
for(lib in c('GO_biological_process','GO_molecular_function','GO_cellular_component','KEGG')){
  d<-read_tsv(paste0('research/results/functional_substitute/',lib,'.tsv'));d<-head(d[order(d$fdr_bh,d$p_value),],8);d$library<-lib;enrichment[[lib]]<-d
}
d<-do.call(rbind,enrichment);d$description<-vapply(d$description,function(s)paste(strwrap(s,width=48),collapse='\n'),character(1));d$rank<-factor(paste(d$library,d$term),levels=rev(paste(d$library,d$term)));labels<-setNames(d$description,paste(d$library,d$term))
pdf(file.path(out,'GO_KEGG_enrichment.pdf'),width=12,height=11)
print(ggplot2::ggplot(d,ggplot2::aes(x=-log10(fdr_bh),y=rank,fill=fdr_bh<.05))+ggplot2::geom_col()+ggplot2::facet_wrap(~library,scales='free_y',ncol=2)+ggplot2::scale_y_discrete(labels=labels)+ggplot2::scale_fill_manual(values=c('FALSE'='grey75','TRUE'='#2166AC'),name='FDR < .05')+ggplot2::theme_bw()+ggplot2::labs(x='-log10(BH FDR)',y=NULL,title='GO and KEGG enrichment against tested-protein background'));dev.off()
net<-read_tsv('research/results/drug_network_substitute/STRING_network.tsv');mapping<-read_tsv('research/results/drug_network_substitute/STRING_mapping.tsv');relations<-read_tsv('research/results/drug_network_substitute/all_drug_target_relations.tsv')
nodes<-unique(mapping[c('stringId','preferredName')]);candidate_ids<-mapping$stringId[mapping$queryItem %in% targets$uniprot];drug_ids<-mapping$stringId[mapping$queryItem %in% relations$uniprot]
nodes$role<-ifelse(nodes$stringId %in% candidate_ids,ifelse(nodes$stringId %in% drug_ids,'Follow-up and drug component','Follow-up protein'),'Drug target/component')
graph<-igraph::graph_from_data_frame(net[c('stringId_A','stringId_B','score')],directed=FALSE,vertices=nodes)
set.seed(20261017);layout<-igraph::layout_with_fr(graph,weights=igraph::E(graph)$score)
colors<-c('Follow-up protein'='#B2182B','Drug target/component'='#2166AC','Follow-up and drug component'='#762A83')
pdf(file.path(out,'STRING_functional_network.pdf'),width=12,height=10)
par(mar=c(1,1,4,1));plot(graph,layout=layout,vertex.color=colors[igraph::V(graph)$role],vertex.frame.color=NA,vertex.size=6,vertex.label=igraph::V(graph)$preferredName,vertex.label.cex=.7,vertex.label.color='black',edge.color='grey75',edge.width=1.5*igraph::E(graph)$score,main='STRING 12.5 functional network (score > 0.4)')
legend('topleft',legend=names(colors),pch=19,col=colors,bty='n',cex=.8);mtext('Network edges are functional associations; drug target/component relations are tabulated separately.',side=3,cex=.75);dev.off()
record<-list(status='figures_generated',inputs=setNames(as.list(vapply(c(script,inputs),sha,character(1))),c(script,inputs)),outputs=setNames(as.list(vapply(file.path(out,c('cell_and_donor_expression.pdf','GO_KEGG_enrichment.pdf','STRING_functional_network.pdf')),sha,character(1))),c('cell_and_donor_expression.pdf','GO_KEGG_enrichment.pdf','STRING_functional_network.pdf')))
jsonlite::write_json(record,file.path(out,'figures.json'),pretty=TRUE,auto_unbox=TRUE)
