source('R/common.R')
library(Seurat)
library(Matrix)
out<-lesson_output('07_single_cell')
set.seed(20261017)
# 每个样本在 QC 前最多取 500 个细胞，所有基因保留；正常/皮损各 8/14 个样本。
con<-gzfile('data/single_cell/counts.mtx.gz');counts<-as(Matrix::readMM(con),'CsparseMatrix')
features<-read_tsv('data/single_cell/features.tsv');cells<-read_tsv('data/single_cell/cells.tsv')
rownames(counts)<-features$gene_id;colnames(counts)<-cells$cell_id;rownames(cells)<-cells$cell_id
stopifnot(nrow(counts)==nrow(features),ncol(counts)==nrow(cells),length(unique(cells$gsm))==22L)
qc<-cells
qc$n_counts<-Matrix::colSums(counts);qc$n_features<-Matrix::colSums(counts>0)
mt<-toupper(as.character(features$mitochondrial))=='TRUE'
qc$percent_mito<-100*Matrix::colSums(counts[mt,,drop=FALSE])/qc$n_counts
qc$keep<-with(qc,n_features>=500 & n_features<=6000 & n_counts>=1000 & n_counts<=30000 & percent_mito<=10)
qc$keep[is.na(qc$keep)]<-FALSE
write_tsv(qc,file.path(out,'cell_qc.tsv'))
write_tsv(aggregate(keep~gsm+donor_id+condition,qc,sum),file.path(out,'sample_qc.tsv'))
obj<-CreateSeuratObject(counts[,qc$keep],meta.data=cells[qc$keep,])
obj<-NormalizeData(obj,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
obj<-FindVariableFeatures(obj,selection.method='vst',nfeatures=2000,verbose=FALSE)
# 仅对 2000 个高变基因 scale，避免生成全基因的稠密矩阵。
obj<-ScaleData(obj,features=VariableFeatures(obj),verbose=FALSE)
obj<-RunPCA(obj,features=VariableFeatures(obj),npcs=30,seed.use=20261017,verbose=FALSE)
obj<-FindNeighbors(obj,dims=1:30,k.param=20,verbose=FALSE)
obj<-FindClusters(obj,resolution=.8,algorithm=1,random.seed=20261017,verbose=FALSE)
obj<-RunUMAP(obj,dims=1:30,n.neighbors=30,min.dist=.3,metric='cosine',seed.use=20261017,verbose=FALSE,n.threads=2)
# 新的聚类必须依据标记重读，不能移用全量分析的 cluster 编号。
panels<-list(Keratinocytes=c('KRT14','KRT5','KRT10','KRT1','SBSN','KLK7'),
 Fibroblasts=c('DCN','LUM'),Endothelial=c('CLDN5','VWF'),Immune=c('HLA-DPB1','CD3D','CD3E'),
 Melanocytes=c('MLANA','DCT'),Pericytes=c('ACTA2','TAGLN'))
lognorm<-LayerData(obj,layer='data');cluster<-as.character(Idents(obj));cl<-sort(unique(cluster))
marker_rows<-list();scores<-matrix(0,length(cl),length(panels),dimnames=list(cl,names(panels)))
for(type in names(panels)) {
 ids<-features$gene_id[features$symbol %in% panels[[type]]];ids<-intersect(ids,rownames(lognorm))
 means<-sapply(cl,function(k)Matrix::rowMeans(lognorm[ids,cluster==k,drop=FALSE]))
 # 逐标记在聚类间标准化后平均，使高表达的单一基因不独占注释。
 z<-t(scale(t(means)));z[!is.finite(z)]<-0;scores[,type]<-colMeans(z)
 for(g in ids)marker_rows[[paste(type,g)]]<-data.frame(cell_type=type,gene=features$symbol[match(g,features$gene_id)],cluster=cl,
   average=as.numeric(means[g,]),percent=sapply(cl,function(k)100*mean(lognorm[g,cluster==k]>0)))
}
annotation<-data.frame(cluster=cl,suggested_type=colnames(scores)[max.col(scores,ties.method='first')],
  marker_score=apply(scores,1,max))
annotation$cell_type<-ifelse(annotation$marker_score>0,annotation$suggested_type,'Unassigned')
write_tsv(annotation,file.path(out,'marker_annotation.tsv'))
write_tsv(do.call(rbind,marker_rows),file.path(out,'marker_expression.tsv'))
embedding<-as.data.frame(Embeddings(obj,'umap'));names(embedding)<-c('UMAP1','UMAP2')
embedding$cell_id<-colnames(obj);embedding$cluster<-cluster
embedding$condition<-obj$condition;embedding$cell_type<-annotation$cell_type[match(cluster,annotation$cluster)]
write_tsv(embedding,file.path(out,'umap.tsv'))
p<-ggplot2::ggplot(embedding,ggplot2::aes(UMAP1,UMAP2,colour=cell_type))+
 ggplot2::geom_point(size=.3,alpha=.7)+ggplot2::facet_wrap(~condition)+
 ggplot2::scale_colour_manual(values=c(Keratinocytes='#D9954C',Fibroblasts='#628EA0',Endothelial='#7772A5',Immune='#378F78',Melanocytes='#B66D82',Pericytes='#8A9C50',Unassigned='#AAAAAA'))+
 ggplot2::labs(title='A compact atlas across 22 skin samples',subtitle='Labels are marker-based suggestions for this sampled dataset',colour=NULL)+theme_course()+
 ggplot2::theme(axis.ticks=ggplot2::element_blank(),axis.text=ggplot2::element_blank())
save_plot(p,out,'cell_atlas',10,5.2)
m<-do.call(rbind,marker_rows)
p<-ggplot2::ggplot(m,ggplot2::aes(gene,cluster,size=percent,colour=average))+
 ggplot2::geom_point()+ggplot2::scale_colour_viridis_c()+ggplot2::labs(x=NULL,y='Cluster',size='% expressed',colour='Mean log-expression',title='Marker evidence for cluster annotation')+theme_course()+
 ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
save_plot(p,out,'marker_dotplot',11,7)
# 原方法的细胞级 Wilcoxon：在每类细胞中比较疾病组，仅作该抽样的描述性计算。
# 检验预先固定的 7 个候选，BH 家族是“7候选 × 可评价细胞类型”，不等同于全基因家族。
obj$teaching_type<-embedding$cell_type
selected_genes<-c('IL12B','SGSH','TNFAIP3','DDX58','IL7R','PRSS53','UBLCP1')
gene_map<-features[features$symbol %in% c(selected_genes,'RIGI'),];gene_map$symbol[gene_map$symbol=='RIGI']<-'DDX58'
cell_de<-list();counter<-0L
for(type in sort(unique(obj$teaching_type))) for(j in seq_len(nrow(gene_map))) {
 gid<-gene_map$gene_id[j];sel<-obj$teaching_type==type
 a<-as.numeric(lognorm[gid,sel & obj$condition=='psoriasis']);b<-as.numeric(lognorm[gid,sel & obj$condition=='normal'])
 counter<-counter+1L
 cell_de[[counter]]<-data.frame(gene=gene_map$symbol[j],cell_type=type,n_psoriasis=length(a),n_normal=length(b),
   p=if(length(a)>2 && length(b)>2 && length(unique(c(a,b)))>1)wilcox.test(a,b,exact=FALSE)$p.value else NA_real_)
}
cell_de<-do.call(rbind,cell_de);cell_de$comparison<-'psoriasis_vs_normal';cell_de$BH<-p.adjust(cell_de$p,'BH');write_tsv(cell_de,file.path(out,'candidate_cell_wilcoxon.tsv'))
# 第二类比较：同一个候选在某一细胞类型与所有其余细胞中的表达。
rest_de<-list();counter<-0L
for(type in sort(unique(obj$teaching_type))) for(j in seq_len(nrow(gene_map))) {
 gid<-gene_map$gene_id[j];sel<-obj$teaching_type==type
 a<-as.numeric(lognorm[gid,sel]);b<-as.numeric(lognorm[gid,!sel]);counter<-counter+1L
 rest_de[[counter]]<-data.frame(gene=gene_map$symbol[j],cell_type=type,n_type=length(a),n_rest=length(b),
  p=if(length(a)>2 && length(b)>2 && length(unique(c(a,b)))>1)wilcox.test(a,b,exact=FALSE)$p.value else NA_real_)
}
rest_de<-do.call(rbind,rest_de);rest_de$comparison<-'type_vs_rest';rest_de$BH<-p.adjust(rest_de$p,'BH')
write_tsv(rest_de,file.path(out,'candidate_type_vs_rest_wilcoxon.tsv'))

# 供体层面：使用全量数据预聚合的 raw counts，每个样本/细胞类型只有一列。
# 这里保留每个生物样本，而不是把成千上万细胞当成独立重复。
pb<-read_tsv('data/single_cell/pseudobulk_counts.tsv.gz');pm<-read_tsv('data/single_cell/pseudobulk_samples.tsv')
pf<-read_tsv('data/single_cell/pseudobulk_features.tsv');pc<-as.matrix(pb[,-1]);rownames(pc)<-pb$gene_id
stopifnot(identical(colnames(pc),pm$sample_id))
pb_results<-list()
for(type in sort(setdiff(unique(pm$cell_type),'Unassigned'))) {
 use<-pm$cell_type==type & pm$n_cells>=10;meta<-pm[use,]
 group<-factor(meta$condition,levels=c('normal','psoriasis'));if(any(table(group)<3))next
 design<-model.matrix(~group);y<-edgeR::DGEList(pc[,use,drop=FALSE])
 keep<-edgeR::filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=FALSE]
 y<-edgeR::calcNormFactors(y,method='TMM');y<-edgeR::estimateDisp(y,design,robust=FALSE)
 fit<-edgeR::glmQLFit(y,design,robust=FALSE);fit<-edgeR::glmQLFTest(fit,coef='grouppsoriasis')
 tab<-edgeR::topTags(fit,n=Inf,sort.by='none')$table;tab$gene_id<-rownames(tab)
 tab$symbol<-pf$symbol[match(tab$gene_id,pf$gene_id)];tab$cell_type<-type;pb_results[[type]]<-tab
}
write_tsv(do.call(rbind,pb_results),file.path(out,'donor_pseudobulk.tsv'))
write_tsv(data.frame(preQC=ncol(counts),postQC=sum(qc$keep),genes=nrow(counts),clusters=nrow(annotation)),file.path(out,'summary.tsv'))
print(data.frame(preQC=ncol(counts),postQC=sum(qc$keep),clusters=nrow(annotation)))
