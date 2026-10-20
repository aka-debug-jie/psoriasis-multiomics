suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(presto))
source("research/scripts/sc_de_helpers.R")
counts <- Matrix::Matrix(rbind(zero=rep(0,12),tied=c(0,0,1,1,2,3,0,0,0,1,1,2),separated=c(5,4,6,3,4,6,0,0,1,0,2,1)),sparse=TRUE)
x<-counts;x@x<-log1p(x@x)
groups<-rep(c("case","control"),each=6)
f<-data.frame(gene_id=rownames(x),symbol=rownames(x))
d<-sc_contrast_table(x,groups,"case","fixture","psoriasis_vs_normal",f)
stopifnot(is.na(d$p_value[1]),is.na(d$fdr_bh[1]))
for(i in 2:3){
 expected<-wilcox.test(as.numeric(x[i,1:6]),as.numeric(x[i,7:12]),exact=FALSE,correct=TRUE)$p.value
 stopifnot(isTRUE(all.equal(d$p_value[i],expected,tolerance=1e-10)))
}
fc<-log2((rowSums(as.matrix(counts[,1:6]))+1)/6)-log2((rowSums(as.matrix(counts[,7:12]))+1)/6)
stopifnot(isTRUE(all.equal(d$log2FC,unname(fc),tolerance=1e-10)))
stopifnot(isTRUE(all.equal(d$fdr_bh,p.adjust(d$p_value,"BH",n=3),tolerance=1e-10)))
flipped<-sc_contrast_table(x,groups,"control","fixture","reverse",f)
stopifnot(isTRUE(all.equal(d$log2FC,-flipped$log2FC,tolerance=1e-10)))
cat("PASS: sparse Wilcoxon agrees with independent R tests; ties/zero variance/BH/direction and Seurat fold-change definition checked\n")
