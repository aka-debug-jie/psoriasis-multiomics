# Presto implements tie- and continuity-corrected asymptotic Wilcoxon tests.
# Fold changes use the installed Seurat 5.3 definition, independently of Presto logFC.
sc_contrast_table <- function(normalized, groups, focal, cell_type, comparison, features, wilcox=NULL) {
  stopifnot(length(groups)==ncol(normalized), all(rownames(normalized)==features$gene_id))
  chosen <- groups==focal
  n1 <- sum(chosen); n0 <- sum(!chosen)
  stopifnot(n1>=3, n0>=3)
  if (is.null(wilcox)) wilcox <- presto::wilcoxauc(normalized, groups)
  w <- wilcox[wilcox$group==focal, ]
  w <- w[match(features$gene_id, w$feature), ]
  stopifnot(!anyNA(w$feature))
  linear <- normalized
  linear@x <- expm1(linear@x)
  sum1 <- Matrix::rowSums(linear[,chosen,drop=FALSE])
  sum0 <- Matrix::rowSums(linear[,!chosen,drop=FALSE])
  p <- w$pval
  p[!is.finite(p)] <- NA_real_
  p[sum1==0 & sum0==0] <- NA_real_
  fdr <- p.adjust(p, method="BH", n=nrow(features))
  fc <- log2((sum1+1)/n1)-log2((sum0+1)/n0)
  data.frame(gene_id=features$gene_id,symbol=features$symbol,cell_type=cell_type,
             comparison=comparison,p_value=p,log2FC=fc,fdr_bh=fdr,
             n_group1=n1,n_group0=n0,mean_normalized_group1=sum1/n1,mean_normalized_group0=sum0/n0,
             mean_log_group1=Matrix::rowMeans(normalized[,chosen,drop=FALSE]),
             mean_log_group0=Matrix::rowMeans(normalized[,!chosen,drop=FALSE]),
             pct_group1=w$pct_in/100,pct_group0=w$pct_out/100,
             test_status=ifelse(is.finite(p),"tested","constant_or_undefined_rank_variance"),
             BH_family_size=nrow(features), stringsAsFactors=FALSE)
}
