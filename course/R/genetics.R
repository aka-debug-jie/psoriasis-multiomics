# 等位基因方向与广义 IVW。第 03、11 课共用，核心矩阵运算保留在此。
allele_sign <- function(ea, oa, target_ea, target_oa) {
  a <- toupper(c(ea, oa, target_ea, target_oa))
  if (anyNA(a) || any(!a %in% c('A','C','G','T')) || a[1]==a[2] || a[3]==a[4]) return(NA_real_)
  if (paste(sort(a[1:2]),collapse='') %in% c('AT','CG')) return(NA_real_)
  if (identical(a[1:2],a[3:4])) return(1)
  if (identical(a[1:2],rev(a[3:4]))) return(-1)
  complement <- c(A='T',T='A',C='G',G='C')[a[3:4]]
  if (identical(a[1:2],unname(complement))) return(1)
  if (identical(a[1:2],unname(rev(complement)))) return(-1)
  NA_real_
}
load_ld <- function() {
  p <- read_tsv('data/genetics/within_protein_ld.tsv')
  setNames(p$r, paste(p$i,p$j,sep=':'))
}
ld_for <- function(d, ld) {
  if (nrow(d)==1) return(matrix(1,1,1))
  ix <- outer(d$LD_index,d$LD_index,paste,sep=':')
  matrix(unname(ld[ix]),nrow(d)) * outer(d$LD_sign,d$LD_sign)
}
freeze_instruments <- function(x, ld, max_abs_r=.95) {
  # 先按暴露 F 排序，同一蛋白/实际变异只保留一行，尚未读取疾病效应。
  x <- x[order(-x$F_statistic,x$cohort,x$assay,x$rsid),]
  x <- x[!duplicated(x[c('uniprot','variant_key')]),]
  out <- lapply(split(x,x$uniprot),function(d) {
    d <- d[order(-d$F_statistic,d$cohort,d$assay,d$rsid),]
    keep <- 1L; reason <- 'single_unique_variant'
    if(nrow(d)>1) {
      r <- ld_for(d,ld)
      if(any(d$LD_index<0)|any(!is.finite(r))) reason <- 'single_IV_missing_or_unresolved_LD' else {
        keep <- integer()
        for(j in seq_len(nrow(d))) if(!length(keep)||all(abs(r[j,keep])<=max_abs_r)) keep<-c(keep,j)
        z <- r[keep,keep,drop=FALSE]
        if(min(eigen(z,symmetric=TRUE,only.values=TRUE)$values)<=0 || kappa(z,exact=TRUE)>10000) {
          keep<-1L;reason<-'single_IV_unstable_LD'
        } else reason<-if(length(keep)<nrow(d))'LD_pruned' else 'all_unique_variants'
      }
    }
    d$selected <- seq_len(nrow(d)) %in% keep;d$selection_reason <- reason;d
  })
  do.call(rbind,out)
}
givw <- function(x,y,sy,r) {
  # Cov(beta_Yi,beta_Yj)=SE_Yi * SE_Yj * signed_LD_ij
  S <- outer(sy,sy)*r
  wx <- solve(S,x)
  denominator <- sum(x*wx)
  beta <- sum(wx*y)/denominator
  se <- sqrt(1/denominator)
  residual <- y-beta*x
  Q <- as.numeric(t(residual)%*%solve(S,residual))
  df <- length(x)-1L
  se_mre <- se*sqrt(if(df>0)max(1,Q/df) else 1)
  data.frame(beta=beta,se=se,p=2*pnorm(-abs(beta/se)),
    CI_log_low=beta-qnorm(.975)*se,CI_log_high=beta+qnorm(.975)*se,
    Q=Q,Q_df=df,Q_p=if(df>0)pchisq(Q,df,lower.tail=FALSE) else NA_real_,
    se_mre=se_mre,p_mre=2*pnorm(-abs(beta/se_mre)),n_variants=length(x))
}
analyse_outcome <- function(plan,y,ld) {
  stopifnot(!anyDuplicated(y$rsid))
  ans<-list();audits<-list()
  for(uid in unique(plan$uniprot)) {
    d<-plan[plan$uniprot==uid & plan$selected,]
    a<-y[match(d$rsid,y$rsid),]
    sign<-mapply(allele_sign,d$EA,d$OA,a$EA,a$OA)
    d$y<-a$beta*sign;d$sy<-a$se
    d$included<-is.finite(d$y)&is.finite(d$sy)&d$sy>0&!grepl('^ambiguous_aliases',d$selection_status)
    d$orientation<-sign;audits[[uid]]<-d
    v<-d[d$included,]
    row<-data.frame(uniprot=uid,gene=d$gene[1],frozen_variants=nrow(d),observed_variants=nrow(v))
    if(nrow(v)) row<-cbind(row,givw(v$beta,v$y,v$sy,ld_for(v,ld)))
    ans[[uid]]<-row
  }
  list(results=as.data.frame(data.table::rbindlist(ans,fill=TRUE)),harmonized=do.call(rbind,audits))
}
