source('R/common.R')
source('R/genetics.R')
out <- lesson_output('02_single_mr')
x <- read_tsv('data/genetics/exposures.tsv')
y <- read_tsv('data/genetics/stuart.tsv')
# IL7R 的 UKB-PPP cis 工具；选择完全基于暴露证据。
x <- x[x$gene=='IL7R' & x$cohort=='UKB-PPP',]
x <- x[which.max(x$F_statistic),]
y <- y[match(x$rsid,y$rsid),]
sign <- allele_sign(x$EA,x$OA,y$EA,y$OA)
stopifnot(nrow(x)==1,is.finite(sign),y$se>0)
beta_y <- sign*y$beta
beta_mr <- beta_y/x$beta
se_mr <- y$se/abs(x$beta) # 一阶 NOME 标准误
se_delta <- sqrt(y$se^2/x$beta^2 + beta_y^2*x$se_exposure^2/x$beta^4)
CI <- beta_mr + c(-1,1)*qnorm(.975)*se_mr
result <- data.frame(gene='IL7R',rsid=x$rsid,exposure_EA=x$EA,outcome_EA=y$EA,
  orientation=sign,beta_exposure=x$beta,beta_outcome_aligned=beta_y,
  beta=beta_mr,se=se_mr,se_delta=se_delta,OR=exp(beta_mr),
  OR_low=exp(CI[1]),OR_high=exp(CI[2]),p=2*pnorm(-abs(beta_mr/se_mr)),F=(x$beta/x$se_exposure)^2)
write_tsv(result,file.path(out,'wald_ratio.tsv'))
# 一个等位基因编码练习：双方同时换到另一等位基因，ratio 不变。
stopifnot(isTRUE(all.equal((-beta_y)/(-x$beta),beta_mr)))
p <- ggplot2::ggplot(result,ggplot2::aes(x=OR,y=gene))+
  ggplot2::geom_vline(xintercept=1,linetype=2,colour='#8A969D')+
  ggplot2::geom_segment(ggplot2::aes(x=OR_low,xend=OR_high,yend=gene),linewidth=.9,colour='#197D8D')+
  ggplot2::geom_point(size=3,colour='#197D8D')+ggplot2::scale_x_log10()+
  ggplot2::labs(x='Odds ratio per unit genetically predicted protein',y=NULL,
    title='IL7R: one cis instrument',subtitle='UKB-PPP exposure; Stuart psoriasis outcome; 95% confidence interval')+theme_course()
save_plot(p,out,'wald_ratio',7,3.5)
print(result)
