source('R/common.R');source('R/genetics.R')
# 固定版本的公开 TwoSampleMR 方法原码；许可及出处见 R/vendor/README.md。
source('R/vendor/mr.R');source('R/vendor/mr_mode.R');source('R/vendor/add_rsq.R')
out<-lesson_output('04_sensitivity')
x<-read_tsv('data/genetics/exposures.tsv');ld<-load_ld()
x$F_recalculated<-(x$beta/x$se_exposure)^2
x$in_MHC_region<-x$chr==6 & x$pos>=25000000 & x$pos<=34000000
write_tsv(x[,c('cohort','uniprot','gene','rsid','build','F_recalculated','in_MHC_region')],file.path(out,'instrument_strength_MHC.tsv'))
# 比较较宽松和严格的 LD 选择；疾病效应不参与工具排序。
branches<-list(main=.95,strict=sqrt(.1));comparisons<-list()
for(b in names(branches)) {
 plan<-freeze_instruments(x,ld,branches[[b]])
 fit<-analyse_outcome(plan,read_tsv('data/genetics/mvp.tsv'),ld)$results
 fit$branch<-b;comparisons[[b]]<-fit
}
write_tsv(do.call(rbind,comparisons),file.path(out,'LD_and_dispersion.tsv'))
# Steiger 的人群患病率 K 与 GWAS 抽样病例比例不同；比较四个明确情景。
s<-read_tsv('data/genetics/steiger_effects.tsv')
N<-c(Iceland=35892,'UKB-PPP'=46218,Fenland_Olink_substitute_485=485)
s$N<-unname(N[s$cohort]);direction<-list();j<-0L
for(K in c(.01,.02,.03,.05)) for(method in c('standardized_effect','bse_N')) {
 rY<-get_r_from_lor(s$beta_outcome,s$eaf_outcome,15635,429252,K,model='logit',correction=FALSE)
 r2X<-if(method=='standardized_effect')2*s$eaf_exposure*(1-s$eaf_exposure)*s$beta_exposure^2 else
  (s$beta_exposure/s$se_exposure)^2/((s$beta_exposure/s$se_exposure)^2+s$N-2)
 j<-j+1L;direction[[j]]<-data.frame(instrument_id=s$instrument_id,gene=s$gene,K=K,method=method,
  R2_exposure=r2X,R2_outcome=rY^2,forward=ifelse(is.finite(r2X)&r2X>=0&r2X<1,r2X>rY^2,NA))
}
write_tsv(do.call(rbind,direction),file.path(out,'steiger_grid.tsv'))
# 从已经 LD clump 且协调的真实工具表开始，分别计算 cis / cis+trans。
# 学习重点：加权中位数、众数和 Egger 的不同假设与工具数要求。
d<-read_tsv('data/genetics/multi_instrument_harmonized_instruments.tsv')
seeds<-read_tsv('data/genetics/method_seeds.tsv');parameters<-default_parameters();parameters$nboot<-1000L
rows<-list();j<-0L
for(i in seq_len(nrow(seeds))) {
 a<-seeds[i,];v<-d[d$cohort==a$cohort & d$uniprot==a$uniprot & d$mode==a$mode,]
 if(!nrow(v))next
 for(method in c('IVW','weighted_median','weighted_mode','MR_Egger')) {
  j<-j+1L;r<-data.frame(cohort=a$cohort,uniprot=a$uniprot,gene=v$gene[1],mode=a$mode,method=method,nsnp=nrow(v),seed=a$seed)
  if(nrow(v)<3 && method!='IVW') {r$status<-'requires_at_least_three_IVs';rows[[j]]<-r;next}
  set.seed(a$seed)
  if(nrow(v)==1) {
   fit<-list(b=v$beta_outcome_aligned/v$beta,se=v$se/abs(v$beta));fit$pval<-2*pnorm(-abs(fit$b/fit$se));r$method<-'single_IV_Wald'
  } else {
   fun<-switch(method,IVW=mr_ivw,weighted_median=mr_weighted_median,weighted_mode=mr_weighted_mode,MR_Egger=mr_egger_regression)
   fit<-fun(v$beta,v$beta_outcome_aligned,v$se_for_analysis,v$se,parameters=parameters)
  }
  for(k in intersect(c('b','se','pval','Q','Q_df','Q_pval','b_i','se_i','pval_i'),names(fit)))r[[k]]<-fit[[k]]
  r$status<-if(is.finite(fit$b))'computed' else 'unavailable';rows[[j]]<-r
 }
}
result<-as.data.frame(data.table::rbindlist(rows,fill=TRUE));write_tsv(result,file.path(out,'multi_instrument_methods.tsv'))
p<-ggplot2::ggplot(result[result$status=='computed' & result$nsnp>=3,],ggplot2::aes(b,gene,colour=method))+
 ggplot2::geom_vline(xintercept=0,linetype=2,colour='grey65')+
 ggplot2::geom_pointrange(ggplot2::aes(xmin=b-1.96*se,xmax=b+1.96*se),position=ggplot2::position_dodge(.6))+
 ggplot2::facet_grid(cohort~mode,scales='free_y')+ggplot2::labs(x='MR log odds ratio (95% CI)',y=NULL,colour='Method',title='What changes across MR assumptions?')+theme_course()
save_plot(p,out,'method_comparison',11,9)
print(table(result$status))
