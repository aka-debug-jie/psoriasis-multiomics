source('R/common.R');source('R/genetics.R')
out<-lesson_output('08_pgs');plink<-tool_path('plink')
x<-read_tsv('data/pgs/IL7R_source_cis_loci.tsv')
bim<-read.table('data/pgs/EUR_teaching.bim',col.names=c('chr','rsid','cm','pos','A1','A2'))
# 权重来自 IL7R 蛋白的同一 UKB-PPP 暴露研究；不使用疾病 beta。
# 数据起点是来源已标注的 cis 条件关联位点。去掉非 SNP / 方向不明回文位点。
x$valid_alleles<-mapply(function(a,b)is.finite(allele_sign(a,b,a,b)),x$EA,x$OA)
x$reference_available<-x$variant %in% bim$rsid
x$eligible<-x$neg_log10_p > -log10(5e-8) & x$valid_alleles & x$reference_available
write_tsv(x,file.path(out,'source_locus_selection.tsv'))
w<-x[x$eligible,c('variant','EA','beta')];names(w)<-c('SNP','A1','weight')
stopifnot(nrow(w)==2L)
write_tsv(w,file.path(out,'weights.tsv'))
writeLines(w$SNP,file.path(out,'score_variants.txt'))
# 验证两个评分位点在本批 EUR 基因型中 r² < .01。
run_tool(plink,c('--bfile','data/pgs/EUR_teaching','--extract',file.path(out,'score_variants.txt'),
 '--recode','A','--recode-allele',file.path(out,'weights.tsv'),'--keep-allele-order','--out',file.path(out,'genotypes'),
 '--threads','2','--memory','1500'),file.path(out,'genotypes_console.log'))
g<-read.table(file.path(out,'genotypes.raw'),header=TRUE,check.names=FALSE)
dosage<-as.matrix(g[,paste(w$SNP,w$A1,sep='_'),drop=FALSE])
complete<-complete.cases(dosage)
r2<-cor(dosage[,1],dosage[,2],use='complete.obs')^2
stopifnot(r2<.01)
write_tsv(data.frame(variants=ncol(dosage),people=nrow(dosage),complete=sum(complete),r2=r2),file.path(out,'input_summary.tsv'))
# 每人评分是效应等位基因剂量（0/1/2）乘以蛋白效应，再求和。
pgs<-as.numeric(dosage[complete,,drop=FALSE] %*% w$weight)
score<-data.frame(FID=g$FID[complete],IID=g$IID[complete],PGS=pgs,PGS_z=as.numeric(scale(pgs)))
write.table(score[,c('FID','IID')],file.path(out,'complete_people.txt'),row.names=FALSE,col.names=FALSE,quote=FALSE)
# 两个位点均完整的样本已选出，因此 PLINK 在这些评分位点不发生均值填补。
run_tool(plink,c('--bfile','data/pgs/EUR_teaching','--keep',file.path(out,'complete_people.txt'),
 '--score',file.path(out,'weights.tsv'),'1','2','3','header','sum',
 '--out',file.path(out,'plink_score'),'--threads','2','--memory','1500'),file.path(out,'score_console.log'))
check<-read.table(file.path(out,'plink_score.profile'),header=TRUE)
score$PLINK_sum<-check$SCORESUM[match(score$IID,check$IID)]
# PLINK 文本仅输出有限位数；这里对照其打印精度，而非把文本当完整双精度。
stopifnot(max(abs(score$PGS-score$PLINK_sum))<1e-5)
write_tsv(score,file.path(out,'individual_pgs.tsv'))
# 前五个遗传 PC 来自另外的 5000 个全基因组 LD 修剪 SNP。
run_tool(plink,c('--bfile','data/pgs/EUR_teaching','--extract','data/pgs/pc_variants.txt',
 '--pca','5','header','tabs','--out',file.path(out,'ancestry'),'--threads','2','--memory','1500'),file.path(out,'pca_console.log'))
pc<-read.table(file.path(out,'ancestry.eigenvec'),header=TRUE);meta<-read_tsv('data/pgs/sample_populations.tsv')
pc$population<-meta$pop[match(pc$IID,meta$sample)];write_tsv(pc,file.path(out,'genetic_pcs.tsv'))
p<-ggplot2::ggplot(score,ggplot2::aes(PGS_z))+ggplot2::geom_histogram(bins=25,fill='#197D8D',colour='white')+
 ggplot2::labs(x='Standardized IL7R protein score',y='Individuals',title='A target-specific protein score',subtitle='1000 Genomes Phase 3 EUR; two cis variants; actual genotype dosages')+theme_course()
save_plot(p,out,'protein_score',7,4.5)
p<-ggplot2::ggplot(pc,ggplot2::aes(PC1,PC2,colour=population))+ggplot2::geom_point(size=1.4,alpha=.8)+
 ggplot2::labs(title='Genetic structure from 5000 genome-wide variants',colour='Population')+theme_course()
save_plot(p,out,'ancestry_pcs',7,5)
print(data.frame(people=nrow(score),r2=r2,max_R_PLINK_difference=max(abs(score$PGS-score$PLINK_sum))))
