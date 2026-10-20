check_course_results <- function(output='results') {
  read<-function(lesson,file)read.delim(file.path(output,lesson,file),check.names=FALSE)
  tests<-list();test<-function(name,condition){tests[[name]]<<-data.frame(check=name,passed=isTRUE(condition));if(!isTRUE(condition))stop('检查失败：',name)}
  source('R/genetics.R')
  test('allele_swap',allele_sign('A','G','G','A')==-1)
  test('strand_complement',allele_sign('A','G','T','C')==1)
  test('palindromic_unresolved',is.na(allele_sign('A','T','A','T')))
  w<-givw(.4,.08,.02,matrix(1,1,1));test('Wald_analytic',abs(w$beta-.2)<1e-12 && abs(w$se-.05)<1e-12)
  b<-read('01_bulk','all_gene_differential_expression.tsv')
  test('bulk_testing_family',nrow(b)==37065 && sum(b$significant)==4266)
  a<-read('03_proteome_mr','all_proteins.tsv');candidates<-read('03_proteome_mr','candidates.tsv')
  test('full_proteome_family',nrow(a)==2709)
  test('candidate_set',setequal(candidates$gene,c('IL12B','SGSH','TNFAIP3','DDX58','IL7R','PRSS53','UBLCP1')))
  h<-read('03_proteome_mr','harmonization.tsv');h<-h[h$included,]
  test('one_GWAS_record_per_variant',!anyDuplicated(h[c('outcome','uniprot','variant_key')]))
  smr<-read('05_smr','two_probe_results.tsv')
  test('SMR_valid_statistics',all(is.finite(smr$p_SMR_teaching)) && all(is.finite(smr$p_HEIDI_teaching)) &&
       all(smr$p_SMR_teaching>=0 & smr$p_SMR_teaching<=1) && all(smr$p_HEIDI_teaching>=0 & smr$p_HEIDI_teaching<=1))
  test('SMR_example_decisions',identical(smr$nominal_SMR_and_HEIDI_rule,c(FALSE,FALSE)) && all(smr$topSNP_absolute_difference==0))
  cl<-read('06_colocalization','coloc_summary.tsv');post<-as.matrix(cl[paste0('PP.H',0:4,'.abf')])
  test('coloc_probability_sum',max(abs(rowSums(post)-1))<1e-10 && all(post>=0&post<=1))
  test('coloc_region_not_significant_only',all(cl$n_variants>1000))
  sc<-read('07_single_cell','cell_qc.tsv')
  test('sample_units',length(unique(sc$gsm))==22 && nrow(sc)==10708 && sum(sc$keep)==7938)
  pg<-read('08_pgs','individual_pgs.tsv');test('R_PLINK_score',nrow(pg)==503 && max(abs(pg$PGS-pg$PLINK_sum))<1e-5)
  pc<-read('08_pgs','genetic_pcs.tsv');test('five_genetic_PCs',all(paste0('PC',1:5)%in%names(pc)) && nrow(pc)==503)
  # 同时变更 dosage 与 beta 时原始评分平移，标准化评分不变。
  d<-read.table(file.path(output,'08_pgs/genotypes.raw'),header=TRUE,check.names=FALSE)
  weights<-read('08_pgs','weights.tsv');dose<-as.matrix(d[,paste(weights$SNP,weights$A1,sep='_')]);raw<-dose%*%weights$weight
  flipped<-(2-dose)%*%(-weights$weight)
  test('score_allele_coding',max(abs(flipped-raw+2*sum(weights$weight)))<1e-10 && max(abs(scale(flipped)-scale(raw)))<1e-10)
  lung<-read.delim('data/cox/lung_ncctg228.tsv');test('Cox_event_definition',sum(lung$status==2)==165 && sum(complete.cases(lung[c('time','age','sex','ph.ecog')]))==227)
  cox<-read('09_cox','cox_model.tsv');test('HR_scale',max(abs(cox$HR-exp(cox$beta)))<1e-10)
  enrichment<-read('10_functional','enrichment_all_families.tsv')
  test('GO_Reactome_families',identical(as.integer(table(factor(enrichment$family,levels=c('BP','MF','CC','Reactome')))),c(1888L,338L,282L,342L)))
  test('zero_overlap_terms_retained',all(enrichment$p_value[enrichment$overlap==0]==1))
  f<-read('11_clinical_mr','clinical_all_374.tsv');z<-read('11_clinical_mr','clinical_candidates_154.tsv')
  test('clinical_full_multiplicity',nrow(f)==374 && nrow(z)==154 && sum(z$BH_original_family_374<.05)==21 && sum(z$BH_candidate_family_exploratory<.05)==22)
  # 第12课已经对作者全部29分与缺失项进行了逐项断言。
  result<-do.call(rbind,tests);write.table(result,file.path(output,'checks.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  print(result,row.names=FALSE);invisible(result)
}
if(sys.nframe()==0L)check_course_results(Sys.getenv('TEACHING_OUTPUT','results'))
