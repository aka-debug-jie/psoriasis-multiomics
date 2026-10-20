source("R/common.R")

out <- lesson_output("05_smr")
data_dir <- "data/smr"
smr <- tool_path("smr")
besd_prefix <- file.path(out, "eqtl_two_probes")
besd_log <- file.path(out, "make_besd.log")
smr_log <- file.path(out, "smr_heidi.log")

# 固定的 SMR/HEIDI 阈值；使用 2 个线程。
parameters <- list(
  peqtl_smr = 5e-8,
  heidi_skip_threshold = 1,
  peqtl_heidi = 0.0015654,
  heidi_ld_upper = 0.9,
  heidi_ld_lower = 0.05,
  cis_window_kb = 2000,
  frequency_difference = 0.2,
  frequency_exclusion_fraction = 0.05,
  threads = 2L
)

required <- c("eqtl_flist.txt", "mvp_regional.ma", "EUR_two_regions.bed",
              "EUR_two_regions.bim", "EUR_two_regions.fam", "reference_results.tsv")
stopifnot(all(file.exists(file.path(data_dir, required))))

# 1. 使用 SMR 官方 --make-besd 将每个 probe 的 ESD 转成 BESD。
run_tool(smr, c("--make-besd", "--eqtl-flist", file.path(data_dir, "eqtl_flist.txt"),
                "--out", besd_prefix), besd_log)
stopifnot(all(file.exists(paste0(besd_prefix, c(".besd", ".epi", ".esi")))))

# 2. 主分析使用从公开 ESD 输入转换的两 probe BESD。
# ESD 转换会按位置重新排序；HEIDI 对并列 LD 剪枝顺序敏感。
args <- c(
  "--bfile", file.path(data_dir, "EUR_two_regions"),
  "--gwas-summary", file.path(data_dir, "mvp_regional.ma"),
  "--beqtl-summary", file.path(data_dir,"public_eqtl"),
  "--peqtl-smr", parameters$peqtl_smr,
  "--heidi-skip-thresh", parameters$heidi_skip_threshold,
  "--peqtl-heidi", parameters$peqtl_heidi,
  "--ld-upper-limit", parameters$heidi_ld_upper,
  "--ld-lower-limit", parameters$heidi_ld_lower,
  "--cis-wind", parameters$cis_window_kb,
  "--diff-freq", parameters$frequency_difference,
  "--diff-freq-prop", parameters$frequency_exclusion_fraction,
  "--thread-num", parameters$threads,
  "--extract-probe", file.path(data_dir, "probes.txt"),
  "--out", file.path(out, "two_probe")
)
run_tool(smr, args, smr_log)
raw <- read_tsv(file.path(out, "two_probe.smr"))
reference <- read_tsv(file.path(data_dir, "reference_results.tsv"))
stopifnot(all(c("probeID", "p_SMR", "p_HEIDI", "b_SMR", "se_SMR", "topSNP") %in% names(raw)),
          setequal(raw$probeID, c("ENSG00000164308", "ENSG00000107798")),
          setequal(reference$probeID, raw$probeID))

# 3. 对照固定输入的 Linux SMR 1.4.3 参考运行；HEIDI 缺失保留为 NA。
columns <- c("topSNP", "b_SMR", "se_SMR", "p_SMR", "p_HEIDI")
comparison <- merge(raw[c("probeID", columns)], reference[c("probeID", columns)],
                    by = "probeID", suffixes = c("_teaching", "_reference"), sort = TRUE)
for (column in columns) {
  left <- comparison[[paste0(column, "_teaching")]]
  right <- comparison[[paste0(column, "_reference")]]
  comparison[[paste0(column, "_absolute_difference")]] <- if (is.numeric(left)) abs(left - right) else as.numeric(left != right)
}
write_tsv(comparison, file.path(out, "reference_comparison.tsv"))

# NA HEIDI is not evidence of passing the heterogeneity check.
comparison$heidi_status <- ifelse(is.na(comparison$p_HEIDI_teaching), "not_computed", "computed")
comparison$nominal_SMR_and_HEIDI_rule <- ifelse(
  is.na(comparison$p_HEIDI_teaching), NA,
  comparison$p_SMR_teaching < 0.05 & comparison$p_HEIDI_teaching > 0.01
)
write_tsv(comparison, file.path(out, "two_probe_results.tsv"))
print(comparison[c("probeID", "topSNP_teaching", "b_SMR_teaching", "p_SMR_teaching",
                   "p_HEIDI_teaching", "heidi_status", "nominal_SMR_and_HEIDI_rule")], row.names = FALSE)

# 所有平台都检查数据有效性和本例判定；版本间数值另列。
stopifnot(all(is.finite(raw$p_SMR)),all(raw$p_SMR>=0 & raw$p_SMR<=1),
 all(is.finite(raw$p_HEIDI)),all(raw$p_HEIDI>=0 & raw$p_HEIDI<=1),
 all(comparison$topSNP_absolute_difference==0))
fields<-c('b_SMR','se_SMR','p_SMR','p_HEIDI')
within_tolerance<-all(vapply(fields,function(k) {
 a<-comparison[[paste0(k,'_teaching')]];b<-comparison[[paste0(k,'_reference')]]
 all(abs(a-b)<=1e-10+1e-8*abs(b))
},logical(1)))
version_check<-data.frame(platform=.Platform$OS.type,valid_P_and_decisions=TRUE,
 numeric_matches_course_reference=within_tolerance,
 status=if(within_tolerance)'matched_reference' else 'computed_with_reference_difference')
write_tsv(version_check,file.path(out,'version_comparison.tsv'))
stopifnot(within_tolerance)
plotdata<-data.frame(gene=c('ERAP2','LIPA'),
 SMR=c(raw$p_SMR[match('ENSG00000164308',raw$probeID)],raw$p_SMR[match('ENSG00000107798',raw$probeID)]),
 HEIDI=c(raw$p_HEIDI[match('ENSG00000164308',raw$probeID)],raw$p_HEIDI[match('ENSG00000107798',raw$probeID)]))
p<-ggplot2::ggplot(plotdata,ggplot2::aes(-log10(SMR),-log10(HEIDI),label=gene))+
 ggplot2::geom_vline(xintercept=-log10(.05),linetype=2,colour='grey60')+
 ggplot2::geom_hline(yintercept=-log10(.01),linetype=2,colour='grey60')+
 ggplot2::geom_point(size=3,colour='#197D8D')+ggplot2::geom_text(nudge_x=.15,hjust=0)+
 ggplot2::expand_limits(x=4)+ggplot2::labs(x='-log10(SMR P)',y='-log10(HEIDI P)',title='Association and heterogeneity answer different questions')+theme_course()
save_plot(p,out,'SMR_HEIDI',7,4.5)
