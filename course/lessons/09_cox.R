# 第 09 课：Cox 比例风险模型
# 在教学项目根目录运行：Rscript lessons/09_cox.R

source("R/common.R")

required <- c("survival", "ggplot2", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("请先安装课程依赖：", paste(missing, collapse = ", "))

out <- lesson_output("09_cox")
data_file <- "data/cox/lung_ncctg228.tsv"
if (!file.exists(data_file)) stop("找不到课程数据：", data_file)
dat <- read_tsv(data_file)

# 1. 先读懂结局：time 是随访天数；status=2 表示死亡，status=1 表示删失。
# Surv() 的事件变量用 TRUE/1 表示事件发生。不要直接把原始 status 当事件变量。
stopifnot(nrow(dat) == 228L, all(dat$status %in% c(1, 2)),
          all(dat$time > 0), all(dat$sex %in% c(1, 2)))
dat$event <- as.integer(dat$status == 2)
stopifnot(sum(dat$event) == 165L)

# 2. 描述完整病例筛选。模型只使用年龄、性别和医生评定 ECOG。
model_vars <- c("time", "event", "age", "sex", "ph.ecog")
missing_by_variable <- data.frame(
  variable = c("time", "age", "sex", "ph.ecog"),
  missing_n = vapply(dat[c("time", "age", "sex", "ph.ecog")], function(x) sum(is.na(x)), integer(1)),
  stringsAsFactors = FALSE
)
cc <- complete.cases(dat[model_vars])
exclusion <- data.frame(
  group = c("原始样本", "纳入主模型（完整病例）", "因模型变量缺失而排除"),
  n = c(nrow(dat), sum(cc), sum(!cc)),
  stringsAsFactors = FALSE
)
write_tsv(missing_by_variable, file.path(out, "missing_by_variable.tsv"))
write_tsv(exclusion, file.path(out, "sample_flow.tsv"))

analysis <- dat[cc, , drop = FALSE]
stopifnot(nrow(analysis) == 227L, sum(analysis$event) == 164L)
analysis$sex_f <- factor(analysis$sex, levels = c(1, 2), labels = c("男", "女"))

# 3. 构造 Surv(time, event)，拟合多变量 Cox 模型。
# age 每增加 1 岁；ph.ecog 每增加 1 级；女性相对男性。
surv_obj <- survival::Surv(time = analysis$time, event = analysis$event)
fit <- survival::coxph(
  surv_obj ~ age + sex_f + ph.ecog,
  data = analysis,
  ties = "efron",
  x = TRUE
)

# 4. 整理 HR、95% CI 和 Wald P 值；HR = exp(coef)。
coef_table <- summary(fit)$coefficients
ci <- exp(stats::confint(fit))
model_results <- data.frame(
  term = rownames(coef_table),
  beta = unname(stats::coef(fit)),
  HR = exp(unname(stats::coef(fit))),
  CI95_low = ci[, 1],
  CI95_high = ci[, 2],
  p_value = coef_table[, "Pr(>|z|)"],
  row.names = NULL,
  check.names = FALSE
)
stopifnot(all(is.finite(model_results$HR)),
          isTRUE(all.equal(model_results$HR, exp(model_results$beta))),
          all(model_results$CI95_low <= model_results$HR),
          all(model_results$HR <= model_results$CI95_high))
write_tsv(model_results, file.path(out, "cox_model.tsv"))

# 5. 画按 ECOG 分组的 Kaplan-Meier 曲线（仅用于描述；模型中的 ECOG 是连续等级）。
# 0、1、2+ 便于展示；3/4 级合并到 2+，并非主模型中重新分组。
analysis$ecog_group <- factor(ifelse(analysis$ph.ecog >= 2, "2+", as.character(analysis$ph.ecog)),
                              levels = c("0", "1", "2+"))
km <- survival::survfit(survival::Surv(time, event) ~ ecog_group, data = analysis)
km_summary <- summary(km)
km_data <- data.frame(time = km_summary$time, surv = km_summary$surv,
                      strata = sub("^ecog_group=", "ECOG ", km_summary$strata))
km_data <- rbind(data.frame(time=0,surv=1,strata=unique(km_data$strata)),km_data)
km_plot <- ggplot2::ggplot(km_data, ggplot2::aes(x = time, y = surv, colour = strata)) +
  ggplot2::geom_step(linewidth = 0.9) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2),
                              labels = scales::label_percent(accuracy = 1)) +
  ggplot2::scale_colour_manual(values = c("ECOG 0" = "#2878B5", "ECOG 1" = "#E07A2D", "ECOG 2+" = "#48976A")) +
  ggplot2::labs(title = "Survival by baseline ECOG", subtitle = "NCCTG lung cancer example; complete cases; death is the event",
                x = "Follow-up since enrollment (days)", y = "Survival probability", colour = "Baseline ECOG") +
  theme_course()
save_plot(km_plot, out, "km_by_ecog", width = 7.2, height = 5.2)

# 6. 用 Schoenfeld 残差检验比例风险假设。P 值较小提示可能偏离 PH 假设。
ph <- survival::cox.zph(fit, transform = "km")
stopifnot(all(c("age", "sex_f", "ph.ecog", "GLOBAL") %in% rownames(ph$table)))
ph_table <- as.data.frame(ph$table)
ph_table$term <- rownames(ph$table)
rownames(ph_table) <- NULL
write_tsv(ph_table, file.path(out, "ph_test.tsv"))

cat("数据来源：survival::lung；survival ", as.character(utils::packageVersion("survival")), "\n", sep = "")
cat("原始样本：", nrow(dat), "；死亡事件：", sum(dat$event), "；完整病例：", sum(cc),
    "；完整病例死亡事件：", sum(analysis$event), "\n", sep = "")
print(model_results, row.names = FALSE, digits = 4)
cat("\nPH 检验（cox.zph）：\n")
print(ph)
cat("\n结果目录：", normalizePath(out), "\n", sep = "")

model_results$label<-c('Age (per year)','Female vs male','ECOG (per grade)')
p<-ggplot2::ggplot(model_results,ggplot2::aes(HR,label))+
 ggplot2::geom_vline(xintercept=1,linetype=2,colour='grey65')+
 ggplot2::geom_pointrange(ggplot2::aes(xmin=CI95_low,xmax=CI95_high),colour='#197D8D')+
 ggplot2::scale_x_log10()+ggplot2::labs(x='Hazard ratio (95% CI)',y=NULL,title='Cox model: adjusted associations')+theme_course()
save_plot(p,out,'hazard_ratios',7,3.6)
