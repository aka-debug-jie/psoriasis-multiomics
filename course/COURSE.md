# 统计函数与参考资料

<a id="appendix-genetics"></a>

## 附录A｜共享遗传统计函数的完整实现

第二、三、四、六、八、十一章的 `source("R/genetics.R")` 载入下面这组函数。这里直接展示[源码](R/genetics.R)，便于从调用追到运算。阅读顺序是：先看等位基因如何改变符号，再看同一变异怎样去重，最后看广义 IVW 的矩阵公式。正文代码已载入它们，不需要另外复制本附录才能运行。

| 函数 | 输入与作用 |
|---|---|
| `allele_sign()` | 根据两对等位基因返回1、−1或NA；回文位点按当前规则留为无法唯一协调 |
| `load_ld()` / `ld_for()` | 读取蛋白内部的相关系数，并按暴露效应等位基因转换符号 |
| `freeze_instruments()` | 先按暴露F及固定次序去重，再按LD选工具，结果不依赖疾病P |
| `givw()` | 构造疾病效应协方差S，求共同效应、SE、Q和离散度敏感性 |
| `analyse_outcome()` | 对齐某一结局、记录可用工具，再对每个蛋白调用估计器 |

`givw()` 中的 `solve(S, x)` 直接解线性方程。分母 `sum(x * solve(S, x))` 表示当前工具在协方差结构下提供的信息量。当只有一个工具时，广义 IVW 回到 Wald 比值。`se_mre` 使用 `max(1, Q/df)` 放大标准误，因此异质性不会在该规则下把区间变窄；它与主固定效应SE分列。这里的一阶SE未计入暴露估计误差，第二章的delta公式展示了单工具的另一种近似。

```r
allele_sign <- function(ea, oa, target_ea, target_oa) {
    a <- toupper(c(ea, oa, target_ea, target_oa))
    if (anyNA(a) || any(!a %in% c("A", "C", "G", "T")) || a[1] == a[2] || a[3] == a[4]) 
        return(NA_real_)
    if (paste(sort(a[1:2]), collapse = "") %in% c("AT", "CG")) 
        return(NA_real_)
    if (identical(a[1:2], a[3:4])) 
        return(1)
    if (identical(a[1:2], rev(a[3:4]))) 
        return(-1)
    complement <- c(A = "T", T = "A", C = "G", G = "C")[a[3:4]]
    if (identical(a[1:2], unname(complement))) 
        return(1)
    if (identical(a[1:2], unname(rev(complement)))) 
        return(-1)
    NA_real_
}
load_ld <- function() {
    p <- read_tsv("data/genetics/within_protein_ld.tsv")
    setNames(p$r, paste(p$i, p$j, sep = ":"))
}
ld_for <- function(d, ld) {
    if (nrow(d) == 1) 
        return(matrix(1, 1, 1))
    ix <- outer(d$LD_index, d$LD_index, paste, sep = ":")
    matrix(unname(ld[ix]), nrow(d)) * outer(d$LD_sign, d$LD_sign)
}
freeze_instruments <- function(x, ld, max_abs_r = 0.95) {
    x <- x[order(-x$F_statistic, x$cohort, x$assay, x$rsid), ]
    x <- x[!duplicated(x[c("uniprot", "variant_key")]), ]
    out <- lapply(split(x, x$uniprot), function(d) {
        d <- d[order(-d$F_statistic, d$cohort, d$assay, d$rsid), ]
        keep <- 1L
        reason <- "single_unique_variant"
        if (nrow(d) > 1) {
            r <- ld_for(d, ld)
            if (any(d$LD_index < 0) | any(!is.finite(r))) 
                reason <- "single_IV_missing_or_unresolved_LD"
            else {
                keep <- integer()
                for (j in seq_len(nrow(d))) if (!length(keep) || all(abs(r[j, keep]) <= max_abs_r)) 
                  keep <- c(keep, j)
                z <- r[keep, keep, drop = FALSE]
                if (min(eigen(z, symmetric = TRUE, only.values = TRUE)$values) <= 0 || kappa(z, 
                  exact = TRUE) > 10000) {
                  keep <- 1L
                  reason <- "single_IV_unstable_LD"
                }
                else reason <- if (length(keep) < nrow(d)) 
                  "LD_pruned"
                else "all_unique_variants"
            }
        }
        d$selected <- seq_len(nrow(d)) %in% keep
        d$selection_reason <- reason
        d
    })
    do.call(rbind, out)
}
givw <- function(x, y, sy, r) {
    S <- outer(sy, sy) * r
    wx <- solve(S, x)
    denominator <- sum(x * wx)
    beta <- sum(wx * y)/denominator
    se <- sqrt(1/denominator)
    residual <- y - beta * x
    Q <- as.numeric(t(residual) %*% solve(S, residual))
    df <- length(x) - 1L
    se_mre <- se * sqrt(if (df > 0) max(1, Q/df) else 1)
    data.frame(beta = beta, se = se, p = 2 * pnorm(-abs(beta/se)), CI_log_low = beta - qnorm(0.975) * 
        se, CI_log_high = beta + qnorm(0.975) * se, Q = Q, Q_df = df, Q_p = if (df > 0) 
        pchisq(Q, df, lower.tail = FALSE)
    else NA_real_, se_mre = se_mre, p_mre = 2 * pnorm(-abs(beta/se_mre)), n_variants = length(x))
}
analyse_outcome <- function(plan, y, ld) {
    stopifnot(!anyDuplicated(y$rsid))
    ans <- list()
    audits <- list()
    for (uid in unique(plan$uniprot)) {
        d <- plan[plan$uniprot == uid & plan$selected, ]
        a <- y[match(d$rsid, y$rsid), ]
        sign <- mapply(allele_sign, d$EA, d$OA, a$EA, a$OA)
        d$y <- a$beta * sign
        d$sy <- a$se
        d$included <- is.finite(d$y) & is.finite(d$sy) & d$sy > 0 & !grepl("^ambiguous_aliases", 
            d$selection_status)
        d$orientation <- sign
        audits[[uid]] <- d
        v <- d[d$included, ]
        row <- data.frame(uniprot = uid, gene = d$gene[1], frozen_variants = nrow(d), observed_variants = nrow(v))
        if (nrow(v)) 
            row <- cbind(row, givw(v$beta, v$y, v$sy, ld_for(v, ld)))
        ans[[uid]] <- row
    }
    list(results = as.data.frame(data.table::rbindlist(ans, fill = TRUE)), harmonized = do.call(rbind, 
        audits))
}
```

<a id="appendix-common"></a>

## 附录B｜读表、绘图与外部程序

运行入口已经选好 R 和课程包库。公共函数负责读表、保存结果和绘图，统计建模步骤保留在各课代码中。下面的实现与随包源文件一致。

```r
source('R/load_environment.R')
# 在 teaching.Rproj 打开的项目根目录中运行各课脚本。
options(stringsAsFactors = FALSE, scipen = 6)
Sys.setenv(OMP_NUM_THREADS = "2", OPENBLAS_NUM_THREADS = "2", MKL_NUM_THREADS = "2")
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(2)
read_tsv <- function(path) read.delim(if (grepl("\\.gz$", path)) gzfile(path) else path,
  check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
  quote = FALSE, na = "NA", fileEncoding = "UTF-8")
lesson_output <- function(id) {
  path <- file.path(Sys.getenv("TEACHING_OUTPUT", "results"), id)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
theme_course <- function() ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(colour = "#566573"),
    axis.line = ggplot2::element_line(linewidth = .35),
    legend.position = "bottom", plot.margin = ggplot2::margin(12, 16, 12, 12))
save_plot <- function(p, output, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(output, paste0(name, ".pdf")), p, width = width, height = height)
  ggplot2::ggsave(file.path(output, paste0(name, ".png")), p, width = width, height = height,
                 dpi = 320, bg = "white")
  if(isTRUE(getOption("teaching.inline",FALSE))) {
    options(repr.plot.width=width,repr.plot.height=height)
    print(p)
  }
}
source('R/tool_support.R')
run_tool <- function(program, args, log) {
  code <- system2(program, vapply(args, shQuote, ""), stdout = log, stderr = log)
  if (code != 0) stop("外部程序运行失败；请查看 ", log)
  invisible(code)
}

```

`tool_path()` 按系统和处理器选择随包的 PLINK 或 SMR；程序定位代码见 [R/tool_support.R](R/tool_support.R)。结果写入“课程资源/work/results”，日志与结果放在对应课次目录。

<a id="appendix-methods"></a>

## 附录C｜第三方方法、数据引用与后续阅读

第四章调用的是固定版本的 TwoSampleMR 方法源码，保存在 `R/vendor/`。这些函数的调用参数已在正文完整展示；原始实现、版权和版本说明见[方法来源](R/vendor/README.md)。

| 调用 | 本章用途 | 阅读时关注 |
|---|---|---|
| `mr_ivw()` | 多工具IVW估计 | 权重、异质性与标准误定义 |
| `mr_weighted_median()` | 加权中位数估计 | 有效工具权重比例及bootstrap |
| `mr_weighted_mode()` | 加权众数估计 | 比值分布、带宽与最大一致组 |
| `mr_egger_regression()` | 含截距的MR回归 | 截距、工具方向和额外识别假设 |
| `get_r_from_lor()` | 二分类效应换算为相关系数 | 人群患病率、病例/对照数及频率 |
| `coloc::coloc.abf()` | 区域五假设的后验比较 | 区域完整性、效应方差与先验 |
| `survival::coxph()` / `cox.zph()` | Cox估计与PH检验 | 时间轴、事件编码、协变量和风险集 |

包函数的完整参数可以在R中用帮助页面查询，例如 `?limma::voom`、`?coloc::coloc.abf`、`?survival::coxph`。方法引用包括 Hemani 等人的 [MR-Base 论文](https://doi.org/10.7554/eLife.34408)，以及 [SMR 官方说明](https://yanglab.westlake.edu.cn/software/smr/)、[PLINK 文档](https://www.cog-genomics.org/plink/1.9/)。所有数据的具体版本和访问入口保存在[来源清单](data/sources.tsv)，其使用条件沿用各提供方要求。

从[课程起步页](notebooks/00_start.nb.html)选择章节，按每课的代码顺序运行并完成练习。
