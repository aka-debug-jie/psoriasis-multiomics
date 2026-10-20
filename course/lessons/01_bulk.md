# 01｜Bulk RNA-seq 差异表达

## 输入

本课使用 GSE54456 的 NCBI 重处理整数计数矩阵和样本资料。计数矩阵包含 171 个样本；3 个未进入计数矩阵的样本仍保留在样本资料中，并在 `sample_inclusion.tsv` 标记。表达值按 GeneID 对应当前 NCBI 符号，注释文件只包含本课所需的 GeneID、Symbol 和 Synonyms。GeneID 23586 的当前符号为 RIGI，`DDX58` 是其别名。

NCBI 的数据说明确定了 raw counts 单位；代码中的整数检查用于发现输入错误，单靠整数外观不能证明数据来源。

分析比较银屑病皮损与正常皮肤。差异判定阈值为 BH 校正后 P 值小于 0.05，且绝对 log2 fold change 至少为 1。

## 问题

1. 如何依据数据来源和数值检查，区分 raw counts 与 RPKM/TPM？
2. TMM 和 voom 分别解决什么问题？
3. 为什么 BH 校正的分母是所有进入检验的基因数，而不是候选基因数？
4. 这七个候选基因的方向、效应大小和统计显著性是否一致？

## 关键代码

```r
# 检查计数并移除全零基因
stopifnot(all(counts >= 0), all(counts == floor(counts)))
keep <- rowSums(counts) > 0
y <- edgeR::DGEList(counts[keep, , drop = FALSE], group = group)

# TMM、voom、线性模型和经验贝叶斯检验
y <- edgeR::calcNormFactors(y, method = "TMM")
design <- model.matrix(~ group)
v <- limma::voom(y, design)
fit <- limma::eBayes(limma::lmFit(v, design))
result <- limma::topTable(fit, coef = "grouppsoriasis", number = Inf,
                          adjust.method = "BH", sort.by = "none")
```

## 结果解读

分析保留 37,065 个非全零基因进行检验，方向定义为 psoriasis − normal。`logFC` 大于 0 表示皮损组表达较高；小于 0 表示皮损组表达较低。`adj.P.Val` 是在全部检验基因上进行 BH 校正的结果。火山图以校正后 P 值和效应大小共同标记差异基因；候选效应图展示七个候选基因的 log2 fold change 和经验贝叶斯调节后的 95% 区间。

候选基因包括 IL12B、SGSH、TNFAIP3、DDX58、IL7R、PRSS53 和 UBLCP1。DDX58 通过 GeneID 23586 的 NCBI 别名映射。样本 QC 表报告每个样本的文库大小和 TMM 因子；样本 QC 图用于查看文库大小分布，不能替代对样本身份或批次的检查。

本课使用的是同一研究数据的 NCBI 重处理计数，不是作者原始历史计数矩阵。表达差异只描述组间关联，不能单独说明因果关系或临床作用。

## 练习与参考答案

### 练习一：概念判断

若某基因在所有样本中的计数均为 0，它是否进入模型？全零过滤与低表达过滤有何区别？

**参考答案：** 不进入模型。全零过滤只移除完全没有观测计数的基因；低表达过滤还需要额外阈值。本课没有增加低表达过滤条件。

### 练习二：代码阅读

`model.matrix(~ group)` 中 `normal` 为什么是参照组？`grouppsoriasis` 的系数表示什么？

**参考答案：** 分组因子的水平顺序设为 normal、psoriasis，因此 normal 是参照组。`grouppsoriasis` 表示皮损组相对于正常组的 log2 表达变化。

### 练习三：结果解释

在 `candidate_effects.tsv` 中找出 `adj.P.Val < 0.05` 且 `abs(logFC) >= 1` 的候选基因，并说明系数方向。

**参考答案：** IL12B、DDX58、IL7R 和 PRSS53 同时满足两项阈值，logFC 为正，表示在银屑病皮损样本中表达较高。其余候选不满足两项联合判定；阈值结果不等同于机制或因果证据。

运行 `source("lessons/01_bulk.R")`，对照[完整脚本](01_bulk.R)、[火山图](../results/01_bulk/bulk_volcano.png)和[候选效应图](../results/01_bulk/candidate_effects.png)。
