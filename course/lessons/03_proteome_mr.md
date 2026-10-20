# 03｜全蛋白筛查与外部结局

## 研究问题与输入

本课在固定的 2,709 个蛋白家族中筛查与银屑病相关的蛋白，并查看Stuart 外部结局的方向。输入包括暴露数据、蛋白内 LD 矩阵，以及 MVP 和 Stuart 结局汇总统计。暴露行按蛋白和变异去重，再按 F 统计量排序冻结工具；`|r|≤0.95` 的规则用于保留工具。等位基因协调后，使用带 signed-LD 协方差的广义 IVW。

每个结果单位是一种蛋白，而不是一个 SNP：同一蛋白保留的变异先合并为一个广义 IVW 估计。`frozen_variants` 与 `observed_variants` 分别表示冻结工具数和结局中成功协调的工具数；`Q` 及 `se_mre` 记录异质性与离散度调整相关信息。

在 Windows 上打开 `teaching.Rproj`，于 RStudio Console 执行 `source('lessons/03_proteome_mr.R')`。详见 [R 脚本](03_proteome_mr.R)。

## 筛选与结果

发现条件为 MVP P<0.05/2,709；随后要求 Stuart P<0.05，作为外部名义支持。实际结果表含 2,709 个蛋白，7 个候选通过两步筛选：DDX58、IL7R、TNFAIP3、IL12B、SGSH、PRSS53、UBLCP1；七者在外部名义显著的方向均与 MVP 一致。举例，IL12B 的 MVP beta=-0.1726，P=1.27×10⁻¹⁴。MVP 的标准误由报告区间反推，源 P 值口径尚未确认，因此源表 `p_reported` 不一定等于用 beta/SE 推出的正态检验 P。

`candidates.tsv` 仅收录同时达到发现和外部名义条件的蛋白；未通过条件的 2,702 个蛋白仍保留在全量结果表，可用于核查阈值附近的估计。

查看 [跨结局效应图](../results/03_proteome_mr/discovery_and_validation.png)、[全蛋白表](../results/03_proteome_mr/all_proteins.tsv) 和 [候选表](../results/03_proteome_mr/candidates.tsv)。

## 练习与答案

1. 发现阈值分母为何是 2,709？**答案：** 这是预先固定的受检蛋白家族；按显著性缩小会改变检验家族。
2. 如何调整发现严格度？**答案：** 修改 Bonferroni 阈值中的分子或固定家族规模前，需先明确研究设计；提高 P 阈值会扩大候选数，但改变发现规则。
3. 一个蛋白没有某个结局工具时，表中 NA 表示什么？**答案：** 当前数据没有可协调的效应估计，不表示零效应；先查看冻结工具和协调表。也要区分 `frozen_variants=0` 与工具已冻结、但 `observed_variants=0` 的情况。
