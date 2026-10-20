# 04｜工具选择、方向与模型敏感性

## 研究问题与输入

本课查看遗传效应是否随工具选择、患病率假设或 MR 方法而变化。进入本课前的2709蛋白家族已按编码基因位置排除30个MHC蛋白；本课另标记工具SNP是否位于25–34 Mb区域，基因位置与工具位置是两种检查。脚本先用 3,802 条暴露记录计算 F 统计量并标记 MHC 区域；再比较主分析与更严格 LD 规则。多工具方法部分直接使用已按 P<5×10⁻⁶、r²<0.1、10 Mb clumping 得到的固定输入表，不在本课重新筛工具。另一类 PGS 分析使用 P<5×10⁻⁸、r²<0.01，不能混作同一套阈值。

在 Windows 上打开 `teaching.Rproj`，于 RStudio Console 执行 `source('lessons/04_sensitivity.R')`。详见 [R 脚本](04_sensitivity.R)。

## 可见分析与结果

主 LD 分支使用 |r|≤0.95，严格分支使用 |r|≤√0.1。两个分支各有 2,709 个蛋白结果，P<0.05 的行数均为 193；可估计 beta 的行分别为 2,096 和 2,091。输入 F 的最小值为 21.54，中位数 479.61；3,802 行中有 4 行位于脚本标记的 MHC 区域。Steiger 网格含 72 个情景组合，在 K=1%、2%、3%、5% 和两种 R² 算法下，forward 均为 TRUE。多工具方法表 88 行中 55 行完成估计，33 行因少于 3 个工具而不运行中位数、众数或 Egger。

`LD_and_dispersion.tsv` 同时保存每个蛋白的主、严格 LD 结果及异质性调整列；同样的名义 P<0.05 行数不意味着每个蛋白效应完全相同。Steiger 表中的 forward 标记只对应列出的患病率与 R² 算法。多工具比较包括 IVW、中位数、众数和 Egger，单工具行会按脚本改用 Wald 比值。

结果见 [方法比较图](../results/04_sensitivity/method_comparison.png)、[LD/离散度表](../results/04_sensitivity/LD_and_dispersion.tsv)、[Steiger 网格](../results/04_sensitivity/steiger_grid.tsv) 和 [多工具方法表](../results/04_sensitivity/multi_instrument_methods.tsv)。

## 练习与答案

1. 严格 LD 分支为何用 √0.1？**答案：** helper 接收 |r| 上限，而设计阈值是 r²=0.1，所以传入 √0.1。
2. 怎样改变 Steiger 假设？**答案：** 修改 K 情景或 R² 方法后重新比较各组合；结论只适用于列出的输入假设。不要把所有组合均为 forward 解读成对其他患病率也成立。
3. 为什么一些方法标记需至少 3 个工具？**答案：** 加权中位数、众数和 Egger 在过少工具时不能按该流程估计；表中状态记录了此限制，不能用 0 代替结果。可按 `cohort`、`gene`、`mode` 查看哪些输入行可计算。
