# 数据与统计量说明

## 文件与身份

`data/table_catalog.tsv`给出每张数据表的记录数、列数及SHA256；`column_catalog.tsv`按列位置登记字段。TSV中的带引号多行文字按CSV/TSV记录解析，不以物理文本行数计数。

| 文件 | 一行代表什么 | 关键身份 |
|---|---|---|
| `primary_all_proteins.tsv` | 一个固定合格蛋白及当前结果 | UniProt，2709项 |
| `primary_candidates.tsv` / `integrated_evidence.tsv` | 一个当前主候选 | UniProt和规范Ensembl ID，7项 |
| `genetic_effects.tsv` | 一个主候选在一个疾病GWAS的估计 | UniProt＋outcome |
| `primary_sensitivity.tsv` | 一个主候选在一种方法情景下的结果 | UniProt＋branch＋scenario＋model |
| `primary_delta_grid.tsv` | 一个主候选的误差相关性情景 | UniProt＋rho_xx＋rho_xy |
| `steiger_selected_instruments.tsv` | 一个实际主工具的方向诊断情景 | instrument_id＋rsid＋prevalence＋exposure_R2_estimator |
| `multi_instrument_methods.tsv`及对应coverage/输入表 | 一来源、工具范围及方法的敏感性结果 | cohort＋UniProt＋mode＋method；不跨来源合并 |
| `SMR_HEIDI.tsv` | 一个组织中的eQTL probe结果 | gene_id＋tissue＋ProbeID |
| `coloc_mRNA.tsv` / `coloc_protein.tsv` | 一个基因、组织/来源、窗口的共定位 | gene_id＋tissue（若有）＋window_kb |
| `expression_cell.tsv` / `expression_pseudobulk.tsv` | 一个基因在一种细胞类型中的比较 | gene_id＋cell_type |
| `clinical_primary.tsv` | 一个候选与一个临床结局的比较 | UniProt＋phenocode，154项 |
| `enrichment/*/*.tsv` | 一个完整检验家族内的GO/KEGG条目 | 背景分支＋库＋term |
| `MVP_free_control_all.tsv` | 一个蛋白在一种FinnGen主导对照中的结果 | UniProt＋branch＋model；每分支2709项 |
| `drugs/*` | 药物记录、复合体成员、网络边或网络关联 | 保留分子、target、UniProt及STRING身份 |

基因符号用于展示，连接优先使用UniProt、Ensembl和GeneID。DDX58在bulk源表中为GeneID 23586、符号RIGI。表达版本桥接见`expression_gene_identity.tsv`。单细胞字段`donor_id`沿用GEO subject标签，图中使用样本单位；原研究报告14名病例，部分标签后缀没有额外人员身份说明。

## 主结果与历史对照

`beta`为每来源标准化蛋白效应单位对应的log OR；`se`为相关工具固定效应gIVW的一阶标准误。多来源共同效应尺度属于模型假设，Fenland总体重建的估计SD不等于观测SD。

MVP缺少发布的标准误时，将CI按95%对数Wald区间解释，使用`(log(upper)-log(lower))/(2×1.95996398454)`重建SE。发布端尚未确认这一口径。`P_equivalent`用源P强制构造正态双尾等效SE，只用于敏感性；不是恢复的原始SE。FinnGen和Stuart使用所获文件报告的SE。

`CI_compatible_only`仅保留现有舍入模型下P/CI相容的工具，未决或缺少审计结果的工具也不纳入。它是依赖结局统计量的诊断，不能当作独立验证；因此不可计算的候选不记成验证失败。

`primary_selected`是当前主候选标记：MVP P<0.05/2709且至少一个外部GWAS P<0.05。方向另列，没有事后增加“双外部同向”作为选择门槛。当前7项恰好在两套外部数据均同向名义显著。

所有`legacy_independent_weight_*`列属于旧独立权重模型，不能用作当前筛选标记。`author_*`属于作者报告，完整29项另表保存。`reconstructed_total_score`和UKB PGS/Cox字段保持NA。

## 检验家族

- 主发现及FinnGen主导对照均保留2709分母，包含未能计算的合格蛋白。
- 富集主背景为2045个MVP及至少一个外部结局可算蛋白，映射后2043个唯一基因；数据库再按注释覆盖形成各自背景。2096个MVP可算蛋白对应2094个唯一基因，只作背景敏感性。GO BP/MF/CC、KEGG分别BH，包含合格的零重叠条目。
- `BH_original_family_374`保留此前项目重建17×22的家族调整，作为临床图主要标记；这里的original指此前重建，并非论文29项名单。`BH_candidate_family_exploratory`是七项×22的154家族新调整，只作事后探索；MRE两列遵循同样区别。
- 表达及SMR原有全基因/全probe校正不因抽取七项而重算。`SMR_FDR_all_tested_probes`与名义SMR＋HEIDI规则分列；缺失HEIDI不当作通过。
- Steiger只提取实际主工具，保留四个患病率和两种R²估计；N为来源队列总量近似。多工具表保留按来源的cis/cis＋trans、IVW/加权中位数/众数/Egger结果与不足工具数状态，不合并为独立支持次数。

## 缺失、支持与网络

NA表示未取得、未计算或未定义，不能填0。`evidence_matrix.tsv`同时保留计划和可算单位数：有实际支持为`supported`；覆盖不全且无观察到的支持为`partial_no_observed_support`，整体`support`保持NA；完全不可评估为`not_evaluable`。`observed_support`只描述已经可检验的部分。综合表中各类`*_any_support`同样只概括原表实际受检比较，不推断未检验部分。

SMR和共定位使用MVP统计量，`MVP_statistics_used=true`；它们不能当作完全独立于MVP来源问题的证据票数。H4>0.8是所用共定位模型的判定阈值，不是位点已获因果证明的概率。

`drug_direct`是单蛋白靶点记录；`drug_components`是复合体/家族成员；`drug_network_links`是经STRING功能关联连接到药物记录。三者不合并为“直接药物作用”。无返回边是本快照中无满足阈值的边，P11473另记为未映射。现有图是对归档查询的精确诱导子图，不声称重新查询API会逐字返回相同结果。

FinnGen主导对照不使用MVP统计量，但重用了既有结局数据，是事后角色重排。ERP29、PPP1CC和BRAP均使用rs10774624；蛋白数量不等于独立遗传信号数量。
