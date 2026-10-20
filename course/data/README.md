# 公开课程数据

本目录保存 v4 课程实际使用的输入。文件身份见 `manifest.tsv`，来源见 `sources.tsv`；各提供方的许可和引用要求单列，不由课程代码的 MIT 许可覆盖。

| 目录 | 数据与处理起点 |
|---|---|
| bulk | GSE54456 的 NCBI 重处理整数 counts，GRCh38.p13；171 个样本（90 皮损、81 正常）。不是提交者的 RPKM，也不是原论文的历史 counts。 |
| genetics | 2,709 个固定合格蛋白的公开暴露记录；MVP 用于发现，Stuart 用于外部验证。MVP SE 由报告 CI 换算，源 P 单独保留。公开课程不附带 FinnGen 原始统计子集。 |
| smr | eQTL Catalogue 重处理 GTEx 血液 QTD000356 的 ERAP2、LIPA 示例，共 670 人，各变异 an/2 为 636–670；ESD 与固定 BESD、MVP 区域统计及 1000G EUR LD。两基因按暴露 eQTL 强度选作教学示例，不是新发现靶点。 |
| coloc | Fenland Olink 的 IL7R 蛋白和 GTEx 的 ERAP2 mRNA 分别与 MVP 比较；保留可获得、共享且可协调的区域 SNP。1 Mb 请求可能超出源 cis 覆盖，实际边界见 regions.tsv。 |
| single_cell | GSE173706 的 22 个样本，QC 前每样本最多随机取 500 个细胞，seed=20261017；10,708 个细胞。pseudobulk 输入由全量 46,225 个质控细胞按样本及细胞类型聚合；不随课堂重新聚类而改变。 |
| pgs | 1000 Genomes Phase3 的 503 个 EUR 参考样本，IL7R 两个评分 SNP 与另外 5,000 个 PCA SNP。权重来自单一 UKB-PPP 蛋白研究；参考样本没有本课所需的疾病随访。 |
| cox | survival::lung 的 228 个 NCCTG 肺癌样本。time 为天，status 1 为删失、2 为死亡；缺失值保留。它与 PGS 课的样本不是同一队列。 |
| functional | 七个科研候选和 2,043 个 GeneID 的机会背景；GO 与 Reactome 第 92 版的成员表。Reactome 替代公开版未分发的 KEGG 快照。STRING、ChEMBL 的关系与直接药物靶点分别记录。 |
| clinical | 22 个 MVP 欧洲祖源 PheCode 结局的固定工具位点统计。每项编号、表型原名、病例和对照数、与旧面板的差异均在 endpoints.tsv。该教学面板与科研版 FinnGen 面板不是逐项等价的队列替换。 |
| evidence | 原作者 29 个靶点评分标记与本项目七个候选的科研结果快照，分别保存。它们用于学习证据整合；不会因运行前面的小规模公开示例而自动更新。UKB 缺失证据保持 NA。 |

SMR 示例坐标由 GRCh38 唯一映射至 GRCh37；负链变异互补等位基因，重复或不能唯一映射的位点排除。共定位进一步排除不能确定方向的回文或不兼容等位。原始效应、映射规则和来源摘要保存在对应目录的 JSON 文件中。

MVP 临床结局可能共享病例、对照及发现分析中的 GWAS。矩阵中的多格支持不能计为多次独立验证。基因层面的科研结果快照与公开教学输入的变化分别记录；数据许可不明的原始 FinnGen、eQTLGen 统计及 KEGG 成员快照不随公开包分发。
