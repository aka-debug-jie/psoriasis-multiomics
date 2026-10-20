# 银屑病多组学科研复现结果

本目录统一保存公开替代数据范围内的修正版结果。主分析、敏感性分析和下游解释采用同一候选身份与明确的检验家族；历史重建保留在原结果目录。

主分析得到7个候选：IL12B、SGSH、TNFAIP3、DDX58、IL7R、PRSS53、UBLCP1。异质性、严格LD和P等效SE情景分别得到不同的6项、6项和8项集合。MVP原始发布只提供OR、CI与P，其精确统计口径仍未确认；这里交付的是保留该来源限制的科研重建，不是原输入和全部数值完全一致的复现。

- [科研报告](REPORT.md)：主要结果、方法选择、论文对照及未决来源问题。
- [四页图集](FIGURE_ATLAS.pdf)与[图注](FIGURE_CAPTIONS.md)：主遗传结果、表达和覆盖、功能背景、临床及MVP无关对照。
- [七项综合证据表](data/integrated_evidence.tsv)：当前MR、表达、SMR、共定位、药物和临床证据。
- [全部受检蛋白](data/primary_all_proteins.tsv)：2709个固定合格蛋白及计算、选择状态。
- [作者29项候选对照](data/author_all_candidates.tsv)：保留作者报告字段，与本次成员身份分列。
- [数据字典](DATA_DICTIONARY.md)与[表清单](data/table_catalog.tsv)：单位、检验家族、缺失语义及文件哈希。
- [来源与版本](SOURCES.tsv)：数据库/队列身份及来源证据；其中证据文件路径以项目根目录为基准。
- [验证记录](qa/validation.json)：独立统计核算、字段与来源核查、新目录重建比较。

在项目根目录运行统一整合入口：

```bash
research/envs/r442_sc/bin/Rscript --vanilla research/execute.R --module final_research
research/envs/r442_sc/bin/Rscript --vanilla research/execute.R --module final_research --resume
bash research/scripts/check_reproduction.sh
```

主定义位于[primary_analysis.json](../config/primary_analysis.json)，数据整合和出图分别由[Python脚本](../scripts/finalize_research.py)与[R脚本](../scripts/plot_final_research.R)完成。它们使用已运行并核验的科研模块输出；不会在整合时重新拟合bulk、单细胞或SMR。若修改上游数据或方法，须先重跑对应模块，再整合。

需要独立检查重建时，在一个尚不存在的新目录运行：

```bash
research/envs/data/bin/python research/scripts/finalize_research.py --output research/records/my_fresh_check
research/envs/data/bin/python research/tests/validate_final_research.py --repeat research/records/my_fresh_check
```

`run.json`记录分析输入及数据表SHA256；`figures_manifest.json`记录绘图输入、逐面板数据和图件SHA256；最终`manifest.json`登记整个交付目录。教学材料尚未制作。
