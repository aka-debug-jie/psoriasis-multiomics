# 银屑病多组学：研究重建与交互式 R 课程

[English](README.md)

本仓库包含对 [Bioinformatics 论文（DOI: 10.1093/bioinformatics/btag347）](https://doi.org/10.1093/bioinformatics/btag347) 的独立重建，以及一套 00–12 共 13 章的 R 实践课程。

## 选择入口

- **在线课程：**[课程网站](https://aka-debug-jie.github.io/psoriasis-multiomics/)。网站发布后再使用此链接。
- **学生课程包：**[Releases 页面](https://github.com/aka-debug-jie/psoriasis-multiomics/releases)。v4 课程包正在准备并进行原生系统验收，暂未提供下载。
- **科研源码：**从 [research/execute.R](research/execute.R) 开始，分析材料位于 [research/](research/)。这是独立的科研重建，与课程教学示例分开。

## 课程一览

课程共 13 章，编号 00–12，覆盖从项目启动到综合证据的完整方法流程。

| 章节 | 内容 |
|---|---|
| 00 | Jupyter 与 R 入门 |
| 01 | Bulk RNA-seq |
| 02 | 单靶点孟德尔随机化 |
| 03 | 蛋白组发现与验证 |
| 04 | 敏感性分析 |
| 05 | SMR 与 HEIDI |
| 06 | 区域共定位与先验敏感性 |
| 07 | 单细胞分析与候选基因 |
| 08 | IL7R 蛋白遗传评分 |
| 09 | Cox 生存分析 |
| 10 | 功能与药物关联 |
| 11 | 相关疾病 MR |
| 12 | 综合证据 |

中文课程源文件位于 [course/notebooks/](course/notebooks/)；英文 Notebook 译文位于 [translations/en/notebooks/](translations/en/notebooks/)。课程网站和学生包发布后，会分别从网站与 Releases 页面提供入口。

## 首次运行

v4 学生包完成后，请按操作系统和处理器选择压缩包，并解压到较短的目录路径。计划提供 Windows、Intel Mac、Apple silicon Mac 三种系统版本，并分别提供中文和英文包。每个包都将自带 R 4.4.2 和 Jupyter，无需另行安装 R 或 Jupyter。

Windows 双击 **Run course.bat**，Mac 双击 **Run course.command**。只想查看课程时，打开 **View course.html**；要运行和编辑课程，请在 Jupyter 中打开 Notebook，并从第 00 章开始。

课程包仍在进行原生系统验收；完成后会在 Releases 页面列出下载。

## 科研重建与范围

[research/](research/) 下的科研代码独立重建了论文中的部分分析，不是对论文全部输入、代码和结果的精确复现。部分原始科研输入无法取得或已由公开数据替代。本次重建通过分析得到 7 个候选，与论文报告的 29 个候选分别对照。

公开数据的几项选择限定了分析范围：

- 公开教学版使用 MVP 进行发现、Stuart 进行外部验证；保留的科研版本还包含 FinnGen。
- 教学示例使用 eQTL Catalogue 和 GTEx 血液数据（n=670）展示 ERAP2 与 LIPA；它们是教学示例，不属于发现候选。
- 公开教学版使用 Reactome 第 92 版替代 KEGG 通路快照。
- 临床分析中的 22 项 MVP 结局面板与 FinnGen 使用的表型定义不同。
- 未复现 UK Biobank 个体层面的 PGS 与 Cox 分析。

课程将 PGS 和 Cox 示例与上述科研分析分开。PGS 课程使用真实的 1000 Genomes Phase 3 基因型，以 IL7R 蛋白效应权重计算两个位点的评分；Cox 课程使用 R 自带数据 survival::lung。两者是互相独立的教学数据，不会合并，也不代表 UK Biobank 个体分析结果。

## 文件结构

- [research/](research/)：独立科研代码、输入、记录与结果。
- [course/notebooks/](course/notebooks/)：中文课程 Notebook 源文件和渲染页面。
- [translations/en/notebooks/](translations/en/notebooks/)：英文 Notebook 译文。
- [course/COURSE.md](course/COURSE.md)：课程说明与数据来源记录。

## 许可

项目原创代码采用 [MIT License](LICENSE)。原创课程说明、练习和图件采用 [CC BY 4.0](LICENSE-CONTENT.md)。第三方数据、代码、截图、文章和随附软件仍遵守各自许可；项目许可不会替代这些条款。请参阅课程[数据来源清单](course/data/sources.tsv)以及各课程包随附的相关声明。
