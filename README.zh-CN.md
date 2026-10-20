# 银屑病多组学：研究重建与交互式 R 课程

[English](README.md)

本仓库包含对 [Bioinformatics 论文（DOI: 10.1093/bioinformatics/btag347）](https://doi.org/10.1093/bioinformatics/btag347) 的独立重建，以及一套 00–12 共 13 章的 R 实践课程。

## 选择入口

- **在线课程：**[中文课程网站](https://aka-debug-jie.github.io/psoriasis-multiomics/zh/)。打开即可阅读讲解、代码和运行结果。
- **学生课程包：**在下方按电脑类型和语言选择下载。
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

中文课程源文件位于 [course/notebooks/](course/notebooks/)；英文 Notebook 译文位于 [translations/en/notebooks/](translations/en/notebooks/)。网页适合先看讲解和结果，学生包用于逐格运行和修改代码。

## 下载与使用

课程包已包含 R 4.4.2、JupyterLab、所需 R 包和教学数据，无需另外安装软件。

| 电脑类型 | 中文版 | 英文版 |
|---|---|---|
| Windows x64 | [ZIP](https://github.com/aka-debug-jie/psoriasis-multiomics/releases/download/v4.0.0/pmo-zh-CN-windows-x64.zip) | [ZIP](https://github.com/aka-debug-jie/psoriasis-multiomics/releases/download/v4.0.0/pmo-en-windows-x64.zip) |
| Mac · Intel | [ZIP](https://github.com/aka-debug-jie/psoriasis-multiomics/releases/download/v4.0.0/pmo-zh-CN-macos-intel.zip) | [ZIP](https://github.com/aka-debug-jie/psoriasis-multiomics/releases/download/v4.0.0/pmo-en-macos-intel.zip) |
| Mac · Apple Silicon | [ZIP](https://github.com/aka-debug-jie/psoriasis-multiomics/releases/download/v4.0.0/pmo-zh-CN-macos-arm64.zip) | [ZIP](https://github.com/aka-debug-jie/psoriasis-multiomics/releases/download/v4.0.0/pmo-en-macos-arm64.zip) |

1. 完整解压 ZIP，放到本地可写的文件夹。Mac 可在“关于本机”中查看芯片，再选择 Intel 或 Apple Silicon 版。
2. 只阅读时打开 **View course.html**；动手练习时，Windows 打开 **Run course.bat**，Mac 打开 **Run course.command**。首次启动会自动准备内置环境，需要稍等片刻。
3. 从第 **00** 章开始，按 **Shift+Enter** 逐格运行。练习本初次打开没有输出；展示页已经带有运行结果，方便对照。
4. 用 **Ctrl+S** 或 **Cmd+S** 保存 Notebook，结束时在 Jupyter 中选择 **File → Shut Down**。计算结果保存在 `resources/work/results/`。

请保留入口旁的完整 `resources/` 文件夹。课程固定使用内置 R，无需调整电脑上已有的 R 或 Python。Mac 首次打开下载文件时，可能需要按系统提示确认打开。

[发行页面](https://github.com/aka-debug-jie/psoriasis-multiomics/releases/tag/v4.0.0)提供校验值、运行记录、组件许可和对应源码。每个平台的中英文版都已运行全部 13 章，并检查普通用户权限下的断网启动和移动目录后重启。具体测试系统见[运行验证与兼容性](docs/release-validation.md)。

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

- [research/](research/)：独立科研代码与汇总结果；原始输入需另行获取。
- [course/notebooks/](course/notebooks/)：中文课程 Notebook 源文件。
- [translations/en/notebooks/](translations/en/notebooks/)：英文 Notebook 译文。
- [course/COURSE.md](course/COURSE.md)：课程说明与数据来源记录。

## 许可

项目原创代码采用 [MIT License](LICENSE)。原创课程说明、练习和图件采用 [CC BY 4.0](LICENSE-CONTENT.md)。第三方数据、代码、截图、文章和随附软件仍遵守各自许可；项目许可不会替代这些条款。请参阅课程[数据来源清单](course/data/sources.tsv)以及各课程包随附的相关声明。

This product includes software developed by the Delft University of Technology.
