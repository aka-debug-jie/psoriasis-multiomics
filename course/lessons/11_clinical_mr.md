# 11｜相关疾病结局的遗传关联

公开教学面板采用 GWAS Catalog 中 22 个具有区域索引的 MVP 欧洲祖源表型，与科研版 FinnGen 面板分开。心肌梗死、基底细胞癌和甲状腺功能减退等是新增教学示例，宽泛表型按其实际名称展示。输入表列出编号、定义、样本量和差异。这些结局共享队列，部分还与 MVP 发现表型重叠，不能视为独立重复验证。

## 输入与问题

本课以 `family_17.tsv` 中的 17 个蛋白为固定分析家族，使用预先整理的暴露数据、蛋白内 LD 和 22 个 MVP EUR 结局。七个正式候选来自 `primary_candidates.tsv`。结局行、定义、病例数和对照数来自 `endpoints.tsv`。

需要回答：冻结工具后，17 个蛋白在 22 个结局中的效应方向和区间如何？七个正式候选在完整 374 项检验中的校正结果怎样？如果只看 154 个候选检验，探索性 BH 会如何变化？

## 方法

先用 `freeze_instruments()` 按蛋白内 LD 规则冻结工具，再对每个结局调用 `analyse_outcome()`。该函数返回等位基因协调结果和广义 IVW 估计。每个蛋白与结局组合保留一行；未观测到效应时保留 NA。

17 × 22 形成预先确定的 374 项检验。先对完整 374 行的 P 值计算 BH，再筛出七个正式候选形成 154 行，并保留各自原有的 374 项 BH；另计算 154 项候选集的探索性 BH。候选子集的 374 项 BH 来自全家族排序，不能用 154 个 P 值、再把 `n` 设为 374 来代替。

`endpoints.tsv` 保留 MVP 的实际 PheCode 定义。例如局部及系统性狼疮包含两类病例，不能简写为 SLE；表型名称和病例选择共同决定结果的含义。

## 关键代码

```r
plan <- freeze_instruments(exposures, ld, max_abs_r = .95)
fit <- analyse_outcome(plan, y, ld)

# 先对完整 374 项校正
all_results$BH_original_family_374 <-
  p.adjust(all_results$p, method = "BH", n = 374)

# 再取 7 候选，并单独保留探索性 154 项校正
candidate_results <- all_results[all_results$uniprot %in% formal_ids, ]
candidate_results$BH_candidate_family_exploratory <-
  p.adjust(candidate_results$p, method = "BH", n = 154)
```

## 结果解读

`clinical_all_374.tsv` 保存完整分析家族，共 374 个蛋白-结局组合；其中 308 项有可估计 P 值，66 项因没有可协调工具而保留 NA。完整家族中 37 项的 BH 小于 0.05。`clinical_candidates_154.tsv` 保存七个正式候选的 154 项组合，均有可估计结果；其中 21 项在完整 374 项 BH 下小于 0.05，22 项在 154 项候选探索性 BH 下小于 0.05。主热图颜色表示 log odds ratio，实心圆标记 374 项 BH 达到 0.05，空心圆标记只在候选探索性 BH 下达到 0.05。点标记区分校正家族，不表示独立复现。

读取每行时同时查看 `phenocode`、`phenotype`、效应方向、置信区间和 BH 列。不同疾病表型的效应不是同一临床终点；统计关联本身也不证明因果或临床效用。

## 练习与参考答案

### 练习一：多重检验

为什么不能先筛出 154 行，再通过 `p.adjust(p, "BH", n = 374)` 得到七个候选的主分析校正值？

**参考答案：** BH 使用全家族中全部可计算 P 值的排序，并固定 n=374；66 个 NA 保持缺失。仅保留候选会丢失其他可计算检验的排序信息，事后把 n 设为 374 不能恢复这些信息。应先在完整 374 行计算 BH，再保留七候选对应的结果。

### 练习二：表型核对

为什么不能把“局部及系统性狼疮”直接标为 SLE？

**参考答案：** 这个表型包含局部狼疮，病例定义比 SLE 更宽。保留实际来源名称和 PheCode，不能把名称相近当作表型等价。

### 练习三：结果阅读

选择一个正式候选，比较它在 `clinical_candidates_154.tsv` 中的 `BH_original_family_374` 与 `BH_candidate_family_exploratory`。若两者结论不同，应如何表述？

**参考答案：** 以完整 374 项家族的 BH 作为主分析结果，并把 154 项候选集 BH 明确标为探索性敏感性结果；不把候选集较小的校正分母描述成新的确认性证据。
