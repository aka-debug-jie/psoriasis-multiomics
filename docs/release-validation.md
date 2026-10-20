# v4.0.0 validation and compatibility / 运行验证与兼容性

Both language editions were run on each of these systems with bundled R 4.4.2. These are the systems actually tested; other OS versions have not been tested individually.

中英文版都在下表系统中使用内置 R 4.4.2 运行过。表中列的是实际验收环境，其他系统版本尚未逐一测试。

| Package / 课程包 | Tested system / 验收系统 | Native run / 运行记录 |
|---|---|---|
| Windows x64 | Windows Server 2022, build 20348 | [39483044033](https://github.com/aka-debug-jie/psoriasis-multiomics/actions/runs/39483044033) |
| Mac Intel | macOS 15.7.9, x86_64 | [39483043765](https://github.com/aka-debug-jie/psoriasis-multiomics/actions/runs/39483043765) |
| Mac Apple Silicon | macOS 15.7.9, arm64 | [39483044085](https://github.com/aka-debug-jie/psoriasis-multiomics/actions/runs/39483044085) |

## What was checked / 检查内容

Each of the six editions passed:

每个版本都通过了以下检查：

| Check / 检查项 | Per edition / 每个版本 |
|---|---|
| Full chapter execution / 完整章节运行 | 13 |
| Input file identities / 输入文件校验 | 93 |
| Result checks / 结果检查 | 21 |
| Numerical comparisons / 数值比对 | 60 |
| Single-cell comparisons / 单细胞比对 | 7 |

The actual launcher started Jupyter and an R kernel under an ordinary user account with external networking blocked. The test also reopened the running course, shut it down, moved the folder, and started it again while preserving saved work. Paths included Chinese characters and spaces.

测试使用实际启动入口，在普通用户权限、外网连接被阻断的条件下启动 Jupyter 和 R 内核，并检查重复打开、关闭服务、移动文件夹后重启及保留已保存内容。测试路径包含中文和空格。

The downloaded exercise notebooks have empty outputs and execution counters. Reading pages contain executed examples. Single-cell checks compare cells, annotations and statistics; UMAP pixels are not an acceptance criterion.

下载包中的练习本没有预填输出或执行编号，展示页保留了实际运行结果。单细胞检查比较细胞集合、注释和统计结果，不要求 UMAP 的每个像素相同。

## Package identity and sources / 文件身份与源码

The [v4.0.0 release](https://github.com/aka-debug-jie/psoriasis-multiomics/releases/tag/v4.0.0) includes a manifest for each package with its byte count, SHA256, tested OS and result checks. The native course build used commit [`14da69d`](https://github.com/aka-debug-jie/psoriasis-multiomics/tree/14da69d16b71fa5ac2b835b5417a7bb2d7a9d332). Later documentation additions record the previous archive hash and verify that all existing member contents and ZIP attributes remain unchanged.

[v4.0.0 发行页面](https://github.com/aka-debug-jie/psoriasis-multiomics/releases/tag/v4.0.0)为每个包提供清单，记录文件大小、SHA256、测试系统和结果检查。原生课程构建对应提交 [`14da69d`](https://github.com/aka-debug-jie/psoriasis-multiomics/tree/14da69d16b71fa5ac2b835b5417a7bb2d7a9d332)。后续补充说明文件时，保留了原压缩包哈希，并逐项验证已有文件的内容和 ZIP 属性未变。

Component notices are under `resources/licenses/`. `SOURCE_DOWNLOAD.txt` points to the matching source archive index. Windows packages also include `WINDOWS_COMPONENT_REVIEW.json`, which explains the retained collection records, line-ending differences and supplemental toolchain sources. These source downloads are separate from the software students need to run the course.

组件许可位于 `resources/licenses/`，`SOURCE_DOWNLOAD.txt`提供对应源码索引。Windows 包另含 `WINDOWS_COMPONENT_REVIEW.json`，说明保留的采集记录、换行差异和补充工具链源码。学生运行课程不需要另行下载这些源码。
