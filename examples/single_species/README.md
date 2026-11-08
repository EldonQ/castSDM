# 单物种应用示例（castSDM >= 0.28.1）

入口是 `00_RUN_ALL.R`。案例保留数据特有的发生记录质控、分析域和未来文件匹配；
拟合、空间交叉验证、缓存、四角场景比较和归因均调用已安装的 castSDM。
默认生成 9 个联合未来的比较，并导出完整地图和诊断图；逐变量归因分析报告场景。

## 直接在 R / RStudio 中运行

直接运行，无需先重启 R：

```r
source("E:/Package/castSDM_workflows/single_species/00_RUN_ALL.R", encoding = "UTF-8")
```

不需要先 `setwd()`，也不需要手动 `library()` 或设置 `.libPaths()`。
入口以脚本所在位置定位配置，自动选择项目库中最新的 castSDM，包括
`E:/Package/.rlib_cast` 和 `E:/Package/.rlib_cast_versions/<版本>`。
显式设置的 `CAST_R_LIBRARY` 或示例内 `.library` 优先。不会在正式运行中自动联网安装软件。
若仍有缺失或过旧依赖，会一次列出并给出安装命令。已经加载的 castSDM 若与项目版
一致，会直接复用；若是旧版或开发版，入口会只切换 castSDM 到项目安装版，保留工作区
对象及其他已加载包，不再仅因路径不同要求重启。极少数其他包正在导入 castSDM、
无法切换的情况，可改用下面的 Rscript 命令独立运行。

默认运行 `Schizothoracinae`。改物种或另开一个结果目录：

```r
Sys.setenv(CAST_SPECIES = "Callosciurus_erythraeus", CAST_RUN_ID = "analysis_01")
source("E:/Package/castSDM_workflows/single_species/00_RUN_ALL.R", encoding = "UTF-8")
```

恢复默认物种和运行标识用 `Sys.unsetenv(c("CAST_SPECIES", "CAST_RUN_ID"))`。
`01_config.R` 是科学参数和资源预算的编辑入口，其他辅助文件通常不需要改动。

## 安装与路径

本机项目库已经就绪。只有更换机器、更新包或入口确实报告缺依赖时，才需要安装。
开发源码修改后须重新安装，再运行示例。已有任务运行时使用独立版本库，避免改写其加载的包：

```r
lib <- "E:/Package/.rlib_cast_versions/0.28.1"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, "E:/Package/.rlib_cast", .libPaths()))
install.packages(c("remotes", "terra", "sf", "data.table", "ranger", "maxnet",
  "mgcv", "glmnet", "pROC", "blockCV", "FNN", "digest", "filelock",
  "ggplot2", "patchwork", "scales", "ragg"),
  lib = .libPaths()[1L], repos = "https://cloud.r-project.org")
remotes::install_github("EldonQ/castSDM", lib = .libPaths()[1L])
# 或在本机安装源码：
# install.packages("E:/Package/cast", repos = NULL, type = "source", lib = .libPaths()[1L])
```

`01_config.R` 集中定义协议与案例参数。默认路径对应当前工作目录的数据布局：

| 内容 | 默认位置 / 设置 |
|---|---|
| 环境库 | `E:/CHECO26/N-SDM/data`，可用 `CAST_DATA_ROOT` 改写 |
| 发生记录 | `E:/Package/castSDM_workflows/species_occ`，可用 `CAST_OCCURRENCE_ROOT` 改写 |
| 示例根目录 | 脚本所在目录，可用 `CAST_EXAMPLE_ROOT` 改写 |
| R 包库 | 自动识别 `.rlib_cast` 和 `.rlib_cast_versions`；自定义位置可设 `CAST_R_LIBRARY` |
| 物种 | `Schizothoracinae`，可用 `CAST_SPECIES` 选择 `Callosciurus_erythraeus` |
| 运行标识 | `cast028_glu` 或 `cast028_hlu`，可用 `CAST_RUN_ID` 改写 |

GBIF 经纬度数据使用外部 PROJ `cs2cs` 转到环境库 CRS，须把其可执行文件加入
`PATH`；这是本数据适配器的要求。`nsdm_occ` 输入应已为投影米坐标。
所有环境栅格须预先对齐，适配器发现不一致会停止，不会自动重采样。

## 运行

在 PowerShell 中，先确保 `Rscript` 在 PATH（或使用本机 Rscript 的完整路径）：

```powershell
cd E:\Package\castSDM_workflows\single_species
Rscript 00_RUN_ALL.R --dry-run
Rscript 00_RUN_ALL.R
```

本机也可直接使用 Rscript 的完整路径；入口会自动设置项目库：

```powershell
& "D:\R\R-4.4.3\bin\x64\Rscript.exe" .\00_RUN_ALL.R --dry-run
& "D:\R\R-4.4.3\bin\x64\Rscript.exe" .\00_RUN_ALL.R
```

`--dry-run` 只检查依赖、配置和输入清单是否存在，不准备数据、不拟合、不投影。
完整运行前，若环境库发生变化或缺少 `config/*.csv`，显式运行一次：

```powershell
Rscript 00_INVENTORY.R
```

该清单扫描包含覆盖率计算，可能较慢，不随每次拟合重复执行。源码仓库仅跟踪
示例代码，数据、生成的 `config` 清单和 `output` 不上传。当前本地案例已有清单可用。

选择另一个案例或开始独立结果目录：

```powershell
$env:CAST_SPECIES = "Callosciurus_erythraeus"
$env:CAST_RUN_ID = "analysis_01"
Rscript 00_RUN_ALL.R
```

相同运行标识用于续跑，包会校验缓存内容和输出摘要。更改科学输入后，应设置新的
`CAST_RUN_ID`；训练缓存独立于结果目录，可在数据和训练设置一致时复用。
损坏或陈旧的比较栅格不会被当作成功结果使用。

0.28.1 修复初始化函数携带源码引用、导致工作进程提前加载系统旧包的问题。
工作进程在隔离启动后固定加载主进程的安装库，
实际版本、路径和核心代码由主进程核对，并在读取大型栅格前完成预检。
本补丁不改变预测方法或训练设置，兼容 0.28.0 的训练、预测和图件数值指纹。配置和科学脚本未变的
`cast028_*` 案例会自动续跑，版本升级另记入 `maintenance.rds`；训练结果不会被改标。

## 输出与重绘

结果位于 `output/<运行标识>/<物种>/`：

- `training/`：包生成的模型、折外预测、验证指标、阶段记录和 ODMAP 报告。
- `scenario_table.csv`：所选动态变量对应的精确未来文件及年份。
- `comparison/`：四角预测、八项对比、配对模型 SD、支持域及表格。
- `diagnostics/`：单模型与冻结集成的留出评价、训练层评分诊断、空间折组成与 AOA、选择频率、共同支持域和模型方向冲突。
- `threshold_sensitivity.csv`：阈值 0.3、0.5、0.7 的敏感性，预测时累计，之后直接读取计数。
- `occurrence_years.csv`、`case_design.csv`：实际建模记录的年代分布、分类层级、年份字段含义及设计限制。
- `figures/`：完整 PNG（600 dpi）与 PDF 图集，宽 183 mm；沿用中国省界、南海插图与原有配色。
- `figures/source_data/`：每类图对应的数值表、版本与显示设置、可用性说明。
- `figures/figure_inventory.csv`、`figure_status.csv`：全部可生成图的清单、实际导出状态及缺失原因。
- `figures/.data/`：响应曲线、驱动分解、MESS、原生网格范围分类和绘图网格的独立缓存。
- `case.rds`：案例配置和结果路径；`completion.rds` 仅表示本次计算与导出完成。

数据准备保存在 `output/.prepared/<物种>/<内容指纹>/`，训练阶段缓存保存在
`output/.training-cache/`。不要把 `completion.rds` 当作模型或论文质量验收证明。

```powershell
Rscript 04_REPLOT.R output/cast028_glu/Schizothoracinae
```

默认输出 01–10：空间验证、生境变化、环境驱动、支持域与模型分歧、两组情景分解、
所有变量的响应曲线、模型性能、范围变化类别、MESS 外推和全部未来适宜度。
扩展输出包括变量筛选及稳定性、阈值敏感性、模型权重与损失、固定响应剖面、AOA，
以及所有可用情景对比层、配对 SD、支持域、方向冲突和每个动态驱动的地图集。
变量筛选概览保留前 30 名和全部最终变量，所有有效分数另以每页最多 45 个变量的
`S_selection_*` 图集完整输出，避免超长图或缩成小字。数值计算得到的零 p 值在对数图中
使用正规化最小正数作为显示下限，原始分数保存在 CSV 中。
输入不足的图会明确记录原因，不静默省略。

首次为旧结果补齐图集时，需要从已保存模型计算缺失的响应曲线、逐变量地图和 MESS；
这一步可能较慢，完成后单独缓存，不重新拟合模型或改写原始科学栅格。
以后的默认重绘复用校验通过的缓存。严格仅重绘时加 `--render-only`；
如果缓存不完整，该模式会直接提示，不调用模型预测。
原生栅格先按固定阈值分类，绘图网格默认 0.02 度；连续图用双线性、类别图用最近邻，
统计表仍采用原生网格。图上标明显示色限与被压缩到色限的像元数。

在 R 会话中重绘已经恢复图集的旧 `cast026`，可运行：

```r
Sys.setenv(CAST_CASE_DIR = "E:/Package/castSDM_workflows/single_species/output/cast026_glu/Schizothoracinae",
           CAST_RENDER_ONLY = "1")
source("E:/Package/castSDM_workflows/single_species/04_REPLOT.R", encoding = "UTF-8")
Sys.unsetenv(c("CAST_CASE_DIR", "CAST_RENDER_ONLY"))
```

只更新部分图时可设置 `CAST_FIGURES = "02,05,10"`，完成后清除该环境变量。
重绘部分图会保留其他图的导出记录。旧 `cast026` / `cast027` 的完整结果目录也可补图；
`04_REPLOT.R` 要求科学阶段已完成，不对仍在计算的目录操作。新训练使用新目录。

## 可调参数与解释边界

`01_config.R` 的 `protocol` 是训练科学设置的唯一入口。总资源预算默认
12 核、3 个场景进程、8000 MB；比较阶段每进程最多 4 个 RF 线程，训练阶段按
协议最多用 12 个 RF 线程。修改配置开头的 `resources` 即同步控制两阶段资源。
完整图集默认计算并缓存报告场景的逐变量归因。`figure_response_grid` 和
`figure_display_resolution` 只控制响应图网格和地图显示分辨率。

GLU 的中期/末期土地利用快照为 2050/2100，气候代表年为 2055/2085；
HLU 使用与气候相同的 2030/2055/2085 年。映射会记录在场景表中，二者不能在论文
中写成完全相同的时间窗口。时间匹配只使用训练发生记录的年份；没有年份时明确跳过。
裂腹鱼案例为亚科记录集合，不能在论文中称为单一物种。

对比是固定拟合响应下的模型条件场景差值，不识别真实生态因果效应。
`*_main` 为基线简单效应，`*_average` 为平均效应；SD 是模型间差异，不是置信区间。
支持检查使用边际训练分位数。二值结果使用 `>=`；表格计数单位为有效像元，
面积解释须核实等面积投影和像元面积。新版本修复 AOA 最近邻定义，旧 AOA 数值须重算。

新版默认写入 `cast028_*`，不会覆盖已有 `cast026_*` / `cast027_*`。0.27 修正了协议截断在
CV、留出评价和地图之间的一致性，后续指标可能变化，不能给旧指标换上新版标签。
评分偏移只固定加权训练均值；`score_views_oof.csv` 比较同一拟合的原生与偏移评分，
不使用留出集选择策略，也不代表更换 RF 拟合权重后的效果。
0.27.1 同时修复保留 R 源码引用时的缓存校验误失效；旧缓存不删除，按新版指纹重算。

`summary_common_support.csv`、`binary_common_support.csv` 对所有未来使用同一支持域；
原有支持域表仍按每个情景分别计算。`agreement.csv` 报告模型正负方向冲突，不能
解释为显著性。常数预测器、Brier 和 PR-AUC 的解释都依赖发生/背景的抽样比例。
诊断提供冻结阈值遗漏率，以及训练 OOF 发生分数 5%/10% 分位阈值的留出遗漏率；
原来的 `omission_5/10` 是评价样本分数的顶部比例诊断，定义不同。

年份筛选默认关闭。确认 `year` 的含义后，可在新分析配置中设置
`occurrence_year_range = c(2000, 2022)` 和 `occurrence_keep_unknown_year = FALSE`。
该选择会进入数据准备缓存标识，并在 QC 中记录剔除原因，不自动将记录年视作采样年。
科学配置改变时，运行声明会在改写既有结果前要求新的 `CAST_RUN_ID`；诊断阈值和
硬件预算可以在同一次分析中更改。不同数据准备写入不同内容指纹目录。
绘图维护使用合成数据与已保存的真实案例结果验证，不重新拟合完整案例。
