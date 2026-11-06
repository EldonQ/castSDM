# 单物种应用示例（castSDM >= 0.26.0）

入口是 `00_RUN_ALL.R`。案例保留数据特有的发生记录质控、分析域和未来文件匹配；
拟合、空间交叉验证、缓存、四角场景比较和归因均调用已安装的 castSDM。
默认生成 9 个联合未来的比较，逐变量归因默认关闭。

## 安装与路径

先安装包和示例需要的依赖。开发源码修改后须重新安装，再运行示例。

```r
install.packages(c("remotes", "terra", "sf", "data.table", "ranger", "maxnet",
  "mgcv", "glmnet", "pROC", "blockCV", "FNN", "digest", "filelock"))
remotes::install_github("EldonQ/castSDM")
# 或在本机安装源码：
# install.packages("E:/Package/cast", repos = NULL, type = "source")
```

`01_config.R` 集中定义协议与案例参数。默认路径对应当前工作目录的数据布局：

| 内容 | 默认位置 / 设置 |
|---|---|
| 环境库 | `E:/CHECO26/N-SDM/data`，可用 `CAST_DATA_ROOT` 改写 |
| 发生记录 | `E:/Package/castSDM_workflows/species_occ`，可用 `CAST_OCCURRENCE_ROOT` 改写 |
| 示例根目录 | 脚本所在目录，可用 `CAST_EXAMPLE_ROOT` 改写 |
| 物种 | `Schizothoracinae`，可用 `CAST_SPECIES` 选择 `Callosciurus_erythraeus` |
| 运行标识 | `cast026_glu` 或 `cast026_hlu`，可用 `CAST_RUN_ID` 改写 |

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

本次验证安装在 `E:/Package/.rlib_cast`。本机可直接使用以下命令：

```powershell
$env:R_LIBS_USER = "E:/Package/.rlib_cast"
& "D:\R\R-4.4.3\bin\Rscript.exe" .\00_RUN_ALL.R --dry-run
& "D:\R\R-4.4.3\bin\Rscript.exe" .\00_RUN_ALL.R
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

## 输出与重绘

结果位于 `output/<运行标识>/<物种>/`：

- `training/`：包生成的模型、折外预测、验证指标、阶段记录和 ODMAP 报告。
- `scenario_table.csv`：所选动态变量对应的精确未来文件及年份。
- `comparison/`：四角预测、八项对比、配对模型 SD、支持域及表格。
- `threshold_sensitivity.csv`：阈值 0.3、0.5、0.7 的敏感性，读取保存的预测计算。
- `figures/`：报告场景的 PNG/PDF 预览。论文图件仍需按科学问题单独审定。
- `case.rds`：案例配置和结果路径；`completion.rds` 仅表示本次计算与导出完成。

数据准备保存在 `output/.prepared/<物种>/`，训练阶段缓存保存在
`output/.training-cache/`。不要把 `completion.rds` 当作模型或论文质量验收证明。

```powershell
Rscript 04_REPLOT.R output/cast026_glu/Schizothoracinae
```

重绘只读取已保存的比较栅格，不拟合、不预测、不计算响应曲线。
原论文运行目录保留；本入口使用新的结果布局，不读取旧版本缓存或旧的结果对象。

## 可调参数与解释边界

`01_config.R` 的 `protocol` 是训练科学设置的唯一入口。总资源预算默认
12 核、3 个场景进程、8000 MB；比较阶段每进程最多 4 个 RF 线程，训练阶段按
协议最多用 12 个 RF 线程。修改配置开头的 `resources` 即同步控制两阶段资源。
需要逐变量归因时，显式设置 `driver_attribution = TRUE`；仅分析报告场景。

GLU 的中期/末期土地利用快照为 2050/2100，气候代表年为 2055/2085；
HLU 使用与气候相同的 2030/2055/2085 年。映射会记录在场景表中，二者不能在论文
中写成完全相同的时间窗口。时间匹配只使用训练发生记录的年份；没有年份时明确跳过。
裂腹鱼案例为亚科记录集合，不能在论文中称为单一物种。

对比是固定拟合响应下的模型条件场景差值，不识别真实生态因果效应。
`*_main` 为基线简单效应，`*_average` 为平均效应；SD 是模型间差异，不是置信区间。
支持检查使用边际训练分位数。二值结果使用 `>=`；表格计数单位为有效像元，
面积解释须核实等面积投影和像元面积。新版本修复 AOA 最近邻定义，旧 AOA 数值须重算。

本次维护未运行这套完整应用；小型合成数据验证由包测试和 R CMD check 执行。
