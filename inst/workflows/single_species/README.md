# CHELSA 单物种工作流（castSDM 0.42.0）

规范源码在 `cast/inst/workflows/single_species`；外部工作目录是部署副本。
请在重启的 R 会话运行。入口检查并使用对应版本的隔离 R 库；环境库仅被读取。

## 运行

工作目录的 `workflow_settings.rds` 持久保存明确选择的物种及记录文件，重启 R 不会丢失。入口和 Shiny 使用同一个选择；未配置时不会默认运行裂腹鱼。若旧的 `CAST_SPECIES` 与已保存物种冲突，入口在准备数据前停止并说明冲突。

配置或切换物种只需调用一次（`run=FALSE` 只检查输入，不训练）：

```r
castSDM::cast_workflow(
  "E:/Package/castSDM_workflows/single_species",
  "E:/CHECO26/N-SDM/data",
  "E:/Package/castSDM_workflows/species_occ/Callosciurus_erythraeus.csv",
  species = "Callosciurus_erythraeus",
  settings = list(input_type = "nsdm_occ", domain = "national"), run = FALSE)
```

之后每次都使用同一入口，无需额外环境变量：

```r
source("E:/Package/castSDM_workflows/single_species/00_RUN_ALL.R", encoding = "UTF-8")
```

物种名不决定坐标格式或分析域；应按实际数据明确指定。

该文件按环境网格 CRS 的 `lon,lat,year` 读取；来源、采样努力和年份含义未经独立核实，不自动标记为审核后的 GBIF 记录。松鼠默认使用全国陆地掩膜，与鱼类的河流局地范围不同；两者不能仅按 AUC 排名。两例均使用 GLU 当前土地利用指标。
新结果写入 `output/cast042_chelsa_ensmean/<物种>`，不覆盖既有研究结果。
`CAST_DRY_RUN=1` 只做输入检查；正常运行前应取消该环境变量。

## 流程与解释

质控与 1 km 稀疏化 → 明确分析域和缺失策略 → 空间留出 → 训练折内 VIF≤10 → RF、MaxEnt、GAM → 五折空间验证 → 全数据重拟合 → CHELSA 气候预测与环境响应地图。

- 空间留出按整块、分别考虑两类记录的目标比例，报告 `split_class_balance.csv`。若某一空间块集中了超过训练比例的存在记录，目标比例本身不可达，会明确警告；不会拆散空间块凑数。
- VIF 只检查线性冗余，不规定变量数量，也不证明生态相关性。
- 背景点案例设置为 10,000；背景标签不等于观测缺失或真实缺席。
- 默认完整样本策略不插补。验证群体在分割前固定，不在查看留出结果后追加记录。`excluded_cohort_final_predictors.csv` 明确哪些预先排除的记录在最终变量上完整；改变候选库或缺失政策应另开分析。
- 鱼类分析域保持河流相关网格与记录周围 10 km 范围的交集。域外没有预测；网格覆盖面积不等于实际河道水面面积。松鼠的全国背景也不自动构成其生态可达范围。
- 当前气候为 CHELSA 基线，未来仅为 ensmean 的 SSP126、SSP370、SSP585 × 2011–2040、2041–2070、2071–2100；其他变量静态。平均气候上的预测不等于各 GCM 预测的平均，也不包括 GCM 间不确定性。
- 环境地图表示固定模型的适宜度响应。单因子响应不必相加为总变化；边际范围支持不等于联合支持或因果识别。加性 GAM 和 lq MaxEnt 的背景方向一致可能来自模型结构。
- 使用现有空间训练折的拟合，在同一批原始留出位置、真实未来协变量上检查解释敏感性，不再增加五次拟合。`diagnostics/response_spatial_sensitivity.csv` 报告共同支持位置、各折变量入选、全零响应、方向一致性与幅度变化；不是独立重复、置信区间或因果效应。

## 查看与重绘

运行 `05_SHINY.R` 浏览结果。高级设置允许填写有生态依据的空间块边长（米）；留空保留自动估计，不设置人为默认边长。匹配的出版响应栅格会被复用；地图统计使用原生共同支持网格，记录位置统计单列。预览按有效网格聚合，支持域在聚合前应用；点击读取原始网格。

默认保留八张主图、既有地图样式和 PNG/PDF。ALE 使用全部真实参考位置，只默认绘制情景改变的拟合变量。分类图、全部响应曲线及完整图集按需计算：

```r
Sys.setenv(CAST_FIGURES = "all") # 或 responses、maps、drivers
source("E:/Package/castSDM_workflows/single_species/04_REPLOT.R", encoding = "UTF-8")
```

`workflow_timings.csv`、`figures/phase_timings.csv` 分别记录流程、图件数值准备和渲染时间。
`predictor_qualification.csv`、`predictor_missing_qc.csv`、`predictor_semantics.csv` 保留覆盖、记录损失和来源尺度依据；`gam_concurvity.csv` 报告 GAM 响应替代性。
未完成比较不生成完成标记；同配置重跑校验并复用成功阶段。不要通过改名旧结果假装兼容新分析。

裂腹鱼来源为 [SchiSOFT v1](https://doi.org/10.6084/m9.figshare.24638538.v1)，CC BY 4.0。保留原始年份并审计采样日期、区间及替代值，不将不明确日期转换成精确采样年。

原始环境栅格目录和输入清单可以共用；准备数据目录和训练缓存按物种隔离，并继续核对记录文件、分析域、数据内容及计算规则。1 km 稀疏、10,000 背景、95% 覆盖和 VIF 10 是声明的默认设置，不是从物种自动学习的生态最优值；入口会打印，主要设置可在 Shiny 或 cast_workflow() 中修改。
