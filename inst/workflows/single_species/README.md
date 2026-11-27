# CHELSA 单物种工作流（castSDM 0.49.0）

用实际记录和环境库启动一个项目。这个示例适配 CHELSA 集成气候情景；包的 `cast_run()` 与 `cast_compare()` 可直接使用对齐的表格或栅格。

```r
castSDM::cast_workflow(
  path = project_dir, data_dir = environmental_library,
  occurrence_file = occurrence_csv, species = taxon_id,
  settings = list(input_type = "nsdm_occ"), design_only = TRUE)
```

`nsdm_occ` 的坐标使用环境网格的 CRS。`gbif` 使用其标准 WGS84 字段，并按输入的分类信息和分析掩膜质控。项目保存选择的记录和数据路径，Shiny 和命令行共用该设置。

- 默认使用环境库掩膜，可显式传入对齐的 `domain_mask`。分析掩膜定义预测范围，实际生态可达范围需由记录和生态证据判断。
- 默认只按原生栅格去重。背景从可用真实格点无放回抽样，默认数量由域内合格出现格决定。若提供同类调查点，则可用全部域内去重的调查格点。背景标签代表抽样背景。
- 默认从合格环境清单取样，VIF 在各训练折内执行。用户明确声明的 `predictor_file` 适用于已有研究设计。仅有含义未核实的记录年时，保留实际环境图层年代和来源。
- 自动空间分块由点坐标和类别计数确定。先通过 Shiny 的「检查数据与空间折」核对样本分布、距离和训练边际覆盖，再判断能支持多远的预测。

在配置目录中运行 `00_RUN_ALL.R` 或通过 Shiny 一键运行。设置 `CAST_DESIGN_ONLY=1` 可单独检查点位设计。完整结果位于 `output/cast047_chelsa_ensmean/<taxon>`。

`05_SHINY.R` 浏览结果与按需计算单因子响应。重绘时设置 `CAST_CASE_DIR` 指向已完成的案例，`CAST_FIGURES="03_environmental_driver"` 请求单因子图，再运行 `04_REPLOT.R`。
