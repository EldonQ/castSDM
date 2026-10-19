# castSDM 定位、SDM 方法与源码审查

审查日期：2026-10-16。对象为当前工作树，不只审查新增归因函数。
后续状态：2026-10-19栅格/年份缺陷已修，详见RASTER_REFACTOR_20261019.md；
下列严重程度表保留原审查时状态，不覆盖历史证据。
证据包括用户提供的本地论文全文、对应源码、固定版本的开源源码与示例图，
以及实际病例和已完成的模拟验证。以下区分已修复、部分修复和待处理；
R check 通过不等于生态设计合理，也不证明方法新颖。

## 1. 五个问题先说清楚

**科学问题**：在一个明确的环境变化情景下，物种的预测适宜度在哪里增加或减少，
哪些发生变化的环境因子参与了这些变化，哪些空间解释对数据和模型选择敏感？

**要解决的缺口**：通常 SDM 结果给出当前/未来地图、总体面积变化、变量重要性或响应曲线，
但这些产品没有自动回答“某地点的未来变化与哪个实际变化因子有关”。
全局重要性不能直接代替未来变化解释；静态因子重要，也可能在指定情景中贡献为零。

**核心贡献应定位为**：面向实际当前—未来环境对的空间解释工作流。
同一模型、同一预测尺度和冻结权重，统一输出完整变化、单因子未来值变化、
保留当前值的对比，以及范围标记和模型分歧。拟合、调参和通用解释算法由成熟工具承担。

**别人能不能做**：能。预测变量替换、PDP/ICE、SHAP 和不同情景预测已存在。
用户可以用 biomod2、iml、DALEX 或自己的预测循环实现这些对比。
因此不能把“两次预测相减”或“加一个因果名称”作为新的估计方法。

**为什么仍需要 castSDM**：只有当它显著降低空间对齐、预测适用范围、固定比较目标、
大栅格执行和论文结果组织的负担，并在同一输入上给出可核验的一致结果时，
它才有独立的软件价值。当前创新是可发展的工作流贡献，尚不能声称已证明新的归因算法
或相对所有现有解释工具的优势。MEE 工具包论文应证明具体使用困难得到解决，
而不是把功能数量、代码抽象程度或漂亮案例当作创新证据。

建议中文定位：**物种栖息地变化的环境驱动解释与空间归因**。
短解释：**把未来适宜度变化解释到具体因子与具体地点，并显示解释的适用范围和模型分歧。**
不引入因果角色、DAG、因果发现模块。实际观测因果效应需要另一套研究设计；
当前替换环境值只描述拟合模型在指定情景下的响应。

## 2. 实际参考了什么，而不是只看包名

| 参考 | 已核验与阅读的材料 | 方法/工程上的可借鉴点 | 不应照搬的部分 |
|---|---|---|---|
| ENMeval 2.0，MEE，DOI 10.1111/2041-210X.13628 | 元数据、摘要、R/classes.R、R/plotting.R、全流程 vignette；源码 7e44494460207456d37d783d2956143f04fef8fb | 一个主要评估对象；调参组合、逐折评价和训练/验证结果清楚；响应曲线说明其余变量的固定基线 | 单纯换成 S4 对象不会自动改善 castSDM；固定基线响应曲线不是局地未来变化归因 |
| flexsdm，MEE，DOI 10.1111/2041-210X.13874 | R/fit_ensemble.R:400–416、R/sdm_predict.R:223 起、R/p_pdp.R、R/plot_res.R 和实际示例图；源码 f6635814e3387feef27a725e004c0e16d8e2d23d | 简单可组合列表；models/predictors/performance/performance_part；阈值、clamp、预测区域和分块显式化 | 不把所有模型接口再复制进 castSDM；学习参数契约，不照搬整套函数 |
| blockCV，MEE，DOI 10.1111/2041-210X.13107 | R/cv_spatial.R、R/cv_similarity.R、R/cv_plot.R；源码 51e6903fe652411248aa20937c80b122e5e12e09 | 空间块尺度、折号、种子、类别均衡和分组图是评价设计的一部分 | “空间 CV”名称本身不能替代与转移任务相符的块尺度 |
| Wallace，MEE，DOI 10.1111/2041-210X.12945 | README、数据/背景/分组/转移组件源码；b722ef2a6dcae5d52bd3b519ebbfdab5c6dc849b | 每阶段回答明确问题；配置、记录与转移范围可追踪 | 不为 castSDM 再建大型 GUI 和组件平台 |
| adm，MEE，DOI 10.1111/2041-210X.70074 | 本地全文、adm-0.0.2/R/fit_abund_gam.R:74–92,169–196,260–261、adm_eval.R 和实际 PDP 图 | 科学问题是丰度；分布族与响应类型对应；超参数、逐折表现和返回列表直观；响应图附数据覆盖 | 不把丰度分布族与 presence/background 适宜度混用 |
| MaskSDM，MEE，DOI 10.1111/2041-210X.70200 | 本地全文、train_model.py:26–47、training_helpers.py、data_helpers.py、概念图与 main.pdf | 随机遮蔽训练对应“任意变量子集/缺失数据”的明确困难；与插补及子集重训模型比较；训练阶段预处理与空间分组显式 | 它的贡献不是一般“可解释性”；它有专门训练机制和相应基线，不能据此称 castSDM 的替换预测同样新颖 |
| RISDM，Ecography | 本地全文、R/checkInput.R:14–41、predict.isdm.R | 分布过程与 PO/PA/AA 观测过程分开；PO 数据须有偏差模型；预测量的定义与观测类型有关 | castSDM 没有观测过程模型，不能把机会性记录拟合值写成绝对占据概率 |
| N-SDM / N-SDM 2.0，Ecography | 本地全文、N-SDM-2.0.1/scripts/settings/param_grid.psv 与 2_mainREG/reg_B.R | 大尺度/区域尺度的问题说明为什么需要复杂度；显式参数表与评价轨迹 | HPC 与多层模型复杂度不适合 castSDM 当前单一解释问题 |
| biomod2 / iml / DALEX | 固定版本函数审查见 OPEN_SOURCE_COMPARISON.md | 通用响应曲线、PDP/ICE/ALE、SHAP 和预测回调已成熟 | 当前模拟没有直接比较这些包，不能声称总体胜过它们 |

ENMeval、flexsdm、blockCV、Wallace 的已发表全文/全部论文图未取得，
该部分结论以可访问元数据、源码与 vignette 为依据；本地 MaskSDM、adm、RISDM、N-SDM 全文已读。
SDMtune 属于 Ecology and Evolution，不混称为 MEE。
公开代码也会有历史限制；借鉴的是任务与接口一致性，不是把它们当作无缺陷标准。

实际回看Git早期版本：061e14d（0.1.0）18个R文件、2,900行；bbd6429（0.11.0）
25个文件、7,828行；本次起点734744f为29个文件、10,479行。
0.1.0的cast()依次执行prepare→DAG→ATE→screen→roles→fit→evaluate→predict→CATE，
阶段读起来直观，但把DAG/角色/神经网络/效应识别混进同一SDM仍存在科学设计问题。
0.11.0的README围绕三个响应产品叙述更集中，但“shift保持support”也不是一般成立的保证。
应恢复的是直观责任边界与一条清楚的用户链，不能把旧理论假设整体恢复。

## 3. 按严重程度列问题

P1：可能改变科学解释或产生误导结果；P2：影响可靠性、参数契约、可复现性或大数据执行；
P3：输出可读性、维护成本。没有证据证明的疑点不列为已确认 bug。

### P1：论文结论与方法设计

| 编号 / 位置 | 成因与可能影响 | 可操作改进 / 当前状态 |
|---|---|---|
| P1-01 定位、R/scenario-attribution.R::.cast_attribution_block()，R/effect-grid.R | 替换一个因子后预测差值是拟合响应。相关因子、不完整生态过程和采样偏差影响其含义；容易被误写为真实环境因果效应。 | 已明确 current/future、only/retained 两个基线、非加和性和 remainder 的含义。论文不使用“已识别因果效应”。创新证据仍需直接竞品任务对比。 |
| P1-02 single_species/01_config.R:21,92 与 02_pipeline.R::prepare_species() | Schizopygopsis 是属级记录；全国陆地背景及 300 km 出现点缓冲区不等于可达水系。模型可能主要区分鱼类记录区域与远离水系的全国背景，高 AUC 不保证河流生境有效。 | 用户已确定：属级流程示范，另选物种生态验证。候选与协议见 ECOLOGICAL_VALIDATION_PLAN_20261016.md。物种论文须有流域/河网可达区、偏差控制和独立评价；不能静默换掩膜后继续引用旧结果。 |
| P1-03 single_species/02_pipeline.R 未来栈对齐及 select_ensure_climate | 9 个入选因子中只有 bio05 有未来变化，所有其它因子按当前值固定。单因子/保留因子/总变化必然重合；静态因子零贡献不表示生态上不重要。 | 已在主图和读我文件明确这一事实。需要实际未来水文/土地利用等输入或另一个多变化因子的案例；不能为了图丰富而制造情景。 |
| P1-04 core_method_validation/01_attribution_validation.R 与结果 | 局地对比继承底层模型误差。RF 的有限变化响应噪声在多数情形大于全局 PDP；加性 GAM 在联合响应情形失配。稳定森林种子掩盖数据重采样不稳定。 | 已完成 6 情形 ×20 独立数据，保留负结果；分别显示重采样、种子和模型家族差异。不能声称普遍优越。需要与实现同一目标的现成解释工具比较、更多模型和独立真实验证。 |
| P1-05 R/spatial-cv.R::make_spatial_folds() 与案例 CV | 自动 grid 的尺度主要由 k 和样本空间跨度决定；没有保证块大小对应空间自相关或预测转移距离，部分块只有背景。折被跳过后均值不再是完整 5 折结果。 | 新增 fold_ids 直接接入 blockCV 或流域设计，保留原始折号。案例仍记录 4/5 有效折，不伪称完整 5 折。下一步按评价目标确定物理尺度及类覆盖，并预先报告失败折处理。 |
| P1-06 R/utils.R::compute_omission_at():229、compute_mpa():253 | omission_5/10 是取所有评价样本分数最高 5/10% 后计算遗漏；并非 ENMeval 的训练 presence 5/10% 遗漏规则。MPA 是样本比例，样本不是等面积地理网格。把它叫面积会混淆背景设计、样本比例和地理面积。 | 核心结果只解释明确的逐折 AUC/TSS/Boyce；论文不得直接与 OR10 或 km² 对比。这些内部诊断需要改名/补充 sampling_fraction 元数据，真正的 omission 需训练 presence 阈值再评估测试 presence；实际面积由投影栅格与 cellSize 计算。待独立版本调整，保持旧列数值含义。 |
| P1-07 R/ensemble.R::.cast_ensemble_scores():378 及模型缺失处理 | AUC/TSS/Boyce 的有限组件平均不保证是最优权重；不同可用组件使模型比较改变。普通 ensemble 在每个像元重分配缺失模型权重，可能让当前/未来不再比较同一目标。 | 归因已固定有效模型和权重，任一有效模型缺失时相应比较为 NA。公开显式 weights 接受外部验证结果。需要主分析/等权/单模型敏感性；不能把近似相等权重当作科学发现。普通 projection 的逐像元集合仍需在含缺失模型的情形审计。 |
| P1-08 R/predict-spatial.R::.cast_aoa_di():361，spatial-cv.R:529 起 | 所谓 DI 是标准化空间的最近距离，未完全复现经典 AOA 的归一化；不同外折选出不同变量时将其距离分布合并，后续全数据变量空间不一定可比。边际范围盒也不能发现联合组合空洞。 | 已禁止插补行声明观测范围可靠；主图仅称训练边际范围覆盖和模型分歧。AOA 应委托 CAST 或在固定 predictor/transform/reference 下核验，不把当前自定义 DI 等同于文献原方法。待专项校准。 |
| P1-09 single_species/09_convert_schisoft.R 的 year、locality 字段与原始 sampledDate/remarksOfSampledDate | 转换只取日期前4位，丢失发表/接收日期替代采样日期的标记；referenceInEnglish缺失时locality也缺失。不能把所有年份当作采样时间用于时间验证。 | 已新增从原始文件读取的物种候选审查：S. pylzovi有21条日期替代，S. thermalis有656条无可追踪文献。时间验证排除替代日期并逐条追溯；原始记录不删除。属级旧流程仍仅作示范，不冒称时间独立验证。 |
| P1-10 R/nsdm_castsdm_adapter.R::.nsdm_period_dir 与未来非气候层 | 代表年与平均时期必须有来源；把土地利用快照2050/2100当成气候2055/2085的同一时期也会混淆输入含义。 | 已核验E:/CHECO26/N-SDM/00_docs/README_data.md:49–50，数据说明确实给出气候源时期2021–2040/2041–2070/2071–2100，现有标签有本地元数据依据；不是已确认的气候标签bug。土地利用仍是近似配对，当前run保持静态；新物种引入未来层时须附产品来源与配对规则。 |
| P1-11 R/ensemble-raster.R 掩膜读取/几何失败与 ok=complete.cases 分支 | 掩膜不匹配或某块读失败时，警告后预测未掩膜区域。这是旧文档明确的宽松策略，但在范围必须受限的科学工作流中会改变预测域。旧mask约定是非NA，因此0有效并不违反其文档；它与归因mask!=0的契约不同，复用同一0/1掩膜易产生不一致。complete.cases也允许Inf流入旧预测路径。 | 当前示范用同网格NA掩膜和有限输入，验收会检查预测/归因一致。下一版应有明确严格模式，掩膜失败即停止；统一0/1和NA规则需兼容迁移，不能静默改变旧参数。补充失败/0-1/NA与非有限环境测试。不能把这一处改变混入正在运行的计算来源。 |

### P2：参数、源码与执行

| 编号 / 文件、函数或逻辑 | 成因与影响 | 修复建议 / 当前状态 |
|---|---|---|
| P2-01 R/model-fit.R::cast_fit() 与 .cast_tune()；cast.R::.cast_tuned_args() | 指定 rf_mtry 曾在 tune 时被网格替换，untuned 全数据 refit 又可能回到默认；评价和最终模型配置不同。 | 已修：RF 固定参数进入网格；从真实拟合对象读取 refit 配置；CV 显式透传 mtry/depth/MaxEnt 参数/线程。回归验证实际每折 mtry。 |
| P2-02 R/model-fit.R::cast_fit()、R/spatial-cv.R::cast_cv() | fractional k/n_repeat/树数等过去可能被截断；重复/空模型集合没有清楚错误；单类/缺失响应、全非有限预测变量流入 engine。 | 已修：整数与范围前置校验、case-insensitive 重复模型拒绝、拟合响应双类与无缺失要求、无有效预测变量报错。有限训练值中位数用于插补。 |
| P2-03 R/predict-spatial.R::cast_predict():99,169–178 | n_imputed 算出后丢弃；插补后再计算环境覆盖，把缺失环境误称正常支持。 | 已修：输出逐行 n_imputed，原始环境有插补的 MESS/AOA 标记为 NA。普通表预测仍可插补；配对归因与栅格保持缺失，所以不同用途的策略必须说明。 |
| P2-04 R/cast.R::cast() 与 spatial-cv.R::cv_one() | 请求 CV/refit/ensemble 失败后继续发布局部结果；全模型没有有限预测仍可能标记 evaluated。 | 已修：顶层必需阶段失败即停止；无有限模型预测标记 no_predictions；案例完成门控包括 model_error/no_predictions。待改进单模型部分失败的统一状态和错误表，避免仅在日志中追溯。 |
| P2-05 R/model-fit.R::fit_traditional(gam):388 起 | smooth 模型失败后自动重拟合 linear-only 模型，仍称 GAM，没有提示；用户可能解释不存在的非线性结构。 | 部分修：fallback 时明确警告并指出检查 formula。下一步输出 effective_formula/fallback_reason 或改为显式用户参数；此次保留兼容，未直接取消合法线性拟合。 |
| P2-06 R/model-fit.R::.tune_rf/.tune_brt/.tune_maxent | RF OOB、BRT 随机内部 CV、MaxEnt 随机分层 CV 使用不同评分/阈值路径；BRT 默认网格忽略显式 depth/shrinkage，与“所有显式参数固定”的总说明矛盾。 | 已修说明，BRT tree-count fold 数与 tune_folds 对齐，gbm 的 cv.fitted 经 plogis 转为响应尺度。源码核验确认旧 TSS 是排序阈值，单调变换通常不会改变它的最优 TSS，不能夸称该修复改变旧选优。尚不是统一空间调参；长期交给外部工具，短期公开 grid 和实际选择。 |
| P2-07 R/as-fit.R、model-evaluate.R::predict_single_model() | 原 new_cast_fit 构造器只有内部包装格式，不能证明跨工具互操作；添加模型要修改内部 dispatch。 | 已修核心边界：cast_as_fit(models, reference, predict_function)，检查一行一个 [0,1] suitability。已实际调用iml::Predictor，三个配对对比与直接GLM最大差0；blockCV原始3折完全保留。见03_external_interoperability.R。复杂预处理归回调负责，互操作不等于优越性。 |
| P2-08 single_species/01_config.R、02_pipeline.R::run_species() | 案例不显式声明树数、mtry、depth 等，fit/CV/refit 独立使用隐藏默认，未来切换 tune 很容易失配。 | 已修：配置显式列出 4 engine 参数、tune_folds 和线程；拟合、CV、refit 透传；全数据 refit 复用真实训练 fit 配置。默认设置不等于推荐普适生态参数。 |
| P2-09 R/ensemble-raster.R:266,315–320；effect-products.R:324；predict-tiled.R:224 | 预测输入分块但输出仍分配整幅栅格数组；并行 worker 放大 RAM。max_memory_mb 是输入块估算，最小 10 行也可超过预算。 | 未全部修：新归因栅格核使用写入块与在线摘要，旧 projection/effect/tiled 仍需逐一流式改写及实际内存基准。已有默认关掉耗时单一 SD 全图诊断，保留显式 opt-in；不能宣传全包严格内存上限。 |
| P2-10 R/utils.R::.cast_write_raster_atomic() 与 scenario-attribution.R 文件发布 | temp+rename 是正常路径，但 rename 失败 fallback copy/overwrite 不保证原子性；Windows 目标文件存在时发布规则需核验。 | 已修NEWS过强的“all atomic”表述，明确rename与copy fallback的区别。现有 tests 覆盖中断/空文件等，仍需Windows重复覆盖、跨盘和失败恢复专项验证。 |
| P2-11 single_species/03_figures.R；R/workflow_provenance.R::workflow_source_signature() | 数值计算与所有绘图脚本共用签名；仅改字号会迫使重跑拟合。反过来移植新图又容易把旧计算误标为新代码。 | 部分修：重绘接口必须出示与 bundle 完全相同的原计算 receipt，并另外记录 plotting_source。旧数值结果不改来源。下一步数值/绘图/输入分别签名并维护依赖边界，保持源校验。 |
| P2-12 R/model-fit.R 与多个 set.seed() 路径 | RF seed=NULL 仍使用42，且全局 RNG 被修改；研究者以为不设种子就是随机，或其它后续计算序列改变。 | 未修默认：应公开 resolved_seed，或用局部 RNG 上下文；统一 NULL 语义需要兼容性版本及复现检查。不能仅以“有 seed 参数”声称完整复现。 |
| P2-13 R/project-future.R:755–757 | Albers 网格下 current_centroid_lon/lat 等实际是米制 x/y，字段名错误；位移 km 的计算本身正确。 | 记录该限制，不把米制值解释成经纬度。后续增加 native_x/y/CRS 和真正 lon/lat；旧字段迁移需兼容性说明与测试。 |
| P2-14 R/utils.R::.cast_thin_points(order='year_desc') | 缺失year经ifelse变成-Inf，升序排序会优先保留日期未知记录，和“最新优先”契约相反；缺失年记录可能抢占同地点较新记录。属级当前run使用random分支，不受此排序路径影响。 | 已定位，待数值版本修复为未知年最后并增加同位置已知/缺失年的回归检查。记录年代较新也不自动说明坐标较精确；不据年份替代不确定度审查。 |

### P3：用户输出与维护成本

| 编号 / 位置 | 问题与影响 | 当前处理 / 后续 |
|---|---|---|
| P3-01 single_species/03_figures.R 默认报告 | 旧案例 158 张图围绕多个不同目标，没有论文主线，重复面板难读；固定 1SD 对比与实际未来变化混排。 | 已改成 4 主图：评价→变化→变化因子解释→范围与分歧；全 18 情景数值仍输出；legacy audits/permutation 默认关闭但保留入口。原历史数值与 inst 资产保留。 |
| P3-02 R/plot-effect-map.R::.cast_display_raster() | bilinear 用于二值范围图制造 0–1 之间的显示值；巨大 colourbar、反复图例、机器变量名和无说明的全域限值掩盖空间结构。 | 已增加 near 显示方法，二值主图显式使用；共享色标与对称变化色标、记录 p99.5 截断；压缩图例、替换主图标签；原始像元值不改。 |
| P3-03 R/methods-plot.R 效应摘要 | coverage 分母取最大支持记录数而非评估 n，会夸大覆盖；重复色阶替换导致警告和样式漂移。 | 已修分母并回归验证；移除重复色阶，增加/减少统一蓝/橙。 |
| P3-04 cast() / cast_fit() / cast_cv() 与多条效果 API | 实际参数分别为31/17/新增后28；主 API 混合数据、筛选、模型、验证、解释和导出。并不是再加一个统一 config 类就解决。 | 核心使用收敛为已有模型→cast_as_fit→cast_attribution→plot；旧拟合 API 为兼容保留。下一次移除依据调用图、使用频次和迁移清单，而不是整目录删除。|
| P3-05 tools/workflows 镜像与活动脚本 | 镜像过时、个人绝对路径、多个历史入口让用户难知当前步骤；删除 inst 又破坏绘图后续。 | 维护一个单案例 runner 和一个方法验证 runner；同步源码镜像与结果指针。恢复示例 CSV 并保留全部地图资产。历史证据不能按当前代码重新命名。 |
| P3-06 single_species/03_figures.R::short_labels 与外部数据窗口口径 | “Surface-water cover (1 km)”把_f1000窗口参数当作真实圆盘半径/覆盖尺度；vege_f20000同样标成20 km，但真实方格41×41 km。 | 已核验数据说明README_data.md:34–46：f1000为3×3 km方格，f20000为41×41 km方格。应在显示标签与变量字典明确窗口类型，网格分辨率1 km单独声明。当前正在运行的源码保持冻结，结果交付说明记录该标签限制。 |

## 4. 最小源码边界与默认结果契约

最小核心链：`external model + prediction callback + training reference`
→ 对齐当前/未来环境 → 一份配对数值核 → 空间图与清晰摘要。
当前 core 入口为 `cast_as_fit()`、`cast_attribution()` 和 S3 plot。
保留基础 validate/prepare/projection/ODMAP 辅助；covsel/双阶段筛选、
RF knockout、1SD audit、permutation 作为有明确用途的可选诊断，不混成核心创新。
不另建 engine 插件系统、因果角色对象、多级计划类或万能输出管理器。

当前源码实测为34个R文件、11,980行、181个函数定义；模块清单见
source_module_inventory.csv。methods-plot.R 1,112行、var-select.R 1,025行、
ensemble.R 864行、project-future.R 835行。此次并没有实现全源码行数下降：
covsel和归因等此前新增功能仍在工作树。已经完成的是核心入口与默认结果收敛、
参数bug修复和外部工具边界；legacy数值路径退役仍需要单独的调用/兼容性审查。
不把“输出图变少”误称成全部代码质量问题已解决。

已生成function_call_inventory.csv和internal_call_edges.csv，记录181个函数的
源码行号、参数数量、公开/S3状态与直接调用边。静态零调用的内部函数为
.cast_support_fraction、.cast_aoa_di_train和.cast_thin_points；最后一个仍被病例脚本调用。
前两个须进一步检查测试、示例和下游使用后才能退役。动态dispatch/序列化回调不在
静态解析范围内，不能把零直接调用等同于可删。这份清单用于后续责任合并与迁移，
不据此删除inst或公开API。

06_legacy_contract_audit.R已实际复现三个未修路径：0掩膜格在旧预测中有效而在配对归因中
被排除；Inf环境值经旧complete.cases传给预测回调而新归因丢弃；同位置year=NA在
year_desc中优先于2020记录。四格栅格示例、结果CSV与来源清单随审查保留，
不是仅凭代码样式做判断，也不把已定位的问题列为已修。

| 默认结果 | 它应回答的科学问题 | 必须附的证据 |
|---|---|---|
| 01_model_validation | 这个模型的预测评价有多可信？ | 原始空间折号、逐折指标、缺失/失败折、背景设计。AUC 仅证明相对于该背景的排序 |
| 02_habitat_change | 当前/未来预测在哪变化？ | 同尺度连续地图、总变化与各情景净面积；阈值、面积单位、mask 与时期。净增不能代替 gain/loss 明细 |
| 03_environmental_driver | 实际哪些变化因子参与预测变化？ | only/retained 的基线、完整 summary、静态因子的解释、非加和 remainder；不能把 mean_signed 抵消误称无影响 |
| 04_support_and_disagreement | 哪些解释依赖超范围或模型选择？ | 范围覆盖与原始比较分开；模型 SD 不是采样 CI；需另做重采样稳定性 |

原始完整栅格、4个模型、18个情景、CV 表与来源收据保留。
“四张主图”是阅读入口，不意味着只允许四种数据或删除绘图资产。
所有地图显示重投影与色标截断要声明；源图由包内绘图函数生成。

## 5. 论文验证计划与目前证据边界

当前完整模拟包含线性、非线性、联合响应、零因子、相关因子和超范围6种情形，
每种20组独立数据。基线是同拟合模型、同目标因子的全局边际 PDP 和统一平均变化量；
这些基线没有覆盖所有 ALE/SHAP，也不是生态观测因果真值。
参考 2026-10-16 04:14 UTC 已完成运行（complete_production_20261016T041431Z）：
120组、492项数值检查、0运行错误，收据与当前包数值源码一致。
四张实际PDF页面已由R渲染并视觉检查，导出表和视觉QA单独记录。
GAM 局地 RMSE 在单调/非线性/相关/超范围为0.0085/0.0225/0.0257/0.0412；
对应全局 PDP 为0.0250/0.0442/0.0291/0.1012。RF仅联合响应对 PDP 有平均优势，
多数情形更噪。森林种子地图相关0.974–0.996，数据重采样0.008–0.787。
不能只展示成功子集；对未知模型结构的优势需要进一步验证。

下一项论文实验证据应是一个**同任务、同底层模型的跨工具工作流对比**：
手工成熟工具方案与 castSDM 均计算同一对 actual future contrasts，
检验数值一致性、对齐/缺失/权重变化错误的发现率、运行内存和配置负担；
对 SHAP 等不同 estimand 的比较必须先说明问题差异，不能用 RMSE 排名替代。
同时选择明确物种、可达背景、独立时期/区域的案例，检验未来转移与空间解释稳定性。
这些才会回答“别人能做，为什么还需要本包”。

已补第一项同任务实测：05_workflow_contract_comparison.R通过iml预测接口手工计算
2,500行、两个变化因子的配对对比；正确对齐后所有数值最大绝对差0，分块栅格也一致。
未来行打乱后，显式配对手工方案和castSDM仍一致；没有配对的示例最大差0.8598，
这是示例漏做空间检查而非iml缺陷。castSDM拒绝错位栅格与无法匹配坐标。
这支持默认空间契约的工程价值；5次小规模计时不支持速度优势，峰值内存与用户负担
仍未实测。没有把不同估计目标的ALE/SHAP拿来做不公平数值排名。

用户已确定属级示范与物种验证分开。已实际审查原始SchiSOFT候选记录，
暂以S. pylzovi为第一候选：144条、119个1 km网格、40篇可追踪文献，
但排除发表/接收日期替代后，2018年起只有5个网格，尚不足以宣称时间生态验证。
S. thermalis虽然741条，只有97个网格且656条无可追踪文献，不能按记录数选它。
详细协议见core_method_validation/ECOLOGICAL_VALIDATION_PLAN_20261016.md；生态验证未运行。

工程优先顺序：已完成参数契约与互操作边界 → 结果契约与案例 QA →
固定变量空间 AOA/标准 omission → 三条 legacy 栅格写出统一流式实现 →
依据调用图迁移/退役冗余接口。每次只改变一个可验证责任，不靠删除文件数衡量瘦身。

当前代码的短 R check 已验证安装、文档、examples 与全部测试，0 error/0 warning/0 note。
每次来源变更后以新运行 receipt 指向的 evidence 为准；此文不把历史计算冒称最新源码。

## 6. 本轮最终验收状态

Schizopygopsis完整运行已完成：532个保留出现点、10,000背景点、9变量、4模型、
18投影与归因栈。4/5空间折有效；主图情景的阈值化净变化+34.53%，训练边际范围
覆盖67.95%。仅bio05改变，静态变量对比为零；不能把属级全国背景预测写成鱼类扩张。
所有18情景接受抽样投影/归因数值检查，原始栅格和完成收据均保留。

实际PDF视觉检查后修复了无效CV折标注重叠、两个方形窗口标签和二值地图连续图例。
修订图放在独立plot_revision目录，计算和绘图身份分别记录；没有改写旧receipt。
只有methods-plot.R中的plot.cast_cv、plot-effect-map.R和工作流03_figures.R变化，
其余计算源码与完成记录一致；最新R check为0 error/0 warning/0 note。

内存抽样窗口内，7进程合计最大工作集9.87 GiB。采样只覆盖归因阶段的一部分，
不是全程峰值。25 MB块输入预算不应被文档描述为整个模型/进程的内存上限。
尚未修复的legacy缺陷和独立生态验证仍按前表列为未完成，没有据此宣布全面瘦身完成。
