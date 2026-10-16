# 先读这里：环境情景如何改变预测适宜度

本次主图情景：ensmean_ssp585_2071_2100。阈值化适宜面积的净变化为 +34.53%；这是模型情景预测，不是已发生的分布变化。
空间交叉验证实际评估 4/5 折；无效折及逐折指标见图 1 和 cv_fold_status.csv。
该情景的非零变化响应涉及 1 个因子：bio05。训练边际范围覆盖率为 67.9%。

按科学问题顺序看四张图：
1. 01_model_validation：空间分组与逐折预测表现；不要只看一个平均 AUC。
2. 02_habitat_change：当前与未来适宜度、空间变化、全部 18 个情景的面积净变化。
3. 03_environmental_driver：实际未来值变化如何影响模型响应；完整变量名与指标在 environmental_driver.csv。
4. 04_support_and_disagreement：训练范围覆盖与模型间分歧；两者均不是置信区间。

Only one predictor changes: total, driver-only and retained-driver contrasts coincide.
当前 Schizopygopsis 使用属级记录；全国陆地背景和出现点缓冲区尚不能代表鱼类可达河流栖息地。
只有一个环境因子变化时，总变化、单因子变化和保留因子对比相同；其余零值表示没有提供该因子的未来变化，不表示该因子没有生态作用。
图由包内绘图函数生成。重绘保留原始计算来源；重绘历史数值结果不等于重新完成拟合实验。
导出成功与视觉检查分别记录。
