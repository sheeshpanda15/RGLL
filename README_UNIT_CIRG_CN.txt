CIRG-U (Unit-level CIRG) v5
===========================

核心变化
--------
1. 不再对单个 observation/row 聚类。
2. 原始 repeated-measure unit 是不可拆分的最小重分组对象。
3. RI unit embedding:
      pooled residual mean + PCA(unit-wise X means)
4. RS unit embedding:
      ridge unit intercept/slope proxies + PCA(unit-wise X means)
5. K-means 只用于生成初始 partition。
6. 模拟退火直接搜索 whole-unit move / swap / merge / split。
7. 每个 candidate partition 仍然使用 RI/RS 对应的 empirical IMSPE 评价。
8. reference set 只用于 target moments；不使用 test y。

方法名
------
CIRG-U   : 新的 unit-level 方法（建议作为 proposed method）。
CIRG-ROW : 旧 observation-level 方法，仅用于 ablation。

Smoke test
----------
cd /users/k21181837/RGSS
Rscript smoke_test_cases10_11.R

首先检查：
- CIRG-U unit_consistency_mean 应严格等于 1；
- CIRG-ROW 通常小于 1；
- CIRG-U 的 Var_a_hat 不应像 CIRG-ROW 一样系统性接近 0；
- CIRG-U MSPE 不应继续机械地等于 LM。

小规模重复实验（建议先 nloop=10）
----------------------------------
METHODS='ORACLE,OBS,LM,KM,CPF,BLM,CIRG-ROW,CIRG-U' NLOOP=10 Rscript run_case10.R
METHODS='ORACLE,UNIT,OBS,LM,KM,CPF,BLM,CIRG-ROW,CIRG-U' NLOOP=10 Rscript run_case11.R

正式实验默认只跑 CIRG-U，不再默认跑 CIRG-ROW。
