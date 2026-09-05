RGSS/CIRG observation-level comparison v6
=========================================

原则：所有会产生新分组的方法都必须对 observation/row 进行分组。
不允许以原始 unit/group 为不可拆分对象。

主方法：
ORACLE  true observation labels，仅模拟上界
OBS     错误/过细/粗糙的 observed observation labels
LM      不分组
KM      X 上 observation-level K-means
GMM     X 上 observation-level Gaussian mixture（需要 mclust）
CPF     observation-level concave pairwise fusion adaptation
        - 每条训练 observation 一个 latent intercept
        - 小样本使用完整 pairwise graph
        - 大样本为可计算性使用 residual-neighbour sparse graph，因此正式文中必须称 computational adaptation
        - test y 不可用，因此训练 X->CPF-label classifier 给 test observations 分组
MIXREG  observation-level latent-intercept finite mixture regression
        - train 用 (X,y) 识别 latent classes
        - test 只通过训练得到的 X->class classifier 分组
CIRG    你的 observation-level X-only + IMSPE 方法

BLM 为什么删除：
Bonhomme-Lamadon-Manresa (2022) 的 two-step GFE 明确先为 panel individual
构造 informative moments，再对 individuals 做 kmeans。单条 observation 没有对应的
individual-level moment；强行用 row residual 替代后已经不是 BLM。因此不放在主表。

运行：
Rscript smoke_test_observation_v6.R

正式：
MODEL=RI NLOOP=10 Rscript run_case10_observation.R
MODEL=RI NLOOP=10 Rscript run_case11_observation.R
MODEL=RS VAR_B=0.1 NLOOP=10 Rscript run_case10_observation.R
MODEL=RS VAR_B=0.1 NLOOP=10 Rscript run_case11_observation.R

先跑 smoke，再 nloop=10；不要直接 nloop=50。
