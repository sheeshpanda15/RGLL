# 本轮真实数据中的KM

## 方法定义

KM是“仅按X做MiniBatch K-means分组，再拟合RI-LMM”的基线，不是单纯用K-means直接预测y，也不是Bonhomme--Manresa分组固定效应估计器的完整复现。

1. 类别按训练集字典作数值编码；连续和编码后的X按训练集min/max缩放，训练样本落在[-1,1]。新数据沿用同一变换，不重新拟合编码或缩放。
2. 只用训练X聚类，不用训练y决定分组，不用参考/测试响应选择K。
3. 每个候选K调用ClusterR::MiniBatchKmeans，选择最小化下式的可用聚类：

```text
C(K) = n * log(max(WSS(K) / n, 1e-12)) + K * p * log(n)
```

n为训练记录数，p为训练特征数，WSS按实际返回分组内的样本均值重新计算。这是代码自定义的BIC式聚类准则，不是下游RI-LMM似然的标准BIC。

4. 候选K若产生小于tau的簇，逐步减K重做；若连K=2都无法通过，该候选跳过。准则使用实际返回的K，而不是原请求值。所有候选失败则该次KM记失败。
5. 在所得组标签上拟合与其他方法相同的RI模型。测试X按最近的训练质心硬分配到组，再作RI预测；不是PGR的高斯软分配。

## 参数

| 参数 | 本轮设置 |
| --- | --- |
| K搜索范围 | 2到floor(n/tau)，最多取8个等间隔整数候选 |
| 最小训练簇大小 | tau=5(p+1) |
| 聚类初始化 | kmeans++ |
| num_init | 1，即每次调用一次初始化，不是50次重启 |
| max_iters | 25 |
| batch_size | min(1024,n) |
| init_fraction | 1.0，库默认 |
| early_stop_iter | 10，库默认 |
| tol | 1e-4，库默认 |
| 内部seed | 未显式传入，使用ClusterR默认1 |
| RI-LMM拟合 | ML-EM，最多300次，容差1e-5 |

本机真实数据运行使用ClusterR 1.3.2。外层代码虽然调用set.seed，但MiniBatchKmeans有独立seed参数；不能把外层R种子当成已传入内部的种子。

| 数据和设置 | p | tau | 请求的K候选 | 实际选出K的均值 |
| --- | ---: | ---: | --- | ---: |
| Student Full | 30 | 155 | 2,3 | 2.00（49次） |
| Student Excl | 21 | 110 | 2,3,4 | 2.00（50次） |
| NLSY79 Full | 18 | 95 | 2,6,9,13,16,20,23,27 | 6.40（50次） |
| NLSY79 Excl | 14 | 75 | 2,7,11,16,20,25,29,34 | 7.14（50次） |

候选K可因最小簇大小检查而下调，所以最终K不一定是请求网格中的一个值。原变量分组对照直接采用显示的类别或联合类别，不使用KM的聚类搜索和tau筛选。

## 初始化诊断更正

此前对Student Full第48次使用30个不同外层R种子的检查，都沿用了内部seed=1，因此不是30套不同的内部随机初始化。2026-11-11直接改变内部seed=1,...,30并仅检查训练X的簇大小：内部seed=2可找到满足tau=155的两簇分区。这说明原失败属于当时初始化/约束检查的结果，不是“不存在可行分组”的证明。

新诊断只使用训练X，不依据测试MSPE选初始化；未修改实际拟合代码、补齐第48次预测或替换已有结果。表格继续使用Student Full共同49次，其余50次。

## 分组字段与代码位置

- Student：studytime、reason、goout、parent_edu及六个两两交叉；parent_edu=max(2,Medu+Fedu)。交叉按类别元组生成组标签，不是把数值相乘。
- NLSY79：Region=R8496500，Worker class=R7898500，Education=R8497000，以及三个两两交叉和一个三重交叉。
- R7898500是第一份工作的class of worker，并非雇主身份标识。官方说明：https://www.nlsinfo.org/content/cohorts/nlsy79/other-documentation/errata/errata-1979-2010-data-release
- 候选准则和网格：ROG_method_comparison.R中的fit_kmeans_bic、thin_integer_grid。
- 聚类参数、最小簇检查：ROG_aligned_fixed.R中的cluster_x。
- 真实数据调用及测试组分配：real_data_experiment.R中的fit_real_method。
