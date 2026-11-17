# PGR 贝叶斯自适应分组试验

首个配置对比 MFM–LMM、DPM–LMM 和 PPMx–LMM，采用论文 Case 1、真实生成组数 (R=10)、随机截距方差 (0.5) 的 RI 设定，重复 **5 次**。三个方法在每次重复内共用 RGLL 生成的 (7500/2500/2500) 条训练、参考和测试观测。`config.pilot.json` 的 `experiment` 可换成论文 Case 1–7 网格中的另一个单一设定。选择 RI 作为首个设定便于先检查后验混合及运行时间。

## 准备与执行

依赖 Python 3.11、NumPy、SciPy、Pandas、PyTorch，以及原 RGLL 的 R、Rcpp、RcppArmadillo、ClusterR、MASS 和 `sourceCpp` 编译环境。采样使用 CPU，多链可并行执行。`--device auto` 在 CUDA 可用时将硬分组 plug-in IMSPE 的 LMM 拟合和矩阵运算放在 GPU。后验条件 IMSPE 使用 CPU Cholesky。

以下命令从含有 `experiments` 的项目根目录执行。代码会从当前 RGLL 仓库根目录或本地工作区的 `RGLL` 子目录读取原始生成器。也可用 `--rgll` 指定其目录。`RSCRIPT` 应指向本机的 Rscript。

```powershell
$env:RSCRIPT = 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe'
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'

python experiments/bayesian_regrouping/test_bayesian.py
python experiments/bayesian_regrouping/run_experiment.py --output experiments/bayesian_regrouping/results/case1_ri_5
```

恢复时使用相同的命令和输出目录。每条链保存重新启动所需的采样参数、分区和随机数状态。已有完成结果自动复用，输入数据指纹、配置或设备改变时需要另一个输出目录。两条链分别从单一占用组和模拟的先验分区出发。初始分区与真实 (R) 独立。

敏感性配置包括 `mfm_support`（总分量数先验支持上界 128）、`mfm_prior`（几何先验概率 0.15）、`dpm_alpha`（浓度先验 rate 0.5）和 `ppmx_similarity`（协变量模型的组内方差 1.5）。通过 `--variant` 选择一个变体，并指定独立输出目录。每个变体只运行相应模型。

## 结果

- `replications.csv` 和 `comparison.csv` 保存逐次结果、均值、以五次重复为单位的 95% t 型 Monte Carlo 区间、成对差异和有效配对数。
- `posterior_K.csv` 和 `posterior_summary.json` 分开报告占用组数、MFM 总分量数、参数及随机效应尺度。MFM 同时报告总分量数先验与支持上界的后验质量。
- `posterior_trace.csv`、`sampler_diagnostics.json` 和 `new_group_trace.csv` 保存多链标量诊断、随机效应与组大小摘要、DPM 浓度轨迹、PPMx 相似性和 cohesion 轨迹，以及每个 draw 在测试协变量上的平均新组概率。`chain_*/draw_*.npz` 保留完整标签与参数，组编号以单个 draw 为单位解释。
- `imspe_summary.json` 分开报告共同的硬分组 ML-EM plug-in IMSPE，以及在固定后验方差分量下的条件 IMSPE 子样本。`test_predictions.csv` 保存硬分组预测、含新组权重的后验预测均值和标准差。`representative_partition.npz` 包含代表分区及其参考和测试权重。
- 如果代表分区只占用一组，硬分组 ML-EM 指标标记为不可辨识，仍可报告基于后验先验正则化的条件 IMSPE 与预测。数值失败和达到 EM 上限分别标记。资源统计见 `runtime.json`，疑似混合不足的链保留原始轨迹供核查。

模型公式、预测方式和先验说明见 [IMPLEMENTATION_NOTES.md](IMPLEMENTATION_NOTES.md)。首个配置和验收测试已编写，运行和数值核验等待执行。
