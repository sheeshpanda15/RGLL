# 真实数据实验

论文的真实数据表使用 Student 和北京空气质量数据。两张表均比较 RI 和 RS，
按 Excl、Full 分面。Full 另列按观测变量分组的模型。

## 复现

在项目根目录执行：

```bash
DATASETS=student MODELS=RI,RS INCLUDE_ADMIN=true NLOOP=50 N_CORES=4 RESULT_DIR=result/student_paper Rscript run_real_data.R
Rscript beijing_experiment/run_beijing.R
```

Student Full 的第 48 次重复因 KM 和 PGR 拟合失败而排除，表格使用 49 次共同有效
重复。Student 运行脚本保存全部结果后，会因这次已知失败返回非零状态。检查
`failures.csv`，并用所有所列方法和两个模型共同有效的重复比较。北京数据、
依赖和表 2 参考汇总见 [北京实验说明](beijing_experiment/README.md)。

## 数据与历史来源

`real_data/data/` 保留原始的 `student-mat.csv`、`student-por.csv`、`2004.csv`。数据从本项目Git历史 `0f45646` 直接恢复，不替换为相似的外部版本。`prepare_real_data.R` 仅在输入缺失时从该版本恢复CSV，不恢复旧的算法文件；复制到HPC时需同时复制数据目录。

对应的历史实验脚本是 `HPC test/true dataset/student/truedataset_student_full_X.R`、`truedataset_student_excl_X_.R`，以及 `wage/truedataset_wages_full_X_fixedN.R`、`truedataset_wages_excl_X_fixedN.R`。

Student Full/Excl分别使用原脚本的全部/排除分组信息的特征清单；响应为G3，不使用G1/G2。两门课的1044条记录通过历史student-merge.R指定的13个属性形成662个身份签名，同签名记录始终在同一分区。该规则不向算法提供分组，仅用于防止跨课程记录泄漏；它是属性匹配，不是核验后的个人ID。

NLSY79沿用原始负值缺失处理、删除R8378703及完整案例清洗，得到5247条记录；响应为R8495200。Full/Excl分别沿用原特征排除清单。受访者ID只用于划分，不能进入预测特征。

## 方法与评价

- Student 和北京的正文表使用 RI、RS、Full/Excl 和 50 次配对随机划分。
- 外层约70%开发、30%最终测试；开发集再拆成约50%有标签训练、20%独立无标签参考。比例按身份单元计算，Student的记录数比例可略有不同。
- X的类别编码、缩放参数和y的标准化只从训练集估计。参考和测试不参与这些参数的估计；参考y不提供给拟合。所有方法使用同一训练与测试样本，不因某种分组过滤测试样本。
- PGR保留内部标识CIRG，直接使用现有RASC、模型匹配IMSPE和SGA预测；lambda=0、tau=5(p+1)、SA=25、ML-EM=300次上限及容差1e-5。固定lambda=0并不检验残差增广的收益。
- LM、KM、OBS、CPF、BLM与本轮主实验一致；OBS/CPF/BLM共享20个随机虚拟组。BLMX不适用，CPFO不运行。
- Student 的 Full 面板纳入 10 种原始变量分组。北京的 Full 面板纳入站点、季节及站点与季节的交叉分组。
- 旧D-opt/I-opt列不能直接套用到新版IMSPE结果。当前LMM实现是ML-EM，不应标成精确REML。
- 输出标准化响应尺度及原始响应单位的MSPE；误差棒口径为均值±1.96倍MCSE，并保留所有有限但未收敛结果。失败记录保留，不以0替代；未达到50个有效重复的条目不能作为完整50次实验使用。
- 论文每个面板使用其中所有列和行共同成功的重复。Student Full 为 49 次，其余面板为 50 次。

## 输出

- `student_results.rds`、`wage_results.rds`：每个数据集一次最终保存，包括每次划分、训练预处理参数、逐次得分、收敛情况和来源指纹。
- `real_data_results.rds`：合并后的最终得分与统计量，不复制模型对象。
- `real_data_raw.csv`、`real_data_summary.csv`：逐次结果及均值/MCSE/区间/有效重复数/收敛计数。
- `method_comparison_summary.csv`：六方法共同成功重复上的主对比统计，与论文图片一致。
- `paired_comparisons.csv`：同一共同划分上PGR与各对照的配对MSPE差；负值表示PGR更低。
- `comparison_excluded_splits.csv`：共同对比未使用的划分及失败方法。原始记录不删除。
- `failures.csv`：失败记录。完整度标记要求所有方法都达到计划重复次数。
- [论文表格](output/real_data_RI_tables/tables.tex)：当前 Student 与北京的 RI/RS 表。
- `result/real_data_RI/`、`result/real_data_latest/`：此前的 Student 与 NLSY79 结果，供查阅。

预处理的本地单元验证（不是HPC任务）：`Rscript test_real_data.R`。
该检查验证分区独立性、训练集预处理和样本量下界，不保证每次随机聚类都能找到满足最小簇大小的分区；后者的失败必须按实际结果披露。

## 存档 RI 结果与聚类诊断

NLSY79的RI六方法均有50个可评分结果。Student Full第48次，当前内部seed=1的KM和PGR未找到满足tau=155的聚类，RI中产生2条失败记录；这2个条目各有49个有效结果，不能写成50次全部成功。主图和论文表为保证同一比较集合，均只使用其余49次。未降低tau、替换成LM或重新抽取划分。

2026-11-11诊断更正：此前“30个固定初值”的说法不准确，旧检查只改了R外层种子，没有传入ClusterR内部seed，所以实际重复使用内部默认seed=1。真正改变内部seed后，seed=2可以得到满足tau的两簇分区；因此这是当前初始化下的失败，不能表述成数据不存在可行分组。诊断不使用测试响应，且没有据此重跑或替换任何已报告的预测分数。

PGR的199条可评分RI拟合均已达到收敛容差；其他方法仍可能有未收敛记录，详见逐项统计。
