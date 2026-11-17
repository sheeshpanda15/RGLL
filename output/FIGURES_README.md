# 论文图片与PDF报告

本目录基于本轮全部50次重复生成，没有重新拟合模型，没有删除未收敛但有有限预测的记录，也没有创建压缩包。

## 推荐使用

- `pdf/main_CASE1-7_RI_mean_95MC.pdf`：RI组合图，七行对应论文CASE1-7，两列为Var(a)=0.50和2.25。
- `pdf/main_CASE1-7_RS_mean_95MC.pdf`：RS组合图，相同布局，Var(b)=0.10。
- `pdf/CASE01_mean_95MC.pdf` 至 `CASE07_mean_95MC.pdf`：按原论文图1样式制作的七张单CASE四面板图。上排RI、下排RS；左列Var(a)=0.50、右列2.25。每个面板使用两层灰色标题，底部共用六方法图例，保留原论文的坐标名称；不在图内重复加入CASE大标题。CASE编号放在文件名和论文图注中。RS的Var(b)=0.10，RI的Var(b)=0，在图注中注明。
- `pdf/PAPER_COMPARISON_2026-10-28.pdf`：完整中文论文对照报告，共6页，正文可检索，附PDF书签。
- `figures/`：全部图的360 dpi PNG版本。插入论文时优先使用矢量PDF，PNG用于预览或要求位图的系统。
- `paper_figure_captions.tex`：可用于英文论文的图注与插图代码。
- `pdf/real_data_RI_methods.pdf`：真实数据RI六方法四面板对比；`real_data_RI_student.pdf`和`real_data_RI_NLSY79.pdf`为分数据集版本。同名PNG在figures/，英文图注在real_data_RI_captions.tex。Student Full统一用共同成功的49次，其余设置50次；结果与合成实验分开，不覆盖RS历史结果。
- 项目根目录的 `additional/`：其余9个代码场景（CASE1、2、4-10）的同格式单CASE图，包括矢量PDF、360 dpi PNG和独立英文图注。其文件名使用 `CODE_CASE`，不与论文CASE编号混用。

主图与单CASE图是两种组织方式，内容相同。论文可采用两张组合图，或采用单CASE图；不必把两套都重复放入正文。组合图建议全页宽使用，不要缩成单栏。

## 上下界分别表示什么

### 主图：均值的精度

每点是50次独立重复的平均测试MSPE。上下界为：

```text
MCSE = sd(MSPE over 50 repetitions) / sqrt(50)
lower = mean(MSPE) - 1.96 * MCSE
upper = mean(MSPE) + 1.96 * MCSE
```

这是均值的近似95% Monte Carlo区间，不是未来观测的预测区间，也不是“95%的单次实验都在此范围内”。各方法在相同重复之间配对，因此不能仅凭两条均值区间是否重叠来判断配对差异是否显著。

### 补充图：重复之间的波动

`pdf/robustness_CASE1-7_RI_median_Q05-Q95.pdf` 和RS版本使用：

```text
point = median(MSPE over 50 repetitions)
lower = empirical 5th percentile
upper = empirical 95th percentile
```

分位数使用R的默认type=7插值。该范围描述样本中的中间90%重复波动，比均值的置信区间更适合辅助观察稳健性；它不是置信区间，也不能保证未来90%实验一定落在区间内。

两个版本的每个点都有上下界。方法之间略作横向错位，仅为避免标记和误差棒重叠，不代表使用了不同的R值。

## 保留的设置和限制

- 论文CASE1-7严格映射到代码CASE11、12、13、14、15、3、16。
- 保留PGR、KM、BLM、CPF、OBS、LM六种实际运行的方法；BLMX为不适用，不被当作零误差画出。PGR是结果数据中CIRG方法的展示名称，仅修改图例，不改变代码内部标识、数值、颜色、点形或线型。
- 本轮lambda=0；这些图不能作为非零残差增强收益的证据。
- 在这七场景内，PGR最终模型RI有1/2100条、RS有8/2100条未收敛；其余方法也有未收敛记录。所有有限预测均按原结果保留，不能把图注写成“所有拟合已收敛”。
- RS中的BLM退化为K=1，预测几乎与LM重合；横向错位保留两者标记，不制造数值差异。
- 组合图中不同CASE使用不同纵轴范围，同一CASE的两列Var(a)共享纵轴；单CASE图的四个面板各自使用纵轴范围。图中均标出真实刻度，没有断轴。
- 图片内部仅保留必要的模型、Var设置、坐标、图例和数据（组合图另保留CASE行标签），不放结论性说明。区间定义和拟合限制应写进图注或正文。

## 图表统计材料

- `source/paper_figure_statistics.rds`：504个参数/方法点的均值、MCSE、上下界、分位数、收敛计数与来源指纹。
- `source/figure_manifest.rds`：11张图的文件清单、点数和区间类型。
