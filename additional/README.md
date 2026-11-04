# 其余场景的补充图片

本目录包含本轮实验中未纳入论文CASE1-7的其余9个场景。每个场景一张四面板图，格式与论文图1一致。

## 编号与文件

文件名中的 `CODE_CASE` 指代码的原始CASE编号，不是论文CASE编号，也不重新编号为论文CASE8-16。

| 代码场景 | PDF文件（在pdf/下） |
| --- | --- |
| CASE1 | CODE_CASE01_mean_95MC.pdf |
| CASE2 | CODE_CASE02_mean_95MC.pdf |
| CASE4 | CODE_CASE04_mean_95MC.pdf |
| CASE5 | CODE_CASE05_mean_95MC.pdf |
| CASE6 | CODE_CASE06_mean_95MC.pdf |
| CASE7 | CODE_CASE07_mean_95MC.pdf |
| CASE8 | CODE_CASE08_mean_95MC.pdf |
| CASE9 | CODE_CASE09_mean_95MC.pdf |
| CASE10 | CODE_CASE10_mean_95MC.pdf |

论文CASE1-7分别对应代码CASE11、12、13、14、15、3、16；这些图仍保留在 `output/` 中。本目录的9个场景加上论文的7个场景，恰好覆盖本轮全部16个场景，没有遗漏或重复。

## 格式与统计口径

- `pdf/`：9张矢量PDF，供论文排版使用。
- `figures/`：9张同名360 dpi PNG，供预览或位图系统使用。
- 上排为RI、下排为RS；左列Var(a)=0.50、右列Var(a)=2.25。
- RI的Var(b)=0，RS的Var(b)=0.10；该固定设置在配套图注中说明。
- CIRG、KM、BLM、CPF、OBS、LM的颜色、点形、线型和坐标名称与正文单CASE图一致。
- 每个点是50次重复的平均测试MSPE；上下界为均值加减1.96倍Monte Carlo标准误（MCSE=样本标准差/sqrt(50)）。这是均值的近似95% Monte Carlo区间，不是单次实验的波动范围或预测区间。
- 每张图72个点，均有上下界；四个面板分别使用自身的纵轴刻度，无断轴。方法略作横向错位以避免重叠，不代表R取值不同。
- 保留全部有限预测，不删除未收敛记录；BLMX在当前设定下不适用，不绘制。
- 图内不加CASE大标题和解释性段落；编号由文件名和论文图注标识。

## 复现与图注

从项目根目录执行：

```bash
Rscript output/source/export_paper_figures.R result/results_main/main_results.rds additional additional
```

命令仅读取本轮最终结果，生成本目录中的派生图片，不拟合模型、不提交任务、不修改原始结果或已生成的七张论文图。

`figure_captions.tex` 提供九张图的LaTeX插图代码和英文图注。`source/additional_figure_statistics.rds` 保存648个绘图点的统计量、原始代码CASE编号、收敛计数、结果文件指纹；`source/figure_manifest.rds` 保存图像清单。
