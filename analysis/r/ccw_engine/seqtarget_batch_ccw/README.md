# SEQTaRget 分批扩展-加权-Cox 引擎

## 模块概述

本模块将 SEQTaRget 的 CCW（Clone-Censor-Weight）扩展步骤与建模步骤解耦，实现**扩展一次落盘、分批计算 PS/IPTW、分批拟合 Cox、最终汇总 HR**的流程。目的是在保留 CCW 方法学严谨性的同时，将内存峰值控制在可接受范围，并可通过小样本验证与 SEQTaRget 官方结果的一致性。

## 文件结构

```
TTE/v26/code/r/seqtarget_batch_ccw/
├── config.R                    # 模块配置（路径、变量、截断、种子、汇总方法）
├── 01_expand_only.R            # 步骤一：读取 parquet 并生成扩展后分区 parquet
├── 02_batch_ps_iptw.R          # 步骤二：按 trial_period 分批拟合 PS、计算 IPTW
├── 03_batch_cox.R              # 步骤三：按 trial_period 分批拟合 IPTW-Cox 并汇总 HR
├── 04_compare_with_seqtarget.R # 步骤四：小样本上与官方 SEQTaRget 对比
├── 05_verify_outputs.R         # 步骤五：验证输出文件格式与字段完整性
└── README.md                   # 本文件
```

## 运行顺序

```r
setwd("TTE/v26/code/r/seqtarget_batch_ccw/")
source("config.R")
source("01_expand_only.R")
source("02_batch_ps_iptw.R")
source("03_batch_cox.R")
source("04_compare_with_seqtarget.R")
source("05_verify_outputs.R")
```

或在终端中：

```bash
cd TTE/v26/code/r/seqtarget_batch_ccw/
Rscript 01_expand_only.R
Rscript 02_batch_ps_iptw.R
Rscript 03_batch_cox.R
Rscript 04_compare_with_seqtarget.R
Rscript 05_verify_outputs.R
```

## 输入输出

### 输入

- `TTE/v26/output/sequential_trials.parquet`：主分析级 parquet（19 核心列）
- 可选 `TTE/v26/output/sequential_trials_augmented.parquet`：增强协变量集

### 输出

| 文件 | 路径 | 说明 |
|:---|:---|:---|
| 扩展后分区 parquet | `TTE/v26/output/seqtarget_batch_ccw/expanded/trial_period=*/part-0.parquet` | 每期一个分区 |
| 加权后分区 parquet | `TTE/v26/output/seqtarget_batch_ccw/weighted_batches/trial_period=*/part-0.parquet` | 每期 PS/IPTW 结果 |
| 扩展摘要 | `TTE/v26/output/seqtarget_batch_ccw/results/expansion_summary.{csv,json}` | 克隆数、事件数、内存 |
| PS/IPTW 摘要 | `TTE/v26/output/seqtarget_batch_ccw/results/batch_ps_iptw_summary.csv` | 每期 PS/IPTW 分布 |
| 按期 HR | `TTE/v26/output/seqtarget_batch_ccw/results/seqtarget_batch_hr_period_specific.csv` | 每期 HR、95% CI |
| 汇总 HR | `TTE/v26/output/seqtarget_batch_ccw/results/seqtarget_batch_hr_summary.csv` | fixed/random/period_specific 汇总 |
| 方法学说明 | `TTE/v26/output/seqtarget_batch_ccw/results/seqtarget_batch_ccw_methodology.md` | 论文补充材料 |
| 对比报告 | `TTE/v26/output/seqtarget_batch_ccw/results/seqtarget_vs_batch_engine_comparison.csv` | 与 SEQTaRget 对比 |
| 对比图 | `TTE/v26/output/seqtarget_batch_ccw/results/FigS_batch_vs_seqtarget_hr.{pdf,png}` | HR 比较图 |
| 验证报告 | `TTE/v26/output/seqtarget_batch_ccw/results/verification_report.json` | 输出文件格式与字段完整性验证结果 |

## 与 SEQTaRget 的主要差异

| 维度 | SEQTaRget `SEQuential()` | 本引擎 |
|:---|:---|:---|
| 扩展策略 | 随访纵向扩展（30-100 倍） | 每条 patient-period 记录作为一个 clone |
| 权重模型 | 时变累积 IPCW（ITT 无 LTFU/Visit 时自动退化为无权重） | 每期独立横截面 IPTW |
| 权重截断 | 单侧 99% | 双侧 1%/99% |
| 时间尺度 | 离散 period | 连续 followup_time（天） |
| 结局模型 | 参数化生存 / pooled logistic | 每期 IPTW-Cox + Meta 合并 |
| 缺失处理 | 完全案例 | 每期独立完整案例 |
| 内存管理 | 全内存 | 扩展落盘 + 分批读取 |
| Bootstrap | 包内并行；小样本可能失败并回退到无 bootstrap | 当前版本暂未实现 Bootstrap（后续可扩展） |

本引擎目前聚焦 **ITT 策略** 与标准 IPTW-Cox，尚未覆盖 PP/As-Treated 的跨期累乘权重。

## 小样本验证说明

`04_compare_with_seqtarget.R` 默认在 n=1,000 与 n=2,000 患者子集上同时运行 SEQTaRget 与本引擎。由于 SEQTaRget 在小样本、低事件数 period 上 bootstrap 可能失败，脚本会自动回退到无 bootstrap 模式以保留点估计；此时 `hr_lower`/`hr_upper` 为 `NA`。详细结果与方法学解释见 `results/seqtarget_batch_ccw_methodology.md`。

## 依赖

- R ≥ 4.1
- `arrow`, `data.table`, `survival`, `jsonlite`
- 可选：`meta`（随机效应合并）、`ggplot2`（对比图）、`SEQTaRget`（对比验证）

## 配置项

关键参数位于 `config.R`：

- `PS_COVARIATES`：PS 模型协变量集
- `IPTW_TRUNC_UPPER` / `IPTW_TRUNC_LOWER`：IPTW 截断分位数
- `RANDOM_SEED`：随机种子
- `POOL_METHOD`：HR 汇总方法（`fixed` / `random` / `period_specific`）
- `COMPARE_SAMPLE_SIZES`：小样本验证患者数
