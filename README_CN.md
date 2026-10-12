# ICU低磷血症患者静脉补磷治疗：目标试验模拟研究

本仓库包含一项目标试验模拟研究的分析代码，评估静脉补磷对ICU低磷血症患者28天死亡率的影响，数据来源于MIMIC-IV数据库（2008–2022年）。

## 研究概述

- **设计**：序贯目标试验模拟（克隆-删失-加权框架）
- **人群**：ICU轻中度低磷血症患者（血磷 0.32–0.80 mmol/L）
- **暴露**：低磷血症检测后24小时内静脉补磷
- **主要结局**：28天全因死亡率
- **分析方法**：通过SEQTaRget R包进行ITT分析，辅以PP和CCW敏感性分析
- **主要发现**：未观察到显著的总体死亡率获益（HR=0.982；95% CI, 0.897–1.076；p=0.695）。基线血磷 <0.50 mmol/L 层的探索性结果（HR=0.765；交互 p=0.0097）因单位依赖性与未做多重校正，降级为假设生成；stay 级敏感性分析不显著（HR=0.746；p=0.117）。

## 仓库结构

```
.
├── README.md                              # 英文说明
├── README_CN.md                           # 本文件（中文说明）
├── data_extraction/
│   ├── python/                            # Python 数据提取脚本（PostgreSQL → Parquet）
│   │   ├── extract_event_times.py         #   KDIGO 分期 + PRBC 输注事件时间（次要结局）
│   │   ├── build_sequential_trials.py     #   序贯试验起点构建（全分析的上游）
│   │   ├── augment_full_dataset.py        #   试验级扩展数据集（含基线协变量）
│   │   └── …                              #   其他预处理脚本
│   └── sql/pure_control/                  # 纯对照 post-hoc 分析的 SQL（M0_40–M0_48）
├── analysis/
│   └── r/
│       ├── main_multi.R                   # 主 ITT（多次入组序贯目标试验）
│       ├── piecewise_multi.R              # 分段风险比 + PP（带 IPCW 删失权重）
│       ├── subgroup_multi.R               # 亚组异质性检验（多次入组）
│       ├── secondary_multi.R              # 次要结局（二分类 + 连续，MI-bootstrap）
│       ├── death90_multi.R                # 90 天死亡率（patients.dod 核验）
│       ├── fig2_fullkm.R                  # 完整随访加权 KM（不在出院处截断）
│       ├── binary_continuous.R            # stay 级敏感性（二分类 + 连续次要）
│       ├── daily_threshold.R              # 日阈值扫描（stay 级）
│       ├── stay_common.R / multi_common.R / common.R   # 公共模块
│       ├── generate_all_figures.R         # 发表级图表（读取落盘结果表）
│       └── pure_control/                  # IV vs 无补充 post-hoc 机制分析（PC01–PC12b）
└── .gitignore
```

## 分析代码

仓库内全部 R 脚本均实现**期刊修回版**分析（2026-10-10 定版），采用重建的**多次入组（trial-period 级）序贯目标试验模拟**，自然日历日级时间线，为稿件的确定性分析。公共模块：`common.R`、`multi_common.R`、`stay_common.R`。源码中仍嵌入绝对路径（`/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/...`）指向作者工作站的本地位置——仓库内文件名已去掉版本号前缀，但本地路径前缀保留以便在工作站上继续运行。

**数据库口令**：所有脚本均不内嵌 PostgreSQL 口令。运行连接本地 MIMIC-IV 数据库的脚本前，请通过环境变量 `MIMIC_PG_PASSWORD` 设置（例如 `MIMIC_PG_PASSWORD=xxx Rscript secondary_outcomes.R`）。

**患者级数据（parquet、csv）不在仓库中**——脚本产出的分析表保存在本地；持 PhysioNet 认证的作者可通过 `data_extraction/python/` 脚本从 MIMIC-IV 重建。

## 依赖

### R包
- **SEQTaRget** (≥0.1.0)：序贯目标试验模拟框架
- **survival**：Cox比例风险模型
- **mice** (≥3.16.0)：多重插补
- **glmnet**：LASSO变量筛选
- **ggplot2**, **survminer**, **forestplot**：发表级可视化
- **cobalt**：协变量平衡评估
- **EValue**：E-value敏感性分析
- **arrow**：Parquet数据读取
- **sandwich**, **lmtest**：聚类稳健标准误

### Python包
- **pandas**, **pyarrow**：Parquet数据处理
- **SQLAlchemy**：PostgreSQL数据库连接

## 数据来源

本研究使用 **MIMIC-IV v3.1** 数据库（PhysioNet, https://physionet.org/content/mimiciv/3.1/）。访问需要：
1. PhysioNet认证访问权限
2. 完成CITI数据或样本研究培训
3. 签署数据使用协议

分析数据集通过 `data_extraction/` 中的Python脚本从9个MIMIC-IV表构建。

## 重现分析

### 步骤1：数据提取
```bash
cd data_extraction/python
python build_sequential_trials.py
python build_analysis_dataset.py
```

### 步骤2：运行分析
```bash
cd analysis/r
Rscript misc/run_all_analyses.R
```

## 研究注册

研究方案已在Open Science Framework回顾性注册：
https://osf.io/yc6kb/

## 引用

如使用本代码，请引用：

[待发表后补充引用信息]

## 许可

本代码仅供研究使用。MIMIC-IV数据库受单独的数据使用协议约束。