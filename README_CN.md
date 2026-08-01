# ICU低磷血症患者静脉补磷治疗：目标试验模拟研究

本仓库包含一项目标试验模拟研究的分析代码，评估静脉补磷对ICU低磷血症患者28天死亡率的影响，数据来源于MIMIC-IV数据库（2008–2022年）。

## 研究概述

- **设计**：序贯目标试验模拟（克隆-删失-加权框架）
- **人群**：ICU低磷血症患者（血磷 <0.80 mmol/L）
- **暴露**：低磷血症检测后24小时内静脉补磷
- **主要结局**：28天全因死亡率
- **分析方法**：通过SEQTaRget R包进行ITT分析，辅以PP和CCW敏感性分析
- **主要发现**：未观察到显著的总体死亡率获益（HR=0.987；95% CI, 0.893–1.092；p=0.802），在严重低磷血症患者（血磷 <0.50 mmol/L）中观察到潜在保护效应（HR=0.806；p=0.022）

## 仓库结构

```
.
├── README.md                          # 英文说明
├── README_CN.md                       # 本文件（中文说明）
├── data_extraction/
│   └── python/                        # 数据提取Python脚本（PostgreSQL → Parquet）
├── analysis/
│   └── r/
│       ├── config/                    # 配置文件（协变量、路径、参数）
│       ├── main/                      # 主要ITT分析
│       ├── pp_analysis/               # 符合方案分析
│       ├── sensitivity/               # 敏感性分析（CCA、宽限期、截断、MNAR、LASSO vs DAG）
│       ├── secondary/                 # 次要结局分析
│       ├── subgroup/                  # 亚组分析
│       ├── interaction/               # 交互检验
│       ├── external_validation/       # eICU外部验证
│       ├── ccw_engine/                # 分批CCW引擎（SEQTaRget的内存友好替代）
│       ├── utils/                     # 辅助函数和诊断
│       └── tests/                     # 单元测试
└── docs/                              # 文档
```

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