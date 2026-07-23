# SEQTaRget 兼容分批引擎 V2

## 设计目标

忠实复现 SEQTaRget 的 CCW 扩展逻辑（SEQexpand），通过分批处理实现内存可控，逐期验证结果一致性。

## 与 V1（seqtarget_batch_ccw）的核心差异

| 维度 | V1 (batch_ccw) | V2 (本引擎) |
|:---|:---|:---|
| **扩展深度** | 1 patient-period = 1 clone | **日粒度展开**: 每条记录展开为 day 0→floor(surv_time) |
| **PS 模型** | **Per-period**: 每期独立拟合 PS | **Pooled**: 全量数据一次性拟合 PS |
| **结局模型** | Per-period Cox + Meta 合并 | **Pooled logistic** + followup² + trial² + cluster SE |
| **权重策略** | Per-period 独立 IPTW | **Pre-expansion IPTW** (weight.preexpansion=TRUE) |
| **验证机制** | 最终结果对比 | **累加式收敛验证**: 逐期检查 HR 是否收敛 |

## 解决的核心问题

V1 与 SEQTaRget 之间数据不一致的 **3 个根因**:

1. **扩展深度不匹配**（最大影响）: V1 每条 patient-period 生 1 条 clone，SEQTaRget 生成 floor(surv_time)+1 条日记录。展开倍率约 10x，信息量差异巨大。
2. **PS 模型策略不匹配**: V1 每期独立拟合 PS（小样本期不稳定），SEQTaRget 使用 pooled PS（weight.preexpansion=TRUE），利用了全量数据信息。
3. **结局模型不匹配**: V1 用 per-period Cox + Meta 合并（仅用 4-6 个 period），SEQTaRget 用 pooled logistic + 时间多项式（利用全部日粒度数据）。

## 分批策略

### 批次分割规则

按 `trial_period` 自然分割，每期一个独立批次:

| Period | 原始记录 | 展开后预估行数 | 事件数 |
|:---:|:---:|:---:|:---:|
| 1 | 22,685 | ~218,910 | 3,628 |
| 2 | 7,262 | ~76,614 | 1,409 |
| 3 | 1,800 | ~20,502 | 421 |
| 4 | 451 | ~5,466 | 114 |
| 5 | 103 | ~1,191 | 30 |
| 6+ | < 30 或事件<5 | 自动跳过 | - |

单批次最大约 22 万行（period 1），**不需要进一步细分**。

### 批次处理顺序

```
Period 1 → Period 2 → Period 3 → ... → Period K
(累加)       (累加)       (累加)
  ↓           ↓           ↓
HR_1        HR_1-2      HR_1-3    ... → HR_full
```

每期读取展开数据后**累加到 cumulative 数据集**，拟合 cumulative pooled logistic。这样:
- 早期（period 1-3）即可获得稳定的 HR 估计（>5000 事件）
- 后期 period 的加入验证 HR 是否收敛
- 若 HR 在 period 3 后已稳定，说明后期数据对结果影响有限

### 批次内存控制

- Period 1 展开约 22 万行，约 30MB Parquet，**单次加载无内存风险**
- Cumulative 数据最多约 32 万行（全量展开），约 50MB，单次加载可行
- 若需进一步拆分（如 period 1 内再分），可在 `02_batched_seqexpand.R` 中调整

## 管线结构

```
seqtarget_batch_v2/
├── config_v2.R              # 统一配置（路径、变量、参数）
├── 01_pooled_ps.R            # Phase 1: Pooled PS + IPTW
├── 02_batched_seqexpand.R    # Phase 2: 分批日粒度展开
├── 03_batched_logistic.R     # Phase 3: 累加式 pooled logistic
├── 04_quality_check.R        # Phase 4: 质量检查 + 对比
├── run_all.R                 # 完整管线入口
└── README.md                 # 本文档
```

## 运行方法

```bash
cd TTE/v26/code/r/seqtarget_batch_v2/
Rscript run_all.R
```

或逐阶段运行（支持断点续跑）:

```bash
Rscript 01_pooled_ps.R        # Phase 1: ~30秒
Rscript 02_batched_seqexpand.R  # Phase 2: ~2分钟
Rscript 03_batched_logistic.R   # Phase 3: ~5分钟
Rscript 04_quality_check.R      # Phase 4: ~1分钟
```

## 输出文件

| 文件 | 路径 | 说明 |
|:---|:---|:---|
| **Phase 1** | | |
| `pooled_ps_weights_summary.csv` | results/ | 每期权重摘要 |
| `pooled_ps_diagnostics.csv` | results/ | PS 模型诊断（AIC、ESS） |
| `01_pooled_ps_metadata.json` | results/ | Phase 1 元数据 |
| | weighted/ | 按 period 分区的权重 Parquet |
| **Phase 2** | | |
| `seqexpand_summary.csv` | results/ | 展开摘要（原始/展开行数、倍率） |
| `seqexpand_event_by_day.csv` | results/ | 事件按日分布 |
| `02_seqexpand_metadata.json` | results/ | Phase 2 元数据 |
| | expanded/ | 按 period 分区的展开后 Parquet |
| **Phase 3** | | |
| `cumulative_hr_by_period.csv` | results/ | **核心输出**: 逐期累加 HR |
| `final_hr_summary.csv` | results/ | 全量数据最终 HR |
| `Fig_cumulative_hr_convergence.pdf` | results/ | HR 收敛曲线图 |
| `03_logistic_metadata.json` | results/ | Phase 3 元数据 |
| **Phase 4** | | |
| `quality_report.json` | results/ | 综合质量报告 |
| `v2_vs_seqtarget_comparison.csv` | results/ | 与 SEQTaRget 对比（如可用） |

## 关键参数（config_v2.R）

| 参数 | 默认值 | 说明 |
|:---|:---|:---|
| `PS_COVARIATES` | DAG12 核心变量 | PS 模型协变量集 |
| `IPTW_TRUNC_UPPER` | 0.99 | 上侧截断分位数 |
| `IPTW_TRUNC_LOWER` | 0.01 | 下侧截断分位数 |
| `SEQEXPAND_MAX_DAY` | 28 | 最大展开天数 |
| `PERIOD_FILTER_MIN_EVENTS` | 5 | 最小事件数阈值 |
| `PERIOD_FILTER_MIN_N` | 30 | 最小样本量阈值 |

## 验证机制

### 1. 逐期累加收敛验证（核心）
```
Period 1-1: HR=0.89 (CI 0.83-0.96)
Period 1-2: HR=0.88 (CI 0.84-0.92)  ← 已稳定
Period 1-3: HR=0.87 (CI 0.85-0.90)  ← 收敛判定: 后50%波动 < 10%
```

### 2. 格式完整性校验
- 每期权重文件: 必需列存在、值域合理
- 每期展开文件: day 范围 ≤ 28, y_daily 编码正确
- 结果文件: HR > 0, SE > 0, CI 合理

### 3. 跨版本对比
- V2 vs V1: 理解差异来源（扩展深度 + PS 策略 + 结局模型）
- V2 vs SEQTaRget: 小样本子集上直接对比（需 SEQTaRget 包）

### 4. 结构化记录
- 所有 Phase 输出包含 JSON 元数据文件
- 质量报告包含 PASS/FAIL 判定和详细说明
- 关键参数在 config_v2.R 中集中管理
