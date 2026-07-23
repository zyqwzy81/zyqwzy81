# ============================================================================
# 构建 SEQTaRget v25 分析数据集
# 目的: 将 35 候选变量 (mice_input_35vars.parquet) 合并到序贯试验结构
# (sequential_trials.parquet), 生成完整的分析输入文件
# ============================================================================

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
import json
from datetime import datetime

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
OUTPUT_DIR = PROJECT_ROOT / "TTE" / "v25" / "output"

# ============================================================
# 1. 加载数据
# ============================================================
print("=== 步骤1: 加载数据 ===")

seq_trials = pd.read_parquet(OUTPUT_DIR / "sequential_trials.parquet")
mice_vars = pd.read_parquet(OUTPUT_DIR / "mice_input_35vars.parquet")
augmented = pd.read_parquet(OUTPUT_DIR / "sequential_trials_augmented.parquet")

print(f"  序贯试验数据: {seq_trials.shape[0]} 行, {seq_trials.shape[1]} 列")
print(f"  变量数据: {mice_vars.shape[0]} 行, {mice_vars.shape[1]} 列")
print(f"  增强数据: {augmented.shape[0]} 行, {augmented.shape[1]} 列")

# ============================================================
# 2. 合并数据
# ============================================================
print("\n=== 步骤2: 合并数据 ===")

# mice_input_35vars 可能有重复 stay_id (每个 trial period 一行)
# 先去重, 保留每个 stay_id 的第一行
dup_count = mice_vars['stay_id'].duplicated().sum()
if dup_count > 0:
    print(f"  ⚠️ 发现 {dup_count} 个重复 stay_id, 正在去重...")
    mice_vars = mice_vars.drop_duplicates(subset=['stay_id'], keep='first')
    print(f"  去重后: {mice_vars.shape[0]} 行")

# 移除与 seq_trials 重叠的列 (保留 stay_id)
overlap_cols = set(mice_vars.columns) & set(seq_trials.columns)
cols_to_drop = overlap_cols - {"stay_id"}
mice_for_merge = mice_vars.drop(columns=list(cols_to_drop))
print(f"  从 mice_vars 移除重叠列: {cols_to_drop}")

# 左连接: 保留序贯试验的所有行
merged = seq_trials.merge(mice_for_merge, on="stay_id", how="left")

print(f"  合并后: {merged.shape[0]} 行, {merged.shape[1]} 列")
print(f"  合并后 stay_id 唯一数: {merged['stay_id'].nunique()}")

# 🆕 v25.8: 合并增强变量 (S/F等)
# 增强数据集每个stay_id有多个trial period行，先去重
augmented_dedup = augmented.drop_duplicates(subset=['stay_id'], keep='first')
aug_cols = set(augmented_dedup.columns) - set(merged.columns)
aug_cols.add("stay_id")
aug_for_merge = augmented_dedup[list(aug_cols)]
merged = merged.merge(aug_for_merge, on="stay_id", how="left")
print(f"  合并增强后: {merged.shape[0]} 行, {merged.shape[1]} 列 (新增 {len(aug_cols)-1} 列)")

# ============================================================
# 3. 变量分类 (time_varying vs fixed)
# ============================================================
print("\n=== 步骤3: 变量分类 ===")

# 时变变量: 仅 baseline_phosphate (每个 trial period 可能不同)
TIME_VARYING = ["baseline_phosphate"]

# 基线变量 (fixed): 25 变量标准方案中的固定变量
FIXED_VARS = [
    "age", "gender_female",
    "baseline_sofa_t0", "admission_sofa", "apsiii", "sepsis3", "egfr",
    "invasive_vent_at_t0", "vasoactive_at_t0", "albumin_min",
    "charlson_comorbidity_index", "liver_disease",
    "insulin_at_t0", "diuretic_at_t0", "steroid_at_t0",
    "nahco3_infusion_at_t0", "rrt_at_t0",
    "bmi", "gcs",
    "diabetes", "heart_failure", "hypertension", "renal_disease",
    "hypokalemia", "hemoglobin_min", "lactate_t0", "bicarbonate_worst",
    # 🆕 氧合指标 (S/F替代P/F, v11四源整合)
    "pao2fio2ratio", "sf_ratio", "sf_estimated_pf", "sf_fio2_defaulted",
    # 🆕 V25.9 营养代谢指标 (来源: baseline_dim7_nutrition → augmented)
    "refeeding_risk", "nutrition_support_type", "en_kcal_24h", "tpn_kcal_24h",
    "en_start_hours_from_t0", "gnri", "pni",
]

# 验证变量存在
missing_tv = [v for v in TIME_VARYING if v not in merged.columns]
missing_fixed = [v for v in FIXED_VARS if v not in merged.columns]
if missing_tv:
    print(f"  ⚠️ 时变变量缺失: {missing_tv}")
if missing_fixed:
    print(f"  ⚠️ 基线变量缺失: {missing_fixed}")

print(f"  时变变量: {len(TIME_VARYING)} 个")
print(f"  基线变量: {len(FIXED_VARS)} 个 (匹配 {len([v for v in FIXED_VARS if v in merged.columns])} 个)")

# ============================================================
# 4. 保存分析数据集
# ============================================================
print("\n=== 步骤4: 保存数据集 ===")

# 保存完整分析数据集
analysis_path = OUTPUT_DIR / "sequential_trials_v25_analysis.parquet"
merged.to_parquet(analysis_path, index=False)
print(f"  已保存: {analysis_path}")

# 保存 CSV 备份 (R 回退兼容)
csv_path = OUTPUT_DIR / "sequential_trials_v25_analysis.csv"
merged.to_csv(csv_path, index=False)
print(f"  已保存: {csv_path}")

# ============================================================
# 5. 保存变量列表 JSON
# ============================================================
print("\n=== 步骤5: 保存变量列表 ===")

# 标准方案 25 变量
var_list_standard = {
    "time_varying": TIME_VARYING,
    "fixed": [v for v in FIXED_VARS if v in merged.columns]
}

var_list_path = OUTPUT_DIR / "var_list_standard_25.json"
with open(var_list_path, "w", encoding="utf-8") as f:
    json.dump(var_list_standard, f, indent=2, ensure_ascii=False)
print(f"  已保存: {var_list_path}")

# 紧凑方案 15 变量 (V25.9: 扩展为18变量，纳入关键营养指标)
COMPACT_VARS = [
    "age", "gender_female",
    "baseline_sofa_t0", "apsiii", "sepsis3", "egfr",
    "invasive_vent_at_t0", "vasoactive_at_t0", "albumin_min",
    "charlson_comorbidity_index", "liver_disease",
    "insulin_at_t0", "diuretic_at_t0", "nahco3_infusion_at_t0",
    # 🆕 V25.9 营养核心变量(3个)
    "refeeding_risk", "nutrition_support_type", "en_kcal_24h",
]

var_list_compact = {
    "time_varying": TIME_VARYING,
    "fixed": [v for v in COMPACT_VARS if v in merged.columns]
}

var_list_compact_path = OUTPUT_DIR / "var_list_compact_15.json"
with open(var_list_compact_path, "w", encoding="utf-8") as f:
    json.dump(var_list_compact, f, indent=2, ensure_ascii=False)
print(f"  已保存: {var_list_compact_path}")

# 扩展方案 35 变量 (V25.9: 扩展为42+变量，纳入完整营养指标)
EXTENDED_EXTRA = [
    "malignant_cancer", "creatinine_t0", "glucose_max", "ph_worst",
    "wbc_max", "platelets_min", "hypocalcemia", "ppi_use",
    # 🆕 V25.9 营养扩展变量(6个，紧凑已含3个核心变量)
    "tpn_kcal_24h", "en_start_hours_from_t0", "gnri", "pni",
]

var_list_extended = {
    "time_varying": TIME_VARYING,
    "fixed": [v for v in (FIXED_VARS + EXTENDED_EXTRA) if v in merged.columns]
}

var_list_extended_path = OUTPUT_DIR / "var_list_extended_35.json"
with open(var_list_extended_path, "w", encoding="utf-8") as f:
    json.dump(var_list_extended, f, indent=2, ensure_ascii=False)
print(f"  已保存: {var_list_extended_path}")

# ============================================================
# 6. 保存元数据
# ============================================================
print("\n=== 步骤6: 保存元数据 ===")

metadata = {
    "timestamp": datetime.now().isoformat(),
    "n_rows": int(merged.shape[0]),
    "n_columns": int(merged.shape[1]),
    "n_unique_stay_ids": int(merged['stay_id'].nunique()),
    "time_varying_count": len(TIME_VARYING),
    "fixed_count": len([v for v in FIXED_VARS if v in merged.columns]),
    "standard_25": len([v for v in FIXED_VARS if v in merged.columns]),
    "compact_15": len([v for v in COMPACT_VARS if v in merged.columns]),
    "extended_35": len([v for v in (FIXED_VARS + EXTENDED_EXTRA) if v in merged.columns]),
    "source_files": {
        "sequential_trials": "sequential_trials.parquet",
        "variables": "mice_input_35vars.parquet"
    }
}

metadata_path = OUTPUT_DIR / "sequential_trials_v25_analysis_metadata.json"
with open(metadata_path, "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2, ensure_ascii=False)
print(f"  已保存: {metadata_path}")

print("\n=== 完成! ===")
print(f"输出目录: {OUTPUT_DIR}")
print(f"分析数据集: sequential_trials_v25_analysis.parquet")
print(f"变量列表: var_list_standard_25.json, var_list_compact_15.json, var_list_extended_35.json")