# -*- coding: utf-8 -*-
"""
增强数据集扩展脚本 (V25.4)
============================
1. 从 mimiciv_icu.icustays 提取 first_careunit (ICU类型)
2. 合并 supplementary_variables.parquet (补充变量)
3. 输出: sequential_trials_augmented_full.parquet (全变量增强数据集)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2

# 路径
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
OUTPUT_DIR = PROJECT_ROOT / "TTE" / "v26" / "output"
AUGMENTED_PARQUET = OUTPUT_DIR / "sequential_trials_augmented.parquet"
SUPP_PARQUET = OUTPUT_DIR / "supplementary_variables.parquet"
OUTPUT_PARQUET = OUTPUT_DIR / "sequential_trials_augmented_full.parquet"
OUTPUT_CSV = OUTPUT_DIR / "sequential_trials_augmented_full.csv"

print("=" * 60)
print("  增强数据集扩展 V25.4")
print("=" * 60)

# ============================================================
# 1. 加载现有增强数据集
# ============================================================
print("\n[1/5] 加载现有增强数据集...")
df_aug = pd.read_parquet(AUGMENTED_PARQUET)
print(f"  现有列数: {len(df_aug.columns)}, 行数: {len(df_aug):,}")
print(f"  现有列: {list(df_aug.columns)}")

# ============================================================
# 2. 提取ICU类型
# ============================================================
print("\n[2/5] 提取ICU类型 (first_careunit)...")
conn = psycopg2.connect(host='localhost', port=5432, dbname='mimiciv3', user='postgres', password='1')

icu_query = """
    SELECT stay_id, first_careunit AS icu_type
    FROM mimiciv_icu.icustays
    WHERE stay_id IN (SELECT DISTINCT stay_id FROM public.study_cohort)
"""
df_icu = pd.read_sql(icu_query, conn)
conn.close()
print(f"  ICU类型记录: {len(df_icu):,} unique stays")
print(f"  ICU类型分布:\n{df_icu['icu_type'].value_counts().to_string()}")

# 合并
df_aug = df_aug.merge(df_icu, on='stay_id', how='left')
print(f"  合并后缺失: {df_aug['icu_type'].isna().sum()}")

# ============================================================
# 3. 合并补充变量
# ============================================================
print("\n[3/5] 合并补充变量...")
df_supp = pd.read_parquet(SUPP_PARQUET)
print(f"  补充变量: {len(df_supp.columns)} 列, {len(df_supp):,} 行")

# 选择要合并的列（排除重复的标识列）
# 使用 (stay_id, time_zero) 复合键进行1:1合并
supp_merge_cols = [c for c in df_supp.columns 
                   if c not in ['subject_id', 'hadm_id']]

df_supp_subset = df_supp[supp_merge_cols].copy()

# 合并 (stay_id, time_zero) 确保1:1匹配
df_aug = df_aug.merge(df_supp_subset, on=['stay_id', 'time_zero'], how='left')
print(f"  合并后总列数: {len(df_aug.columns)}")
print(f"  合并后总行数: {len(df_aug):,}")

# ============================================================
# 4. 数据质量检查
# ============================================================
print("\n[4/5] 数据质量检查...")

# 4.1 新增变量缺失率
new_vars = ['icu_type', 'magnesium_t0', 'phosphate_nadir', 'alcohol_abuse', 
            'peep_t0', 'fio2_t0', 'ventilation_duration_hours', 'crp_t0',
            'tidal_volume_observed_t0', 'plateau_pressure_t0', 'ventilator_mode_t0',
            'sofa_t0', 'cardiovascular_event_2comp', 'cardiovascular_event_3comp',
            'cardiovascular_event_comprehensive', 'new_arrhythmia_rhythm',
            'vasopressor_escalation', 'myocardial_injury', 'ecg_abnormal']

print(f"\n  {'变量':<35s} {'非空数':>8s} {'缺失%':>7s} {'均值':>10s} {'备注'}")
print(f"  {'-'*80}")
for col in new_vars:
    if col not in df_aug.columns:
        print(f"  {col:<35s} {'N/A':>8s} {'N/A':>7s} {'N/A':>10s}  列不存在!")
        continue
    s = df_aug[col]
    n_nonnull = s.notna().sum()
    miss_pct = round(100.0 * s.isna().sum() / len(s), 1)
    if s.dtype in ['float64', 'int64']:
        mean_val = f"{s.mean():.2f}" if n_nonnull > 0 else 'N/A'
    else:
        mean_val = 'cat'
    print(f"  {col:<35s} {n_nonnull:>8,} {miss_pct:>6.1f}% {mean_val:>10s}")

# 4.2 行数一致性
assert len(df_aug) == 32330, f"行数不一致: {len(df_aug)} != 32330"
print(f"\n  行数一致性: PASS ({len(df_aug):,})")

# ============================================================
# 5. 导出
# ============================================================
print(f"\n[5/5] 导出...")

# 确保时间戳列
for col in df_aug.columns:
    if 'time' in col.lower() and df_aug[col].dtype == 'object':
        try:
            df_aug[col] = pd.to_datetime(df_aug[col])
        except:
            pass

df_aug.to_parquet(OUTPUT_PARQUET, index=False)
print(f"  Parquet: {OUTPUT_PARQUET}")
print(f"  维度: {df_aug.shape[0]:,} 行 × {df_aug.shape[1]} 列")

df_aug.to_csv(OUTPUT_CSV, index=False)
print(f"  CSV: {OUTPUT_CSV}")

# 列名列表
print(f"\n  全部列名 ({len(df_aug.columns)}列):")
for i, c in enumerate(df_aug.columns):
    print(f"    {i+1:3d}. {c}")

print("\n=== 完成 ===")