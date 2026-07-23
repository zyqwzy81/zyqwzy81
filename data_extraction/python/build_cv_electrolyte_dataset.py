"""
build_cv_electrolyte_dataset.py
===============================
将心血管事件和电解质数据合并到 sequential_trials.parquet，
输出 sequential_trials_cv.parquet + CSV + metadata JSON。

合并逻辑：按 stay_id + time_zero LEFT JOIN，
确保输出 32330 行，保留所有 19 核心列。
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# ── 路径配置 ──────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
SQL_DIR = PROJECT_ROOT / "sql" / "result"
MAIN_PATH = PROJECT_ROOT / "TTE" / "v25" / "output" / "sequential_trials.parquet"
OUT_DIR = PROJECT_ROOT / "TTE" / "v25" / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EXPECTED_ROWS = 32330
MAIN_CORE_COLS = [
    "stay_id", "subject_id", "trial_num", "time_zero", "itt_group",
    "baseline_phosphate", "age", "gender", "baseline_sofa", "egfr",
    "sepsis3", "invasive_vent_at_t0", "vasoactive_at_t0", "charlson",
    "admission_type", "death_28day_num", "surv_time",
    "death_time_from_t0", "dischtime_from_t0",
]


# ═══════════════════════════════════════════════════════════
# 1. 加载主数据集
# ═══════════════════════════════════════════════════════════
print("=" * 70)
print("1. Loading main dataset …")
df_main = pd.read_parquet(MAIN_PATH)
print(f"   Rows: {len(df_main)}, Cols: {len(df_main.columns)}")
if len(df_main) != EXPECTED_ROWS:
    raise SystemExit(f"ERROR: Expected {EXPECTED_ROWS} rows, got {len(df_main)}.")
if set(MAIN_CORE_COLS) - set(df_main.columns):
    raise SystemExit(f"ERROR: Missing core columns: {set(MAIN_CORE_COLS) - set(df_main.columns)}")

JOIN_KEYS = ["stay_id", "time_zero"]
# time_zero 必须统一时区（与join的来源datasets对比）
print(f"   itt_group distribution:\n{df_main['itt_group'].value_counts().to_string()}")
print()

# ═══════════════════════════════════════════════════════════
# 2. 加载心血管分表，选择性提取所需列
# ═══════════════════════════════════════════════════════════
print("2. Loading CV sub-tables …")

# ── 2a. Composite ──
cv_comp = pd.read_parquet(SQL_DIR / "M0_14d_composite.parquet", columns=JOIN_KEYS + ["cardiovascular_event_comprehensive"])
# 不做 rename（后面统一更名），先暂存为临时列
cv_comp = cv_comp.copy()
print(f"   cv_composite: {len(cv_comp)} rows")

# ── 2b. Arrhythmia ──
arrhythmia_cols = JOIN_KEYS + [
    "new_arrhythmia_rhythm",
    "new_af_afl_rhythm",
    "new_vt_rhythm",
    "new_vf_rhythm",
    "cardiac_arrest_icd",
    "rhythm_event_hours_from_t0",
]
cv_arrh = pd.read_parquet(SQL_DIR / "M0_14a_arrhythmia.parquet", columns=arrhythmia_cols).copy()

# 派生：vt OR vf
cv_arrh["_cv_vt_vf"] = ((cv_arrh["new_vt_rhythm"].fillna(0).astype(bool)) | (cv_arrh["new_vf_rhythm"].fillna(0).astype(bool))).astype(int)
print(f"   cv_arrhythmia: {len(cv_arrh)} rows")

# ── 2c. Hemodynamic ──
hemo_cols = JOIN_KEYS + [
    "vasopressor_escalation",
    "map_sustained_hypotension",
    "co_decline_30pct",
]
cv_hemo = pd.read_parquet(SQL_DIR / "M0_16_cv_hemodynamic.parquet", columns=hemo_cols).copy()
print(f"   cv_hemodynamic: {len(cv_hemo)} rows")

# ── 2d. Myocardial ──
myo_cols = JOIN_KEYS + [
    "myocardial_injury",
    "troponin_severity",
]
cv_myo = pd.read_parquet(SQL_DIR / "M0_14c_myocardial.parquet", columns=myo_cols).copy()

# troponin_severity: categorical → ordinal int
_ts_map = {"not_measured": 0, "normal": 1, "mild": 2, "moderate": 3, "severe": 4}
cv_myo["_cv_troponin_severity"] = cv_myo["troponin_severity"].map(_ts_map).fillna(0).astype(int)
print(f"   cv_myocardial: {len(cv_myo)} rows")

# ═══════════════════════════════════════════════════════════
# 3. 加载电解质数据
# ═══════════════════════════════════════════════════════════
print("\n3. Loading electrolyte data …")
elyte = pd.read_parquet(SQL_DIR / "M0_05f_dim6_hypo_electrolyte.parquet",
                        columns=JOIN_KEYS + ["hypokalemia", "hypocalcemia", "hypomagnesemia"]).copy()
print(f"   dim6_hypo_electrolyte: {len(elyte)} rows")
print()

# ═══════════════════════════════════════════════════════════
# 4. 四表级联 LEFT JOIN
# ═══════════════════════════════════════════════════════════
print("4. Cascading LEFT JOINs …")

df = df_main.copy()

# ── join composite ──
df = df.merge(cv_comp[JOIN_KEYS + ["cardiovascular_event_comprehensive"]], on=JOIN_KEYS, how="left")
# ── join arrhythmia ──
_arrh_use = cv_arrh[JOIN_KEYS + [
    "new_arrhythmia_rhythm", "new_af_afl_rhythm", "_cv_vt_vf",
    "cardiac_arrest_icd", "rhythm_event_hours_from_t0",
]]
df = df.merge(_arrh_use, on=JOIN_KEYS, how="left")
# ── join hemodynamic ──
df = df.merge(cv_hemo[JOIN_KEYS + ["vasopressor_escalation", "map_sustained_hypotension", "co_decline_30pct"]], on=JOIN_KEYS, how="left")
# ── join myocardial ──
df = df.merge(cv_myo[JOIN_KEYS + ["myocardial_injury", "_cv_troponin_severity"]], on=JOIN_KEYS, how="left")
# ── join electrolyte ──
df = df.merge(elyte[JOIN_KEYS + ["hypokalemia", "hypocalcemia", "hypomagnesemia"]], on=JOIN_KEYS, how="left")

print(f"   Merged rows: {len(df)}")
assert len(df) == EXPECTED_ROWS, f"Row count changed! Expected {EXPECTED_ROWS}, got {len(df)}."
assert set(MAIN_CORE_COLS).issubset(set(df.columns)), "Core columns lost!"
print()

# ═══════════════════════════════════════════════════════════
# 5. 列重命名 & 类型标准化
# ═══════════════════════════════════════════════════════════
print("5. Renaming & typing new columns …")

RENAME_MAP = {
    "cardiovascular_event_comprehensive": "cv_composite",
    "new_arrhythmia_rhythm":              "cv_arrhythmia_new",
    "new_af_afl_rhythm":                  "cv_af_afl",
    "_cv_vt_vf":                          "cv_vt_vf",
    "cardiac_arrest_icd":                 "cv_cardiac_arrest",
    "rhythm_event_hours_from_t0":         "cv_arrhythmia_hours",
    "vasopressor_escalation":             "cv_vasopressor_escalation",
    "map_sustained_hypotension":          "cv_hypotension",
    "co_decline_30pct":                   "cv_co_decline",
    "myocardial_injury":                  "cv_myocardial_injury",
    "_cv_troponin_severity":              "cv_troponin_severity",
    "hypokalemia":                        "electrolyte_hypokalemia",
    "hypocalcemia":                       "electrolyte_hypocalcemia",
    "hypomagnesemia":                     "electrolyte_hypomagnesemia",
}

df = df.rename(columns=RENAME_MAP)

# 类型标准化：确保所有新列为适当类型
NEW_COLS = list(RENAME_MAP.values())

_int_cols = [c for c in NEW_COLS if c != "cv_arrhythmia_hours" and c != "cv_vasopressor_hours"]
for c in _int_cols:
    if c in df.columns:
        df[c] = df[c].fillna(0).astype(int)

# cv_arrhythmia_hours → float
if "cv_arrhythmia_hours" in df.columns:
    df["cv_arrhythmia_hours"] = df["cv_arrhythmia_hours"].astype(float)

# cv_vasopressor_hours：在 hemodynamic 中不存在对应小时列 → 创建为 NaN
if "cv_vasopressor_hours" not in df.columns:
    df["cv_vasopressor_hours"] = np.nan
    df["cv_vasopressor_hours"] = df["cv_vasopressor_hours"].astype(float)
    NEW_COLS.append("cv_vasopressor_hours")

print(f"   New columns ({len(NEW_COLS)}): {NEW_COLS}")
print()

# ═══════════════════════════════════════════════════════════
# 6. 验证
# ═══════════════════════════════════════════════════════════
print("=" * 70)
print("6. Validation")
print("-" * 70)

print(f"\n   Total rows: {len(df)}")
print(f"   Total columns: {len(df.columns)} (19 core + {len(NEW_COLS)} new)")
print(f"\n   itt_group distribution:")
print(df["itt_group"].value_counts().to_string())

print(f"\n   New column non-null & summary:")
all_empty = []
for c in NEW_COLS:
    non_null = df[c].notna().sum()
    pct = non_null / len(df) * 100
    if pct == 0:
        all_empty.append(c)
    # 对于 0/1 列打印 sum
    if c not in ("cv_arrhythmia_hours", "cv_vasopressor_hours") and pct > 0:
        _sum = df[c].sum()
        print(f"   {c:<35s} non-null={non_null:>5d} ({pct:5.1f}%)  sum={_sum}")
    else:
        _non_null_vals = df[c].dropna()
        if len(_non_null_vals) > 0:
            print(f"   {c:<35s} non-null={non_null:>5d} ({pct:5.1f}%)  "
                  f"range=[{_non_null_vals.min():.1f}, {_non_null_vals.max():.1f}]")
        else:
            print(f"   {c:<35s} non-null={non_null:>5d} ({pct:5.1f}%)  ALL NULL")

if all_empty:
    print(f"\n   ⚠ WARNING: columns completely empty: {all_empty}")
else:
    print(f"\n   ✓ No completely empty columns among new fields.")

# 按 itt_group 分组的新列统计
print(f"\n   New column means by itt_group:")
new_summary = df.groupby("itt_group")[NEW_COLS].mean()
print(new_summary.to_string())

# ═══════════════════════════════════════════════════════════
# 7. 导出
# ═══════════════════════════════════════════════════════════
print(f"\n{'=' * 70}")
print("7. Exporting …")

PARQUET_PATH = OUT_DIR / "sequential_trials_cv.parquet"
CSV_PATH = OUT_DIR / "sequential_trials_cv.csv"
META_PATH = OUT_DIR / "sequential_trials_cv_metadata.json"

# parquet (snappy)
df.to_parquet(PARQUET_PATH, compression="snappy", index=False)
print(f"   ✓ {PARQUET_PATH}")

# CSV
df.to_csv(CSV_PATH, index=False)
print(f"   ✓ {CSV_PATH}")

# ═══════════════════════════════════════════════════════════
# 8. 生成 metadata JSON
# ═══════════════════════════════════════════════════════════
print(f"\n8. Writing metadata JSON …")

meta = {
    "exported_at": datetime.now(timezone.utc).isoformat(),
    "source": str(MAIN_PATH),
    "row_count": int(len(df)),
    "column_count": int(len(df.columns)),
    "core_columns": MAIN_CORE_COLS,
    "new_columns": NEW_COLS,
    "column_dtypes": {c: str(df[c].dtype) for c in df.columns},
    "non_null_fraction": {c: round(float(df[c].notna().sum()) / len(df), 4) for c in df.columns},
    "itt_group_counts": df["itt_group"].value_counts().to_dict(),
    "empty_columns": all_empty,
    "input_files": [
        "sql/result/M0_14d_composite.parquet",
        "sql/result/M0_14a_arrhythmia.parquet",
        "sql/result/M0_16_cv_hemodynamic.parquet",
        "sql/result/M0_14c_myocardial.parquet",
        "sql/result/M0_05f_dim6_hypo_electrolyte.parquet",
    ],
}

with open(META_PATH, "w", encoding="utf-8") as f:
    json.dump(meta, f, indent=2, ensure_ascii=False, default=str)

print(f"   ✓ {META_PATH}")

# ── 最终汇总 ──
print(f"\n{'=' * 70}")
print("DONE — Summary")
print("-" * 70)
print(f"  Input  rows : {EXPECTED_ROWS}")
print(f"  Output rows : {len(df)}")
print(f"  Output cols : {len(df.columns)}  (19 core + {len(NEW_COLS)} cv/electrolyte)")
print(f"  Parquet     : {PARQUET_PATH}")
print(f"  CSV         : {CSV_PATH}")
print(f"  Metadata    : {META_PATH}")
if all_empty:
    print(f"  ⚠ Warnings  : Columns all-null → {all_empty}")
print("=" * 70)