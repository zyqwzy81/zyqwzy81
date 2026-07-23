#!/usr/bin/env python3
"""
构建 24h 洗脱期敏感性分析的时序试验数据集
基于 study_cohort_washout24h 表，输出至 TTE/v28/output/washout24h/
"""
import sys, time, logging, json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# ---- 项目路径 ----
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from framework.config.loader import load_config
from framework.sql_extractor.connection import ConnectionPool

# ---- 输出目录 (v28/washout24h) ----
OUTPUT_DIR = PROJECT_ROOT / "TTE" / "v28" / "output" / "washout24h"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
(OUTPUT_DIR / "logs").mkdir(parents=True, exist_ok=True)

# ---- 日志 ----
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(
            OUTPUT_DIR / "logs" / "build_washout24h.log",
            encoding="utf-8"
        ),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger("washout24h.build")

# ---- 核心分析必需列 ----
REQUIRED_COLUMNS = [
    "stay_id", "subject_id", "trial_num", "time_zero", "itt_group",
    "itt_group_48h",
    "baseline_phosphate", "age", "gender", "baseline_sofa", "egfr",
    "sepsis3", "invasive_vent", "vasoactive", "charlson",
    "admission_type", "pao2fio2ratio", "sf_ratio",
    "bmi", "race_group", "sapsii", "gcs", "liver_disease",
    "lactate_t0", "albumin_min", "hemoglobin_min", "platelets_min", "wbc_max",
    "icu_type", "sofa_t0",
    "death_28day_num", "icu_los_days",
    # 精确时间列 (用于边界精确对位)
    "icu_intime", "icu_outtime", "deathtime", "dischtime",
    "death_time_from_t0", "dischtime_from_t0", "surv_time",
    "enrollment_seq",
]

# ---- 缺失填充列 ----
FILL_ZERO_COLS = [
    "diuretic", "steroid", "insulin", "nahco3",
    "beta2_agonist", "ppi", "alkalosis",
    "transfusion", "cardiac_surgery", "aki", "sepsis3",
    "invasive_vent", "vasoactive",
]


def load_main_data(pool: ConnectionPool) -> pd.DataFrame:
    """Step 1: 加载 study_cohort_washout24h + 基线协变量 + 结局"""
    logger.info("Step 1: 加载 study_cohort_washout24h + baseline_covariates + outcome_data...")
    t0 = time.time()

    with pool.cursor() as cur:
        cur.execute("""
            SELECT sc.stay_id, sc.subject_id, sc.hadm_id, sc.time_zero,
                   sc.baseline_phosphate,
                   sc.itt_group, sc.itt_group_48h, sc.enrollment_seq,
                   sc.first_iv_time, sc.iv_within_12h,
                   bc.age, bc.gender,
                   bc.admission_sofa AS baseline_sofa,
                   bc.sepsis3, bc.invasive_vent_at_t0,
                   bc.charlson_comorbidity_index,
                   bc.egfr,
                   bc.vasoactive_at_t0,
                   bc.admission_type,
                   bc.pao2fio2ratio, bc.sf_ratio, bc.sf_estimated_pf, bc.sf_fio2_defaulted,
                   bc.bmi, bc.race_group,
                   bc.sapsii, bc.gcs, bc.liver_disease,
                   bc.lactate_t0, bc.albumin_min,
                   bc.hemoglobin_min, bc.platelets_min, bc.wbc_max,
                   bc.icu_type,
                   bc.sofa_t0,
                   oc.death_28day,
                   oc.icu_los_days,
                   a.deathtime,
                   a.dischtime,
                   sc.icu_intime, sc.icu_outtime
            FROM study_cohort_washout24h sc
            LEFT JOIN (
                SELECT DISTINCT ON (stay_id, time_zero)
            stay_id, time_zero,
            age, gender,
            admission_sofa,
            sepsis3, invasive_vent_at_t0,
            charlson_comorbidity_index,
            egfr, vasoactive_at_t0,
            admission_type,
            pao2fio2ratio, sf_ratio, sf_estimated_pf, sf_fio2_defaulted,
            bmi, race_group,
             sapsii, gcs, liver_disease,
             lactate_t0, albumin_min,
             hemoglobin_min, platelets_min, wbc_max,
             icu_type,
             sofa_t0
        FROM baseline_covariates
            ) bc ON sc.stay_id = bc.stay_id AND sc.time_zero = bc.time_zero
            LEFT JOIN mimiciv_hosp.admissions a
                ON sc.hadm_id = a.hadm_id
            LEFT JOIN outcome_data oc
                ON sc.stay_id = oc.stay_id AND sc.time_zero = oc.time_zero
        """)
        df = pd.DataFrame(cur.fetchall(), columns=[desc[0] for desc in cur.description])

    logger.info(f"  加载: {len(df):,} 行 × {len(df.columns)} 列 ({time.time() - t0:.1f}s)")
    return df


def check_table_exists(pool: ConnectionPool, table_name: str) -> bool:
    with pool.cursor() as cur:
        cur.execute("SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = %s)", (table_name,))
        return cur.fetchone()[0]


def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    """Step 3: 数据清洗"""
    logger.info("Step 3: 数据清洗...")
    t0 = time.time()

    n_before = len(df)

    # 去重
    dup_mask = df.duplicated(subset=["stay_id", "time_zero"], keep="first")
    n_dup = dup_mask.sum()
    if n_dup > 0:
        df = df[~dup_mask].copy()
        logger.info(f"  去重: {n_dup} 行 (保留 {len(df)} 行)")

    # 重命名 charlson
    if "charlson_comorbidity_index" in df.columns:
        df.rename(columns={"charlson_comorbidity_index": "charlson"}, inplace=True)

    # 二值变量缺失填充0
    for col in FILL_ZERO_COLS:
        if col in df.columns:
            df[col] = df[col].fillna(0)

    # 结局变量: 优先使用 outcome_data 的 death_28day
    # 对于 washout24h 新增的 trial periods (outcome_data 无匹配), 从 admissions 计算
    if "death_28day" in df.columns:
        # 标记 outcome_data 缺失的行
        death_from_oc = df["death_28day"].notna()
        n_oc_match = death_from_oc.sum()
        # 从 admissions 表计算 death_28day
        df["death_28day"] = df["death_28day"].fillna(0)
        # 对于 outcome_data 缺失的行，从 deathtime 计算
        missing_mask = ~death_from_oc
        if missing_mask.any():
            death_from_deathtime = (
                df.loc[missing_mask, "deathtime"].notna()
                & (df.loc[missing_mask, "deathtime"] <= df.loc[missing_mask, "time_zero"] + pd.Timedelta(days=28))
            )
            df.loc[missing_mask, "death_28day"] = death_from_deathtime.astype(int)
            n_death_computed = death_from_deathtime.sum()
            logger.info(f"  outcome_data 匹配: {n_oc_match}, 从deathtime补充: {n_death_computed}")
    else:
        # 从 admissions 直接计算
        df["death_28day"] = (
            df["deathtime"].notna()
            & (df["deathtime"] <= df["time_zero"] + pd.Timedelta(days=28))
        ).astype(int)
    df["death_28day"] = df["death_28day"].astype(int)

    # 精确时间计算
    df["time_zero"] = pd.to_datetime(df["time_zero"])
    df["deathtime"] = pd.to_datetime(df["deathtime"])
    df["dischtime"] = pd.to_datetime(df["dischtime"])
    df["icu_intime"] = pd.to_datetime(df["icu_intime"])
    df["icu_outtime"] = pd.to_datetime(df["icu_outtime"])

    # death_time_from_t0 (小时)
    df["death_time_from_t0"] = (df["deathtime"] - df["time_zero"]).dt.total_seconds() / 3600
    df["dischtime_from_t0"] = (df["dischtime"] - df["time_zero"]).dt.total_seconds() / 3600

    # iv_within_24h
    df["iv_within_24h"] = 0
    iv_mask = df["first_iv_time"].notna()
    df.loc[iv_mask, "iv_within_24h"] = (
        (pd.to_datetime(df.loc[iv_mask, "first_iv_time"]) - df.loc[iv_mask, "time_zero"]).dt.total_seconds() / 3600 <= 24
    ).astype(int)

    # surv_time
    df["surv_time"] = np.where(
        df["death_28day"] == 1,
        df["death_time_from_t0"].clip(lower=0),
        df["dischtime_from_t0"].clip(lower=0)
    )
    df["surv_time"] = df["surv_time"].fillna(28 * 24)
    df["surv_time"] = df["surv_time"].clip(upper=28 * 24)

    # 重命名为 R 兼容格式
    df.rename(columns={"death_28day": "death_28day_num"}, inplace=True)

    # trial_num
    df["trial_num"] = df["enrollment_seq"]

    logger.info(f"  清洗完成: {len(df)} 行 ({time.time() - t0:.1f}s)")
    return df


def report_summary(df: pd.DataFrame):
    """Step 4: 报告摘要"""
    logger.info("Step 4: 数据摘要...")
    total = len(df)
    n_death = int(df["death_28day_num"].sum())
    n_alive = total - n_death
    n_treatment = int((df["itt_group"] == "treatment").sum())
    n_conservative = total - n_treatment
    n_unique_pt = df["subject_id"].nunique()
    n_unique_stay = df["stay_id"].nunique()

    logger.info("=" * 50)
    logger.info("24h 洗脱期敏感性分析数据集摘要")
    logger.info(f"  总试验期: {total:,}")
    logger.info(f"  治疗组: {n_treatment:,} ({100*n_treatment/total:.1f}%)")
    logger.info(f"  对照组: {n_conservative:,} ({100*n_conservative/total:.1f}%)")
    logger.info(f"  28天死亡: {n_death:,} ({100*n_death/total:.1f}%)")
    logger.info(f"  存活: {n_alive:,} ({100*n_alive/total:.1f}%)")
    logger.info(f"  唯一患者: {n_unique_pt:,}")
    logger.info(f"  唯一ICU住院: {n_unique_stay:,}")
    logger.info("=" * 50)


def export_data(df: pd.DataFrame):
    """Step 5-6: 导出 Parquet + CSV + 元数据"""
    logger.info("Step 5: 导出数据...")
    t0 = time.time()

    # 核心列
    avail_cols = [c for c in REQUIRED_COLUMNS if c in df.columns]
    df_core = df[avail_cols].copy()

    # 导出 Parquet
    core_parquet_path = OUTPUT_DIR / "sequential_trials_washout24h.parquet"
    df_core.to_parquet(core_parquet_path, index=False)
    logger.info(f"  核心 Parquet: {core_parquet_path} ({len(df_core):,} 行 × {len(df_core.columns)} 列)")

    core_csv_path = OUTPUT_DIR / "sequential_trials_washout24h.csv"
    df_core.to_csv(core_csv_path, index=False)
    logger.info(f"  核心 CSV: {core_csv_path} ({core_csv_path.stat().st_size / 1024 / 1024:.1f} MB)")

    # 元数据
    meta = {
        "dataset": "sequential_trials_washout24h",
        "description": "24h washout sensitivity analysis dataset",
        "created": datetime.now().isoformat(),
        "n_rows": len(df_core),
        "n_columns": len(df_core.columns),
        "columns": list(df_core.columns),
        "n_treatment": int((df_core["itt_group"] == "treatment").sum()),
        "n_conservative": int((df_core["itt_group"] == "conservative").sum()),
        "n_death_28day": int(df_core["death_28day_num"].sum()),
        "n_unique_subject": int(df_core["subject_id"].nunique()),
        "n_unique_stay": int(df_core["stay_id"].nunique()),
    }
    meta_path = OUTPUT_DIR / "sequential_trials_washout24h_meta.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    logger.info(f"  元数据: {meta_path}")

    logger.info(f"  导出完成 ({time.time() - t0:.1f}s)")


def main():
    config = load_config(version="v26")
    logger.info(f"Config loaded: dbname={config.database.dbname}")

    pool = ConnectionPool(config)

    # Step 1: 加载数据
    df = load_main_data(pool)

    # Step 2: 合并混杂变量 (简化为必需)
    logger.info("Step 2: 合并混杂变量 (跳过, 使用核心变量)")

    # Step 3: 清洗
    df = clean_data(df)

    # Step 4: 报告
    report_summary(df)

    # Step 5-6: 导出
    export_data(df)

    pool.close()
    logger.info("全部完成")


if __name__ == "__main__":
    main()