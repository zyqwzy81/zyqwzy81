"""v26 序贯试验分析数据集构建脚本

从 PostgreSQL 表构建 sequential_trials 分析数据集，直接导出 Parquet。

数据源（9张表）：
  study_cohort + baseline_covariates + outcome_data
  + medication_confounders + acid_base_confounders + procedure_transfusion_confounders
  + nutrition_refeeding_confounders + crrt_confounders + calcium_confounders

输出路径：
  TTE/v26/output/sequential_trials.parquet   （主格式，分析级数据）
  TTE/v26/output/sequential_trials.csv       （CSV 备份，R 脚本回退兼容）

方法学变更（V26修订，2026-06-12）：
  ⚠️ 所有协变量改为以 ICU 入院为基线参考点，而非 T0（timezero）
     - baseline_sofa: 使用 admission_sofa（first_day_sofa）替代 T0±24h 匹配
     - egfr: 基于 ICU 入院24h内 creatinine，回退至 T0 值
     - invasive_vent: ICU 入院24h内启动有创通气（非 T0 是否正在通气）
     - vasoactive: ICU 入院24h内使用血管活性药（非 T0 是否在使用）
     - sepsis3: 本身为 per-stay 判定，不受影响
     - charlson: 本身为 per-hadm_id 判定，不受影响
  理由：同一 ICU stay 内序贯试验共享同一个ICU入院基线值，
        避免同一患者因序贯入组时间不同而获得不同的"基线"协变量值。

使用方式：
  python TTE/v26/code/python/build_sequential_trials.py

依赖：
  - framework.config.loader (load_config)
  - framework.sql_extractor.connection (ConnectionPool)
  - pandas, numpy, pyarrow
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

# ---- 输出目录 ----
OUTPUT_DIR = PROJECT_ROOT / "TTE" / "v26" / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
(OUTPUT_DIR / "logs").mkdir(parents=True, exist_ok=True)

# ---- 日志 ----
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(
            OUTPUT_DIR / "logs" / "v26_build_sequential_trials.log",
            encoding="utf-8"
        ),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger("v26.build")

# ---- 核心分析必需列（19列，SEQTaRget 兼容 + 精确时间列） ----
REQUIRED_COLUMNS = [
    "stay_id", "subject_id", "trial_num", "time_zero", "itt_group",
    "itt_group_48h",  # V26.11: 48h宽限期定义 (敏感性分析保留)
    "baseline_phosphate", "age", "gender", "baseline_sofa", "egfr",
    "sepsis3", "invasive_vent_at_t0", "vasoactive_at_t0",
    "charlson", "admission_type",
    "death_28day_num", "surv_time",
    "death_time_from_t0", "dischtime_from_t0",  # 🆕 P0-M05: 精确时间列
]

# ---- 全协变量增强列（V25.3 可直接用于全协变量手动IPTW分析，含营养变量） ----
AUGMENTED_COLUMNS = [
    "stay_id", "subject_id", "trial_num", "time_zero", "itt_group",
    "itt_group_48h",  # V26.11: 48h宽限期定义 (敏感性分析保留)
    "baseline_phosphate", "age", "gender", "baseline_sofa", "egfr",
    "sepsis3", "invasive_vent_at_t0", "vasoactive_at_t0",
    "charlson", "admission_type",
    "death_28day_num", "surv_time",
    "death_time_from_t0", "dischtime_from_t0",
    "first_iv_time", "iv_within_12h", "iv_within_24h",
    # 药物混杂
    "diuretic_at_t0", "steroid_at_t0", "insulin_at_t0",
    "nahco3_infusion_at_t0", "beta2_agonist_at_t0", "ppi_use",
    # 酸碱平衡
    "ph_worst", "bicarbonate_worst",
    "alkalosis", "resp_alkalosis", "metabolic_alkalosis",
    # 操作/输血
    "net_fluid_balance_24h", "transfusion_at_t0", "cardiac_surgery", "aki_kdigo2_plus",
    # 再喂养风险 + 营养
    "refeeding_risk", "tpn_kcal_24h", "en_kcal_24h", "nutrition_support_type",
    "en_start_hours_from_t0", "tpn_start_hours_from_t0",
    "gnri", "pni", "conut_score",
    # CRRT
    "ever_crrt",
    # 钙剂
    "corrected_calcium_mmol_min",
    # 氧合指标
    "pao2fio2ratio", "sf_ratio", "sf_estimated_pf", "sf_fio2_defaulted",
    # ===== V26.3 补充变量 =====
    # 人口学补充
    "bmi", "race_group",
    # 疾病严重度补充
    "sapsii", "gcs", "liver_disease",
    # 实验室补充
    "lactate_t0", "albumin_min", "hemoglobin_min", "platelets_min", "wbc_max",
    # 电解质补充
    "potassium_t0", "sodium_t0", "bilirubin_t0", "bun_t0", "magnesium_t0",
    # 凝血补充
    "inr_max", "ptt_max",
    # ICU类型
    "icu_type",
    # ===== V26.4: sofa_t0 (T0时刻SOFA) + 心血管事件 =====
    "sofa_t0",
    "cardiovascular_event_2comp", "cardiovascular_event_3comp", "cardiovascular_event_comprehensive",
    "new_arrhythmia_rhythm", "vasopressor_escalation", "myocardial_injury", "ecg_abnormal",
    # ===== V26.4: 口服磷完善 =====
    "po_formulation", "cumulative_po_phosphate_mmol_48h", "po_route_category",
    "first_po_starttime", "first_po_drugname", "first_po_route",
    "po_started_within_grace", "po_administrations",
    # ===== V26.10: IV 剂量/药物/速率 =====
    "cumulative_iv_phosphate_mmol_48h",
    "first_iv_drug_name", "first_iv_hours_from_tz",
    "first_iv_rate", "first_iv_amount",
    "exposure_group",
]

# ---- 需要填充缺失为 0 的二值变量 ----
FILL_ZERO_COLS = [
    "diuretic_at_t0", "steroid_at_t0", "insulin_at_t0",
    "nahco3_infusion_at_t0", "beta2_agonist_at_t0", "ppi_use",
    "alkalosis", "resp_alkalosis", "metabolic_alkalosis",
    "transfusion_at_t0", "cardiac_surgery", "aki_kdigo2_plus",
    "sepsis3", "invasive_vent_at_t0", "vasoactive_at_t0",
]

# ---- 需要合并的混杂变量表 ----
MERGE_TABLES = [
    {
        "table": "medication_confounders",
        "label": "药物混杂",
        "columns": [
            "stay_id", "time_zero", "diuretic_at_t0", "steroid_at_t0",
            "insulin_at_t0", "nahco3_infusion_at_t0", "beta2_agonist_at_t0",
            "ppi_use",
        ],
        "dedup": False,
    },
    {
        "table": "acid_base_confounders",
        "label": "酸碱平衡",
        "columns": [
            "stay_id", "time_zero", "ph_worst", "bicarbonate_worst",
            "alkalosis", "resp_alkalosis", "metabolic_alkalosis",
        ],
        "dedup": True,  # 该表有重复 (stay_id, time_zero)
    },
    {
        "table": "procedure_transfusion_confounders",
        "label": "操作/输血",
        "columns": [
            "stay_id", "time_zero", "net_fluid_balance_24h",
            "transfusion_at_t0", "cardiac_surgery", "aki_kdigo2_plus",
        ],
        "dedup": False,
    },
    # ---- V25.3: 从实际存在的维度表提取（wrapper表从未创建，直接引用源表） ----
    {
        "table": "baseline_dim7_nutrition",
        "label": "再喂养风险+营养(全9变量)",
        "columns": [
            "stay_id", "time_zero",
            "refeeding_risk",
            "tpn_kcal_24h", "en_kcal_24h", "nutrition_support_type",
            "en_start_hours_from_t0", "tpn_start_hours_from_t0",
            "gnri", "pni", "conut_score",
        ],
        "dedup": True,  # dim7 表可能有重复 (stay_id, time_zero)
    },
    {
        "table": "crrt_flow_and_clearance",
        "label": "CRRT",
        "columns": [
            "stay_id", "time_zero", "ever_crrt",
        ],
        "dedup": True,
    },
    {
        "table": "baseline_dim6_hypo_electrolyte",
        "label": "钙剂补充",
        "columns": [
            "stay_id", "time_zero", "corrected_calcium_mmol_min",
        ],
        "dedup": True,
    },
    # ===== V26.3 补充变量合并 =====
    {
        "table": "baseline_labs_t0",
        "label": "实验室补充(电解质+肝肾)",
        "columns": [
            "stay_id", "time_zero",
            "potassium_t0", "sodium_t0", "bilirubin_t0", "bun_t0", "magnesium_t0",
        ],
        "dedup": True,
    },
    {
        "table": "cbc_coagulation",
        "label": "凝血补充(INR+PTT)",
        "columns": [
            "stay_id", "time_zero",
            "inr_max", "ptt_max",
        ],
        "dedup": True,
        "where": "row_type = 'baseline'",
    },
    # ===== V26.4: 心血管事件合并 =====
    {
        "table": "cv_composite",
        "label": "心血管事件复合终点",
        "columns": [
            "stay_id", "time_zero",
            "cardiovascular_event_2comp", "cardiovascular_event_3comp",
            "cardiovascular_event_comprehensive",
        ],
        "dedup": True,
    },
    {
        "table": "cv_arrhythmia",
        "label": "心律失常",
        "columns": [
            "stay_id", "time_zero",
            "new_arrhythmia_rhythm",
        ],
        "dedup": True,
    },
    {
        "table": "cv_hemodynamic",
        "label": "血流动力学",
        "columns": [
            "stay_id", "time_zero",
            "vasopressor_escalation",
        ],
        "dedup": True,
    },
    {
        "table": "cv_myocardial",
        "label": "心肌损伤",
        "columns": [
            "stay_id", "time_zero",
            "myocardial_injury", "ecg_abnormal",
        ],
        "dedup": True,
    },
    # ===== V26.4: 口服磷完善列 (从 exposure_data 提取) =====
    {
        "table": "exposure_data",
        "label": "口服磷+IV暴露",
        "columns": [
            "stay_id", "time_zero",
            # PO 列
            "po_formulation", "cumulative_po_phosphate_mmol_48h",
            "po_route_category", "first_po_starttime",
            "first_po_drugname", "first_po_route",
            "po_started_within_grace", "po_administrations",
            # V26.10: IV 剂量/药物/速率列
            "cumulative_iv_phosphate_mmol_48h",
            "first_iv_drug_name",
            "first_iv_hours_from_tz",
            "first_iv_rate",
            "first_iv_amount",
            "exposure_group",
        ],
        "dedup": True,
    },
]


def load_main_data(pool: ConnectionPool) -> pd.DataFrame:
    """Step 1: 加载主队列 + 基线协变量 + 结局 + 精确死亡/出院时间"""
    logger.info("Step 1: 加载 study_cohort + baseline_covariates + outcome_data + admissions...")
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
                   -- V26.3: 补充变量
                   bc.bmi, bc.race_group,
                   bc.sapsii, bc.gcs, bc.liver_disease,
                   bc.lactate_t0, bc.albumin_min,
                   bc.hemoglobin_min, bc.platelets_min, bc.wbc_max,
                   bc.icu_type,
                   bc.sofa_t0,
                   oc.death_28day,
                   oc.icu_los_days,
                   a.deathtime,
                   a.dischtime
            FROM study_cohort sc
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
            LEFT JOIN outcome_data oc
                ON sc.stay_id = oc.stay_id AND sc.time_zero = oc.time_zero
            LEFT JOIN mimiciv_hosp.admissions a
                ON sc.hadm_id = a.hadm_id
        """)
        df = pd.DataFrame(cur.fetchall(), columns=[desc[0] for desc in cur.description])

    logger.info(f"  加载: {len(df):,} 行 × {len(df.columns)} 列 ({time.time() - t0:.1f}s)")
    return df


def check_table_exists(cur, table_name: str) -> bool:
    """检查表是否存在于 public schema 中"""
    try:
        cur.execute(
            "SELECT EXISTS(SELECT 1 FROM information_schema.tables "
            "WHERE table_schema = 'public' AND table_name = %s)",
            (table_name,)
        )
        return cur.fetchone()[0]
    except Exception:
        return False


def merge_confounders(df: pd.DataFrame, pool: ConnectionPool) -> pd.DataFrame:
    """Step 2: 合并 v26 新增混杂变量（缺失表容错跳过，不填NaN）"""
    for cfg in MERGE_TABLES:
        t0 = time.time()
        try:
            with pool.cursor() as cur:
                # 先检查表是否存在
                if not check_table_exists(cur, cfg["table"]):
                    logger.warning(f"  {cfg['label']}: 表 {cfg['table']} 不存在，跳过（不填NaN）")
                    # ❌ 不填 NaN — 缺失列由下游管线的 MISSING_COLUMNS 列表管理
                    continue

                if cfg["dedup"]:
                    where_clause = ""
                    if cfg.get("where"):
                        where_clause = f" WHERE {cfg['where']}"
                    cur.execute(
                        f"SELECT DISTINCT ON (stay_id, time_zero) {', '.join(cfg['columns'])} "
                        f"FROM {cfg['table']}{where_clause} ORDER BY stay_id, time_zero"
                    )
                else:
                    cur.execute(f"SELECT {', '.join(cfg['columns'])} FROM {cfg['table']}")
                df_merge = pd.DataFrame(
                    cur.fetchall(), columns=[desc[0] for desc in cur.description]
                )

            new_cols = [
                c
                for c in df_merge.columns
                if c not in ["stay_id", "time_zero", "subject_id", "hadm_id"]
            ]
            df = df.merge(
                df_merge[["stay_id", "time_zero"] + new_cols],
                on=["stay_id", "time_zero"],
                how="left",
            )

            # 报告各列合并后的非空率
            for col in new_cols:
                n_valid = df[col].notna().sum()
                pct = 100 * n_valid / len(df)
                logger.info(f"    {col}: {n_valid:,}/{len(df):,} ({pct:.1f}%)")

            logger.info(f"  {cfg['label']}: +{len(new_cols)} 列 ({time.time() - t0:.1f}s)")
        except Exception as e:
            logger.warning(f"  {cfg['label']}: 跳过（查询失败: {e}），不填NaN")
            # ❌ 不再填 NaN — 缺失列由下游管线的 MISSING_COLUMNS 列表管理
            continue
    return df


def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    """Step 3: 数据清洗、去重、缺失填充、类型修正"""
    logger.info("Step 3: 数据清洗...")

    # 去重（study_cohort 有 6,181 条重复，DISTINCT ON (stay_id, time_zero)）
    n_before = len(df)
    df = df.drop_duplicates(subset=["stay_id", "time_zero"], keep="first")
    n_dedup = n_before - len(df)
    if n_dedup > 0:
        logger.info(f"  去重: {n_before:,} → {len(df):,} (-{n_dedup:,})")

    # 缩短长变量名（SEQTaRget 公式缓冲区 ~200 字符限制）
    if "charlson_comorbidity_index" in df.columns:
        df = df.rename(columns={"charlson_comorbidity_index": "charlson"})
        logger.info("  重命名: charlson_comorbidity_index → charlson")

    # 二值变量缺失填充 0
    for col in FILL_ZERO_COLS:
        if col in df.columns:
            df[col] = df[col].fillna(0).astype(int)

    # 结局变量推导
    if "death_28day" in df.columns:
        df["death_28day_num"] = df["death_28day"].fillna(0).astype(int)

    # 🆕 P0-M05: 精确死亡/出院时间计算（从 admissions 原始时间戳推导）
    if "time_zero" in df.columns and "deathtime" in df.columns and "dischtime" in df.columns:
        time_zero_dt = pd.to_datetime(df["time_zero"])
        deathtime_dt = pd.to_datetime(df["deathtime"])
        dischtime_dt = pd.to_datetime(df["dischtime"])

        # death_time_from_t0: T0 到死亡的天数（仅死亡患者）
        df["death_time_from_t0"] = (deathtime_dt - time_zero_dt).dt.total_seconds() / 86400.0
        # 无效值处理：死亡时间早于T0或超出28天范围
        df.loc[df["death_time_from_t0"] < 0, "death_time_from_t0"] = np.nan
        df.loc[df["death_time_from_t0"] > 28, "death_time_from_t0"] = 28.0

        # dischtime_from_t0: T0 到出院的天数（所有患者）
        df["dischtime_from_t0"] = (dischtime_dt - time_zero_dt).dt.total_seconds() / 86400.0
        df.loc[df["dischtime_from_t0"] < 0, "dischtime_from_t0"] = np.nan
        # 出院时间截断至28天
        df["dischtime_from_t0"] = df["dischtime_from_t0"].clip(upper=28.0)

        n_death_time = df["death_time_from_t0"].notna().sum()
        n_disch_time = df["dischtime_from_t0"].notna().sum()
        logger.info(f"  death_time_from_t0: {n_death_time:,}/{len(df):,} 有效 ({100*n_death_time/len(df):.1f}%)")
        logger.info(f"  dischtime_from_t0: {n_disch_time:,}/{len(df):,} 有效 ({100*n_disch_time/len(df):.1f}%)")

        # 🆕 iv_within_24h: 首次补磷时间是否在T0后24h内（宽限期敏感性分析用）
        if "first_iv_time" in df.columns:
            first_iv_dt = pd.to_datetime(df["first_iv_time"])
            time_diff_hours = (first_iv_dt - time_zero_dt).dt.total_seconds() / 3600.0
            df["iv_within_24h"] = ((time_diff_hours <= 24) & (time_diff_hours >= 0)).astype(int)
            # first_iv_time 为 NaT 的行填 0
            df["iv_within_24h"] = df["iv_within_24h"].fillna(0).astype(int)
            n_within_24h = df["iv_within_24h"].sum()
            logger.info(f"  iv_within_24h: {n_within_24h:,}/{len(df):,} ({100*n_within_24h/len(df):.1f}%)")
    else:
        logger.warning("  deathtime/dischtime 列缺失，无法计算精确时间列")
        df["death_time_from_t0"] = np.nan
        df["dischtime_from_t0"] = np.nan
        df["iv_within_24h"] = 0

    if "surv_time" not in df.columns:
        # 🆕 修正 surv_time 推导（替代 icu_los_days 代理法）
        # 优先级: (1) 精确死亡时间 (2) 出院时间 (3) ICU LOS (4) 默认28天
        if "death_time_from_t0" in df.columns and df["death_time_from_t0"].notna().any():
            # 死亡患者：使用实际死亡时间；存活患者：使用出院时间或ICU LOS
            df["surv_time"] = df.apply(
                lambda row: min(row["death_time_from_t0"], 28.0)
                if pd.notna(row["death_time_from_t0"]) and row.get("death_28day_num", 0) == 1
                else (
                    min(row["dischtime_from_t0"], 28.0)
                    if pd.notna(row.get("dischtime_from_t0"))
                    else min(row.get("icu_los_days", 28.0) or 28.0, 28.0)
                ), axis=1
            )
            logger.info("  surv_time: 基于 death_time_from_t0 + dischtime_from_t0 精确计算")
        elif "dischtime_from_t0" in df.columns and df["dischtime_from_t0"].notna().any():
            # 存活患者：最后已知存活时间 + 死亡患者：实际死亡时间
            df["surv_time"] = df.apply(
                lambda row: min(row.get("death_time_from_t0", 28.0) or 28.0, 28.0)
                if row.get("death_28day_num", 0) == 1 and pd.notna(row.get("death_time_from_t0"))
                else min(row["dischtime_from_t0"], 28.0), axis=1
            )
            logger.info("  surv_time: 基于 dischtime_from_t0 计算（存活患者）")
        elif "icu_los_days" in df.columns:
            # 回退方案：使用 ICU LOS
            df["surv_time"] = df["icu_los_days"].clip(lower=0).fillna(28)
            logger.warning("  surv_time: 回退使用 icu_los_days（非精确方案，可能引入偏倚）")
        else:
            df["surv_time"] = 28.0
            logger.warning("  surv_time: 默认28天（无任何时间列可用，严重偏倚风险！）")
        
        # 添加 surv_time 推导来源标记
        df["surv_time_source"] = "derived"

    # trial_num 从 enrollment_seq 推导
    if "trial_num" not in df.columns:
        if "enrollment_seq" in df.columns:
            df["trial_num"] = (
                df.groupby("stay_id")["enrollment_seq"]
                .rank(method="dense")
                .astype(int)
            )
        else:
            df["trial_num"] = 1

    # 类型标准化（确保与 REQUIRED_TYPES 兼容）
    int_cols = ["stay_id", "subject_id", "hadm_id", "enrollment_seq"]
    float_cols = [
        "age", "baseline_sofa", "charlson",
        "egfr", "baseline_phosphate",
    ]

    # ===== V26.3 数据清洗：异常值截断 =====
    # magnesium_t0: 正常范围 0.5-5.0 mmol/L，超出设为缺失
    if "magnesium_t0" in df.columns:
        n_before = df["magnesium_t0"].notna().sum()
        df.loc[df["magnesium_t0"] > 5.0, "magnesium_t0"] = np.nan
        df.loc[df["magnesium_t0"] < 0.5, "magnesium_t0"] = np.nan
        n_after = df["magnesium_t0"].notna().sum()
        n_clipped = n_before - n_after
        if n_clipped > 0:
            logger.info(f"  magnesium_t0: 截断 {n_clipped} 个异常值（正常范围 0.5-5.0 mmol/L）")

    # corrected_calcium_mmol_min: 正常范围 1.0-5.0 mmol/L，超出设为缺失
    if "corrected_calcium_mmol_min" in df.columns:
        n_before = df["corrected_calcium_mmol_min"].notna().sum()
        df.loc[df["corrected_calcium_mmol_min"] > 5.0, "corrected_calcium_mmol_min"] = np.nan
        df.loc[df["corrected_calcium_mmol_min"] < 1.0, "corrected_calcium_mmol_min"] = np.nan
        n_after = df["corrected_calcium_mmol_min"].notna().sum()
        n_clipped = n_before - n_after
        if n_clipped > 0:
            logger.info(f"  corrected_calcium_mmol_min: 截断 {n_clipped} 个异常值（正常范围 1.0-5.0 mmol/L）")

    # en_kcal_24h: Winsorize 99% 分位数截断
    if "en_kcal_24h" in df.columns:
        upper = df["en_kcal_24h"].quantile(0.99)
        n_before = (df["en_kcal_24h"] > upper).sum()
        df["en_kcal_24h"] = df["en_kcal_24h"].clip(upper=upper)
        if n_before > 0:
            logger.info(f"  en_kcal_24h: Winsorize 99% 截断（上限={upper:.0f} kcal），{n_before} 个值受影响")

    # net_fluid_balance_24h: Winsorize 99% 分位数截断
    if "net_fluid_balance_24h" in df.columns:
        upper = df["net_fluid_balance_24h"].quantile(0.99)
        lower = df["net_fluid_balance_24h"].quantile(0.01)
        n_before = (df["net_fluid_balance_24h"].abs() > max(abs(upper), abs(lower))).sum()
        df["net_fluid_balance_24h"] = df["net_fluid_balance_24h"].clip(lower=lower, upper=upper)
        if n_before > 0:
            logger.info(f"  net_fluid_balance_24h: Winsorize 99% 截断（范围=[{lower:.0f}, {upper:.0f}] mL），{n_before} 个值受影响")
    # V26.4 营养支持派生二分类变量（TPN/EN/营养不良）
    # 这些变量在 V26.3 的 R 层已有定义，现上移至数据工程层派生以保证可重复性
    if "tpn_kcal_24h" in df.columns:
        df["tpn_at_t0"] = (df["tpn_kcal_24h"] > 0).astype("int64")
        logger.info(f"  tpn_at_t0: 派生完成，阳性 {(df['tpn_at_t0'] == 1).sum():,} 例")
    if "en_kcal_24h" in df.columns:
        df["en_at_t0"] = (df["en_kcal_24h"] > 0).astype("int64")
        logger.info(f"  en_at_t0: 派生完成，阳性 {(df['en_at_t0'] == 1).sum():,} 例")
    if "albumin_min" in df.columns and "bmi" in df.columns:
        low_albumin = df["albumin_min"] < 3.0
        low_bmi = df["bmi"] < 18.5
        has_albumin = df["albumin_min"].notna()
        has_bmi = df["bmi"].notna()
        df["malnutrition_at_t0"] = np.where(
            low_albumin | low_bmi, 1,
            np.where(has_albumin | has_bmi, 0, np.nan)
        ).astype("float64")  # 保留 NA，后续 MICE 插补
        n_pos = int((df["malnutrition_at_t0"] == 1).sum())
        logger.info(f"  malnutrition_at_t0: 派生完成，阳性 {n_pos:,} 例")

    # V26.4 字符型变量类别压缩（减少 LASSO dummy 数量和 SEQTaRget 内存峰值）
    if "admission_type" in df.columns:
        admission_map = {
            "EW EMER.": "Emergency",
            "URGENT": "Emergency",
            "DIRECT EMER.": "Emergency",
            "OBSERVATION ADMIT": "Observation",
            "EU OBSERVATION": "Observation",
            "DIRECT OBSERVATION": "Observation",
            "AMBULATORY OBSERVATION": "Observation",
            "SURGICAL SAME DAY ADMISSION": "Surgical_Elective",
            "ELECTIVE": "Surgical_Elective",
        }
        df["admission_type"] = df["admission_type"].map(admission_map).astype("str")
        logger.info(f"  admission_type: 压缩为 3 类 -> {df['admission_type'].value_counts().to_dict()}")
    if "race_group" in df.columns:
        race_map = {
            "White": "White",
            "Black": "Black",
            "Hispanic": "Other",
            "Asian": "Other",
            "Other": "Other",
        }
        df["race_group"] = df["race_group"].map(race_map).astype("str")
        logger.info(f"  race_group: 压缩为 3 类 -> {df['race_group'].value_counts().to_dict()}")
    if "icu_type" in df.columns:
        # 先合并 MICE 阶段的小类别
        icu_tab = df["icu_type"].value_counts()
        small_levels = icu_tab[icu_tab < 100].index.tolist()
        if small_levels:
            df["icu_type"] = df["icu_type"].replace(small_levels, "Other_ICU")
        icu_map = {
            "MICU": "Medical",
            "MICU/SICU": "Medical",
            "Medicine": "Medical",
            "SICU": "Surgical",
            "TSICU": "Surgical",
            "Neuro SICU": "Surgical",
            "Surgery/Vascular": "Surgical",
            "CCU": "Cardiac",
            "CVICU": "Cardiac",
            "Neuro Intermediate": "Neuro",
            "Neuro Stepdown": "Neuro",
            "PACU": "Other_ICU",
            "Other": "Other_ICU",
        }
        df["icu_type"] = df["icu_type"].map(icu_map).astype("str")
        logger.info(f"  icu_type: 压缩为 {df['icu_type'].nunique()} 类 -> {df['icu_type'].value_counts().to_dict()}")

    for col in int_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("int64")
    for col in float_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("float64")

    logger.info(f"  清洗后: {len(df):,} 行 × {len(df.columns)} 列")
    return df


def report_summary(df: pd.DataFrame) -> None:
    """Step 4: 输出变量摘要"""
    logger.info("=" * 70)
    logger.info("变量缺失率:")
    for c in sorted(df.columns):
        n_miss = df[c].isna().sum()
        if n_miss > 0:
            logger.info(f"  {c}: {n_miss:,} ({100 * n_miss / len(df):.1f}%)")

    if "death_28day_num" in df.columns and "itt_group" in df.columns:
        logger.info("=" * 70)
        logger.info("结局分布:")
        for grp in sorted(df["itt_group"].dropna().unique()):
            s = df[df["itt_group"] == grp]
            n = int(s["death_28day_num"].sum())
            logger.info(
                f"  {grp}: N={len(s):,}, events={n:,} "
                f"({100 * n / max(len(s), 1):.1f}%)"
            )


def export_data(df: pd.DataFrame) -> None:
    """Step 5-6: 筛选必需列并导出 Parquet + CSV + Augmented"""
    logger.info("=" * 70)

    # ---- 核心数据集（SEQTaRget 兼容，17列） ----
    for col in REQUIRED_COLUMNS:
        if col not in df.columns:
            logger.warning(f"  缺失必需列: {col}，填充 NaN")
            df[col] = np.nan

    df_core = df[REQUIRED_COLUMNS].copy()

    # 🆕 将所有 object/Decimal 类型列转换为合适的数值类型（避免 pyarrow 兼容性问题）
    import decimal
    for col in df_core.columns:
        if df_core[col].dtype == "object":
            # 检查列中的值类型
            sample_vals = df_core[col].dropna()
            if len(sample_vals) > 0:
                first_val = sample_vals.iloc[0]
                if isinstance(first_val, (decimal.Decimal, float, int)):
                    df_core[col] = pd.to_numeric(df_core[col], errors="coerce").astype("float64")
                    logger.info(f"  类型转换: {col} (object → float64)")
                else:
                    df_core[col] = df_core[col].astype(str)
                    logger.info(f"  类型转换: {col} (object → str)")

    pq_path = OUTPUT_DIR / "sequential_trials.parquet"
    table = pa.Table.from_pandas(df_core)
    pq.write_table(table, str(pq_path), compression="snappy")
    logger.info(
        f"核心 Parquet: {pq_path} ({pq_path.stat().st_size / 1024 / 1024:.1f} MB,"
        f" {len(df_core.columns)} 列)"
    )

    csv_path = OUTPUT_DIR / "sequential_trials.csv"
    df_core.to_csv(csv_path, index=False)
    logger.info(
        f"核心 CSV:     {csv_path} ({csv_path.stat().st_size / 1024 / 1024:.1f} MB)"
    )

    # 核心元数据
    meta = {
        "rows": len(df_core),
        "columns": len(df_core.columns),
        "column_names": list(df_core.columns),
        "build_time": datetime.now().isoformat(),
        "parquet_path": str(pq_path),
        "csv_path": str(csv_path),
        "treatment_distribution": (
            df_core["itt_group"].value_counts().to_dict()
            if "itt_group" in df_core.columns
            else {}
        ),
        "event_rate": (
            float(df_core["death_28day_num"].mean())
            if "death_28day_num" in df_core.columns
            else None
        ),
    }
    meta_path = OUTPUT_DIR / "sequential_trials_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    logger.info(f"核心元数据:  {meta_path}")

    # ---- 增强数据集（全协变量，39列，用于V26全协变量手动IPTW） ----
    missing_augmented = []
    for col in AUGMENTED_COLUMNS:
        if col not in df.columns:
            missing_augmented.append(col)

    if missing_augmented:
        logger.warning(
            f"增强数据集缺失列 ({len(missing_augmented)}个): {', '.join(missing_augmented)}"
        )
        logger.warning(
            f"  原因: 对应的数据库表不可用或合并失败。这些列不会被添加到增强数据集中。"
        )

    # 仅保留实际存在的 AUGMENTED_COLUMNS
    available_aug_cols = [c for c in AUGMENTED_COLUMNS if c in df.columns]
    df_aug = df[available_aug_cols].copy()

    # 🆕 增强数据集同样需要类型转换
    for col in df_aug.columns:
        if df_aug[col].dtype == "object":
            sample_vals = df_aug[col].dropna()
            if len(sample_vals) > 0:
                first_val = sample_vals.iloc[0]
                if isinstance(first_val, (decimal.Decimal, float, int)):
                    df_aug[col] = pd.to_numeric(df_aug[col], errors="coerce").astype("float64")
                else:
                    df_aug[col] = df_aug[col].astype(str)

    aug_pq_path = OUTPUT_DIR / "sequential_trials_augmented.parquet"
    aug_table = pa.Table.from_pandas(df_aug)
    pq.write_table(aug_table, str(aug_pq_path), compression="snappy")
    logger.info(
        f"增强 Parquet: {aug_pq_path} ({aug_pq_path.stat().st_size / 1024 / 1024:.1f} MB,"
        f" {len(df_aug.columns)} 列)"
    )

    aug_csv_path = OUTPUT_DIR / "sequential_trials_augmented.csv"
    df_aug.to_csv(aug_csv_path, index=False)
    logger.info(
        f"增强 CSV:     {aug_csv_path} ({aug_csv_path.stat().st_size / 1024 / 1024:.1f} MB)"
    )

    aug_meta = {
        "rows": len(df_aug),
        "columns": len(df_aug.columns),
        "column_names": list(df_aug.columns),
        "build_time": datetime.now().isoformat(),
        "parquet_path": str(aug_pq_path),
        "csv_path": str(aug_csv_path),
        "variable_groups": {
            "core": ["baseline_phosphate", "age", "gender", "baseline_sofa",
                     "egfr", "sepsis3", "invasive_vent_at_t0",
                     "vasoactive_at_t0", "charlson", "admission_type"],
            "medication": ["diuretic_at_t0", "steroid_at_t0", "insulin_at_t0",
                          "nahco3_infusion_at_t0", "beta2_agonist_at_t0", "ppi_use"],
            "acid_base": ["ph_worst", "bicarbonate_worst",
                         "alkalosis", "resp_alkalosis", "metabolic_alkalosis"],
            "procedure_transfusion": ["net_fluid_balance_24h", "transfusion_at_t0",
                                     "cardiac_surgery", "aki_kdigo2_plus"],
            "nutrition": ["refeeding_risk", "tpn_kcal_24h", "en_kcal_24h",
                         "nutrition_support_type",
                         "en_start_hours_from_t0", "tpn_start_hours_from_t0",
                         "gnri", "pni", "conut_score"],
            "renal_crrt": ["ever_crrt"],
            "calcium": ["corrected_calcium_mmol_min"],
        },
    }
    aug_meta_path = OUTPUT_DIR / "sequential_trials_augmented_metadata.json"
    with open(aug_meta_path, "w", encoding="utf-8") as f:
        json.dump(aug_meta, f, indent=2, ensure_ascii=False)
    logger.info(f"增强元数据:  {aug_meta_path}")


def main() -> None:
    start_time = datetime.now()
    logger.info("=" * 70)
    logger.info("v26 sequential_trials 数据集构建")
    logger.info(f"启动时间: {start_time.isoformat()}")
    logger.info(f"输出目录: {OUTPUT_DIR}")
    logger.info("=" * 70)

    config = load_config()
    pool = ConnectionPool(config)

    try:
        df = load_main_data(pool)
        df = merge_confounders(df, pool)
        df = clean_data(df)
        report_summary(df)
        export_data(df)
    finally:
        pool.close()

    elapsed = (datetime.now() - start_time).total_seconds()
    logger.info("=" * 70)
    logger.info(f"完成! 耗时 {elapsed:.0f}s")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()