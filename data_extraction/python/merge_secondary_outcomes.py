#!/usr/bin/env python3
"""
V26 次要结局变量合并脚本
将 outcome_data 表中的次要结局变量合并进 augmented parquet 文件

使用方式:
    python merge_secondary_outcomes.py

数据流:
    outcome_data (PostgreSQL) → 提取 → 合并 → sequential_trials_augmented.parquet (更新)
"""
import logging
import sys
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import yaml

# 项目根目录
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# 加载数据库配置
defaults_path = PROJECT_ROOT / "framework" / "config" / "defaults.yaml"
with open(defaults_path, "r", encoding="utf-8") as f:
    config_data = yaml.safe_load(f)

from framework.config.loader import Config
from framework.sql_extractor.connection import ConnectionPool

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")
logger = logging.getLogger("merge_secondary_outcomes")


# 次要结局变量定义（来自 V26_可用变量清单.md §4）
SECONDARY_OUTCOME_COLS = [
    "stay_id",
    "death_icu",
    "death_hospital",
    "death_90day",
    "vfd_28",
    "total_mv_days",
    "aki_stage2_plus",
    "phosphate_normalized",
    "phosphate_normalization_days",
    "icu_los_from_t0_days",
    "hospital_los_from_t0_days",
]


def extract_outcome_data(config: Config) -> pd.DataFrame:
    """从 PostgreSQL outcome_data 表提取次要结局变量"""
    logger.info("连接 PostgreSQL 并提取 outcome_data...")
    pool = ConnectionPool(config)
    sql = f"""
        SELECT {', '.join(SECONDARY_OUTCOME_COLS)}
        FROM outcome_data
    """
    df = pd.read_sql(sql, pool.conn)
    pool.close()
    logger.info(f"提取完成: {len(df)} 行, {len(df.columns)} 列")
    return df


def merge_and_save(augmented_path: Path, outcome_df: pd.DataFrame, dry_run: bool = False):
    """合并次要结局到 augmented parquet"""
    logger.info(f"读取 augmented parquet: {augmented_path}")
    table = pq.read_table(str(augmented_path))
    df_aug = table.to_pandas()

    # 统计合并前状态
    existing_cols = [c for c in SECONDARY_OUTCOME_COLS if c in df_aug.columns and c != "stay_id"]
    missing_cols = [c for c in SECONDARY_OUTCOME_COLS if c not in df_aug.columns and c != "stay_id"]
    logger.info(f"已存在次要结局: {existing_cols}")
    logger.info(f"缺失次要结局: {missing_cols}")

    if dry_run:
        logger.info("=== 干运行模式，不执行实际合并 ===")
        for col in missing_cols:
            if col in outcome_df.columns:
                logger.info(f"  将合并: {col} ({outcome_df[col].notna().sum()} 非空值)")
        return

    # 执行合并
    if not missing_cols:
        logger.info("所有次要结局变量已存在，跳过合并")
        return

    # 只合并缺失的列
    merge_cols = ["stay_id"] + missing_cols
    merge_subset = outcome_df[merge_cols].copy()
    # 去重: outcome_data 表可能有重复 stay_id，保留第一条
    n_before = len(merge_subset)
    merge_subset = merge_subset.drop_duplicates(subset="stay_id", keep="first")
    logger.info(f"去重: {n_before} → {len(merge_subset)} 行")
    df_aug = df_aug.merge(merge_subset, on="stay_id", how="left", suffixes=("", "_from_outcome"))
    logger.info(f"合并完成: {len(df_aug)} 行, {len(df_aug.columns)} 列")

    # 验证合并结果
    for col in missing_cols:
        if col in df_aug.columns:
            n_valid = df_aug[col].notna().sum()
            logger.info(f"  {col}: {n_valid}/{len(df_aug)} 非空值 ({100*n_valid/len(df_aug):.1f}%)")

    # 类型转换（确保与变量清单一致）
    type_map = {
        "death_icu": "int32",
        "death_hospital": "int32",
        "death_90day": "int32",
        "vfd_28": "float64",
        "total_mv_days": "float64",
        "aki_stage2_plus": "int32",
        "phosphate_normalized": "int32",
        "phosphate_normalization_days": "float64",
        "icu_los_from_t0_days": "float64",
        "hospital_los_from_t0_days": "float64",
    }
    for col, dtype in type_map.items():
        if col in df_aug.columns:
            try:
                df_aug[col] = df_aug[col].astype(dtype)
            except (ValueError, TypeError) as e:
                logger.warning(f"类型转换失败 {col} → {dtype}: {e}")

    # 保存
    backup_path = augmented_path.with_suffix(".parquet.bak")
    logger.info(f"备份原文件: {backup_path}")
    pq.write_table(table, str(backup_path))

    logger.info(f"保存合并后文件: {augmented_path}")
    new_table = pa.Table.from_pandas(df_aug)
    pq.write_table(new_table, str(augmented_path), compression="snappy")
    logger.info("合并完成！")


def main():
    logger.info("=" * 60)
    logger.info("V26 次要结局变量合并脚本")
    logger.info("=" * 60)

    # 路径配置
    version_dir = PROJECT_ROOT / "TTE" / "v26"
    augmented_path = version_dir / "output" / "sequential_trials_augmented.parquet"

    if not augmented_path.exists():
        logger.error(f"augmented parquet 不存在: {augmented_path}")
        sys.exit(1)

    # 加载配置
    config = Config(config_data)

    try:
        outcome_df = extract_outcome_data(config)
        merge_and_save(augmented_path, outcome_df)
    except Exception as e:
        logger.error(f"合并失败: {e}")
        logger.info("提示: 如果 PostgreSQL 不可用，可以使用已备份的 augmented parquet 文件")
        sys.exit(1)

    logger.info("=" * 60)
    logger.info("完成！后续分析脚本将可直接使用次要结局变量。")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()