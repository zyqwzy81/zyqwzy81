#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V28.4 次要结局"事件时间序列"提取（PostgreSQL -> Parquet）

为多次入组（multi-enrollment）目标试验提供两项原始日级事件时间：
  1) KDIGO AKI 分期序列（mimiciv_derived.kdigo_stages）
       字段: stay_id, charttime, aki_stage
       用途: 每个试验起点后"首次 aki_stage>=2"的事件日
  2) 红细胞输注事件（mimiciv_icu.inputevents, itemid 225168/226368 = PRBC）
       字段: stay_id, starttime
       用途: 每个试验起点后"首次 PRBC 输注"的事件日

仅提取研究队列（public.study_cohort）涉及的 stay_id，避免全库扫描。
口令不写入脚本：从环境变量 MIMIC_PG_PASSWORD 读取；未设置则走 .pgpass/信任连接。

输出（V28.4 分析专属，不改动 sql/result 核心索引）:
  TTE/v28.4/output/extraction/kdigo_stages_ts.parquet
  TTE/v28.4/output/extraction/prbc_events.parquet

用法:
  MIMIC_PG_PASSWORD=xxx python extract_event_times.py
"""

import os
import sys
from pathlib import Path

import pandas as pd
import psycopg

# ---- 路径 ----
PROJECT_ROOT = Path(__file__).resolve().parents[4]
OUT_DIR = PROJECT_ROOT / "TTE" / "v28.4" / "output" / "extraction"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---- 连接参数（口令仅环境变量）----
PGPW = os.environ.get("MIMIC_PG_PASSWORD", "")
DSN = (
    "host=localhost port=5432 dbname=mimiciv user=postgres"
    + (f" password={PGPW}" if PGPW else "")
)


def run_query(cur, sql, label):
    print(f"[extract] {label} ...", flush=True)
    cur.execute(sql)
    cols = [d[0] for d in cur.description]
    rows = cur.fetchall()
    df = pd.DataFrame(rows, columns=cols)
    print(f"[extract] {label}: {len(df):,} rows", flush=True)
    return df


def main():
    with psycopg.connect(DSN) as conn:
        # KDIGO 分期序列（仅队列 stay）
        sql_kdigo = """
            SELECT k.stay_id,
                   k.charttime,
                   k.aki_stage_creat,
                   k.aki_stage_uo,
                   k.aki_stage_crrt,
                   k.aki_stage
            FROM mimiciv_derived.kdigo_stages k
            WHERE k.stay_id IN (SELECT DISTINCT stay_id FROM public.study_cohort)
              AND k.charttime IS NOT NULL
        """
        # PRBC 输注事件（仅队列 stay）
        sql_prbc = """
            SELECT ie.stay_id,
                   ie.starttime
            FROM mimiciv_icu.inputevents ie
            WHERE ie.stay_id IN (SELECT DISTINCT stay_id FROM public.study_cohort)
              AND ie.itemid IN (225168, 226368)
              AND ie.starttime IS NOT NULL
        """
        with conn.cursor() as cur:
            dk = run_query(cur, sql_kdigo, "kdigo_stages")
            dp = run_query(cur, sql_prbc, "prbc_events")

    # 类型规整
    dk["stay_id"] = dk["stay_id"].astype("int64")
    dk["charttime"] = pd.to_datetime(dk["charttime"])
    # 分期分量可能为 NULL：保留为可空 Int64（pyarrow 支持）
    for c in ("aki_stage_creat", "aki_stage_uo", "aki_stage_crrt", "aki_stage"):
        dk[c] = pd.to_numeric(dk[c], errors="coerce").astype("Int64")
    dp["stay_id"] = dp["stay_id"].astype("int64")
    dp["starttime"] = pd.to_datetime(dp["starttime"])

    kp = OUT_DIR / "kdigo_stages_ts.parquet"
    pp = OUT_DIR / "prbc_events.parquet"
    dk.to_parquet(kp, index=False, compression="snappy")
    dp.to_parquet(pp, index=False, compression="snappy")
    print(f"[extract] saved {kp}", flush=True)
    print(f"[extract] saved {pp}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
