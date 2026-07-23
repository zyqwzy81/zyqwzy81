#!/usr/bin/env python3
"""检查 outcome_data 表的结构"""
import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from framework.config.loader import load_config
from framework.sql_extractor.connection import ConnectionPool

config = load_config(version="v26")
pool = ConnectionPool(config)

with pool.cursor() as cur:
    # 检查 outcome_data 是否存在
    cur.execute("SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'outcome_data')")
    exists = cur.fetchone()[0]
    print(f"outcome_data exists: {exists}")
    
    if exists:
        # 检查列
        cur.execute("SELECT column_name, data_type FROM information_schema.columns WHERE table_name = 'outcome_data' ORDER BY ordinal_position")
        cols = cur.fetchall()
        print(f"Columns ({len(cols)}):")
        for c in cols:
            print(f"  {c[0]}: {c[1]}")
        
        # 检查行数
        cur.execute("SELECT COUNT(*) FROM outcome_data")
        n = cur.fetchone()[0]
        print(f"Total rows: {n}")
        
        # 检查 death_28day 统计
        if any(c[0] == 'death_28day' for c in cols):
            cur.execute("SELECT death_28day, COUNT(*) FROM outcome_data GROUP BY death_28day")
            stats = cur.fetchall()
            print(f"death_28day distribution:")
            for s in stats:
                print(f"  {s[0]}: {s[1]}")
        
        # 检查与 study_cohort_washout24h 的匹配
        cur.execute("""
            SELECT COUNT(*) FROM study_cohort_washout24h sc
            LEFT JOIN outcome_data od ON sc.stay_id = od.stay_id AND sc.time_zero = od.time_zero
            WHERE od.stay_id IS NULL
        """)
        unmatched = cur.fetchone()[0]
        print(f"study_cohort_washout24h rows without outcome_data match: {unmatched}")

pool.close()