#!/usr/bin/env python3
"""
运行 M0_07_outcome_extraction_optimized.sql 以填充 outcome_data 表
"""
import sys, time, logging
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from framework.config.loader import load_config
from framework.sql_extractor.connection import ConnectionPool

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("run_outcome_extraction")

def main():
    config = load_config(version="v26")
    pool = ConnectionPool(config)

    # 检查 outcome_data 当前状态
    with pool.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM outcome_data")
        n = cur.fetchone()[0]
        logger.info(f"outcome_data 当前行数: {n}")

    # 读取 SQL
    sql_path = PROJECT_ROOT / "sql" / "M0_07_outcome_extraction_optimized.sql"
    sql_content = sql_path.read_text(encoding="utf-8")
    logger.info(f"SQL 文件已读取: {sql_path} ({len(sql_content)} 字符)")

    # 执行 SQL
    t0 = time.time()
    try:
        with pool.cursor() as cur:
            cur.execute("SET statement_timeout = '14400s'")
            cur.execute(sql_content)
        logger.info(f"SQL 执行成功 ({time.time()-t0:.1f}s)")
    except Exception as e:
        logger.error(f"SQL 执行失败: {e}")
        pool.close()
        sys.exit(1)

    # 验证
    with pool.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM outcome_data")
        n = cur.fetchone()[0]
        cur.execute("SELECT death_28day, COUNT(*) FROM outcome_data GROUP BY death_28day")
        stats = cur.fetchall()
        logger.info(f"outcome_data 新行数: {n}")
        for s in stats:
            logger.info(f"  death_28day={s[0]}: {s[1]}")
    
    pool.close()
    logger.info("完成")

if __name__ == "__main__":
    main()