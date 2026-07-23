#!/usr/bin/env python3
"""
运行 24h 洗脱期队列构建 SQL 脚本
读取 M0_02_cohort_building_washout24h.sql 并执行
生成 study_cohort_washout24h 表
"""
import sys, time, logging
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from framework.config.loader import load_config
from framework.sql_extractor.connection import ConnectionPool

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("run_sql_washout24h")

def main():
    # 加载配置 (默认连接 mimiciv 数据库)
    config = load_config(version="v26")
    logger.info(f"Config loaded: dbname={config.database.dbname}")

    # 创建连接池
    pool = ConnectionPool(config)
    
    # 1. 检查 core_cohort 是否存在
    with pool.cursor() as cur:
        cur.execute("SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'core_cohort')")
        exists = cur.fetchone()[0]
        logger.info(f"core_cohort 表存在: {exists}")
    
    if not exists:
        logger.error("core_cohort 表不存在，请先运行 00_core_cohort_V2.sql")
        pool.close()
        sys.exit(1)

    # 2. 读取 SQL 文件
    sql_path = PROJECT_ROOT / "sql" / "M0_02_cohort_building_washout24h.sql"
    if not sql_path.exists():
        logger.error(f"SQL 文件不存在: {sql_path}")
        pool.close()
        sys.exit(1)
    
    sql_content = sql_path.read_text(encoding="utf-8")
    logger.info(f"SQL 文件已读取: {sql_path} ({len(sql_content)} 字符)")

    # 3. 执行 SQL (逐条执行，避免多语句事务问题)
    # 按分号分割，但保留 BEGIN/COMMIT
    start_time = time.time()
    
    try:
        with pool.cursor() as cur:
            # 先设置超时 (4小时)
            cur.execute("SET statement_timeout = '14400s'")
            # 执行整个 SQL
            cur.execute(sql_content)
        logger.info("SQL 执行成功")
    except Exception as e:
        logger.error(f"SQL 执行失败: {e}")
        pool.close()
        sys.exit(1)

    elapsed = time.time() - start_time
    logger.info(f"SQL 执行耗时: {elapsed:.1f} 秒")

    # 4. 验证新表
    with pool.cursor() as cur:
        # 检查表是否存在
        cur.execute("SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'study_cohort_washout24h')")
        exists = cur.fetchone()[0]
        logger.info(f"study_cohort_washout24h 表存在: {exists}")
        
        if exists:
            # 统计行数
            cur.execute("SELECT COUNT(*) FROM study_cohort_washout24h")
            total = cur.fetchone()[0]
            cur.execute("SELECT COUNT(*) FROM study_cohort_washout24h WHERE exclusion_reason IS NOT NULL")
            excluded = cur.fetchone()[0]
            cur.execute("SELECT COUNT(DISTINCT stay_id) FROM study_cohort_washout24h")
            unique_stay = cur.fetchone()[0]
            
            logger.info("=" * 50)
            logger.info("study_cohort_washout24h 统计 (24h洗脱期)")
            logger.info(f"总行数: {total}")
            logger.info(f"排除行数: {excluded}")
            logger.info(f"有效行数: {total - excluded}")
            logger.info(f"唯一stay_id: {unique_stay}")
            logger.info(f"排除率: {100.0 * excluded / total:.1f}%")
            logger.info("=" * 50)

    pool.close()
    logger.info("完成")

if __name__ == "__main__":
    main()