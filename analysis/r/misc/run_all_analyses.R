# =============================================================================
# V26.9 统一分析入口脚本
# 顺序执行所有核心分析，确保各步骤结果一致性
# =============================================================================
# 使用说明：
#   1. 确保 config.R 中参数已正确设置
#   2. 确保 MICE 插补数据已存在 (v26_mice_imputed_imp01-20.parquet)
#   3. 运行本脚本: Rscript v26_9_run_all_analyses.R
# =============================================================================

cat("========================================================\n")
cat("V26.9 统一分析入口 — 顺序执行所有分析\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

# ---- 阶段1: 加载配置 ----
cat("[阶段 1/5] 加载配置...\n")
source("config.R")

# ---- 阶段2: 主分析 ----
cat("\n[阶段 2/5] SEQTaRget ITT 主分析 (DAG14, ×20 MICE, boot=200)...\n")
cat("  预计运行时间: ~8-12小时\n")
source("v26_9_main_analysis.R")
cat("  ✅ 主分析完成\n")

# ---- 阶段3: 敏感性分析 ----
cat("\n[阶段 3/5] 敏感性分析...\n")

# SA1: DAG13 敏感性
cat("  SA1 DAG13 敏感性...\n")
source("v26_9_sensitivity_dag13.R")
cat("  ✅ DAG13 完成\n")

# SA13: 宽限期 24h
cat("  SA13 宽限期 24h...\n")
source("v26_9_sensitivity_grace.R")
cat("  ✅ Grace period 完成\n")

# CCA: 完整个案 LASSO 筛选
cat("  CCA 完整个案 LASSO...\n")
source("v26_9_sensitivity_cca.R")
cat("  ✅ CCA 完成\n")

# IPTW 截断敏感性
cat("  IPTW 截断敏感性...\n")
source("v26_9_sensitivity_trunc.R")
cat("  ✅ Truncation 完成\n")

# ---- 阶段4: 次要结局 ----
cat("\n[阶段 4/5] 次要结局分析...\n")
cat("  预计运行时间: ~4-6小时\n")
source("v26_9_secondary_unified_seqtarget.R")
cat("  ✅ 次要结局完成\n")

# ---- 阶段5: 诊断 + 图表 ----
cat("\n[阶段 5/5] 诊断与图表...\n")

cat("  CCW 诊断...\n")
source("v26_9_ccw_diagnostics.R")
cat("  ✅ CCW 诊断完成\n")

cat("  发表级图表...\n")
source("v26_9_publication_figures.R")
cat("  ✅ 图表生成完成\n")

cat("\n========================================================\n")
cat("V26.9 全部分析完成!\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")
