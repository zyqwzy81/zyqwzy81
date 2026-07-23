#!/bin/bash
# =============================================================================
# V26.9 完整 R 分析管线 — 后台运行脚本
# =============================================================================
# 用途: 按顺序执行所有分析步骤，含监控日志输出
# 使用方法:
#   nohup bash run_v26_9_pipeline.sh > /dev/null 2>&1 &
#   tail -f $(ls -t TTE/v26/output/v26_9/logs/pipeline_*.log | head -1)
# =============================================================================

set -o pipefail

# ---- 路径配置 ----
PROJECT_ROOT="/data/小米云盘/文档/论文/血磷"
CODE_DIR="${PROJECT_ROOT}/TTE/v26/code/r"
LOG_DIR="${PROJECT_ROOT}/TTE/v26/output/v26_9/logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
PIPELINE_LOG="${LOG_DIR}/pipeline_${TIMESTAMP}.log"
START_TIME=$(date +%s)

# ---- 阶段计数器 ----
STEP=0
TOTAL_STEPS=6

# 创建日志目录
mkdir -p "${LOG_DIR}"

# ---- 日志函数 ----
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*"
}
log_step() {
    STEP=$((STEP + 1))
    echo ""
    echo "======================================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [步骤 ${STEP}/${TOTAL_STEPS}] $*"
    echo "======================================================================"
    echo ""
}
log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*"
}
log_ok() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [OK] $*"
}

# 重定向所有输出到日志文件（同时保持终端显示）
exec > >(tee -a "${PIPELINE_LOG}") 2>&1

# =============================================================================
# 开始
# =============================================================================
echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║          V26.9 完整 R 分析管线 (12h ITT 定义)                       ║"
echo "╠══════════════════════════════════════════════════════════════════════╣"
echo "║  启动时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "║  日志文件: ${PIPELINE_LOG}"
echo "║  工作目录: ${CODE_DIR}"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

# 切换到 R 脚本目录
cd "${CODE_DIR}" || {
    log_error "无法切换到目录: ${CODE_DIR}"
    exit 1
}

# ---- 步骤 0: 重建 sequential_trials_imputed_v26.parquet ----
log_step "0/${TOTAL_STEPS}: 重建 sequential_trials_imputed_v26.parquet"

if Rscript v26_9_rebuild_imputed_parquet.R; then
    log_ok "imputed parquet 重建完成"
else
    log_error "imputed parquet 重建失败! 管线终止。"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ═══════ 管线异常终止 ═══════"
    exit 1
fi

# ---- 步骤 1: MICE + LASSO 变量筛选 ----
log_step "1/${TOTAL_STEPS}: MICE 多重插补 + MITS-LASSO 变量筛选"
log_info "预计运行时间: 2-4 小时 (m=20, maxit=20)"
log_info "数据源: sequential_trials_augmented.parquet (12h ITT)"

START_SUB=$(date +%s)
if Rscript _archive/v24_v25_legacy/v25_mice_lasso_variable_selection.R; then
    ELAPSED=$(( $(date +%s) - START_SUB ))
    log_ok "MICE + LASSO 完成 (耗时: $((ELAPSED/3600))h $((ELAPSED%3600/60))m)"
else
    log_error "MICE + LASSO 失败! 管线终止。"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ═══════ 管线异常终止 ═══════"
    exit 1
fi

# ---- 步骤 2: MICE → Parquet 格式转换 ----
log_step "2/${TOTAL_STEPS}: MICE 插补结果 → Parquet 格式"
log_info "预计运行时间: ~10 分钟"

START_SUB=$(date +%s)
if Rscript _archive/v26_legacy/v26_4_save_mice_parquet.R --force; then
    ELAPSED=$(( $(date +%s) - START_SUB ))
    log_ok "MICE Parquet 转换完成 (耗时: $((ELAPSED/60))m)"
    
    # 验证输出
    MICE_DIR="${PROJECT_ROOT}/sql/result/mice_imputed"
    echo ""
    echo "  [验证] MICE Parquet 文件:"
    ls -lh "${MICE_DIR}"/v26_mice_imputed_imp*.parquet 2>/dev/null | head -5
    echo "  ..."
    ls -lh "${MICE_DIR}"/v26_mice_imputed_imp*.parquet 2>/dev/null | tail -3
    echo "  [验证] 元数据:"
    cat "${MICE_DIR}/mice_metadata.json" 2>/dev/null | python3 -m json.tool 2>/dev/null | head -10
else
    log_error "MICE Parquet 转换失败! 管线终止。"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ═══════ 管线异常终止 ═══════"
    exit 1
fi

# ---- 步骤 3: SEQTaRget ITT 主分析 (×20 MICE, boot=200) ----
log_step "3/${TOTAL_STEPS}: SEQTaRget ITT 主分析 (DAG14, ×20 MICE, boot=200)"
log_info "预计运行时间: 8-12 小时"
log_info "每插补集 ~25-35 min, 共 20 集"
log_info "输出: v26_9/output/v26_9/tables/ + rds/ + figures/"

START_SUB=$(date +%s)
if Rscript v26_9_main_analysis.R; then
    ELAPSED=$(( $(date +%s) - START_SUB ))
    log_ok "SEQTaRget 主分析完成! (耗时: $((ELAPSED/3600))h $((ELAPSED%3600/60))m)"
    
    # 输出关键结果
    echo ""
    echo "  [结果] Rubin 合并结果:"
    cat "${PROJECT_ROOT}/TTE/v26/output/v26_9/tables/v26_9_seqtarget_rubin_results.csv" 2>/dev/null
    echo ""
else
    log_error "SEQTaRget 主分析失败! 管线终止。"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ═══════ 管线异常终止 ═══════"
    exit 1
fi

# ---- 步骤 4: 敏感性分析 ----
log_step "4/${TOTAL_STEPS}: 敏感性分析"
log_info "SA1: DAG13, SA2: 宽限期24h, SA3: CCA, SA4: IPTW截断"
log_info "预计运行时间: 2-4 小时"

START_SUB=$(date +%s)

# SA1: DAG13
log_info "[SA1/4] DAG13 敏感性分析..."
if Rscript v26_9_sensitivity_dag13.R; then log_ok "DAG13 完成"; else log_error "DAG13 失败!"; fi

# SA2: 宽限期 24h
log_info "[SA2/4] 宽限期 24h 敏感性分析..."
if Rscript v26_9_sensitivity_grace.R; then log_ok "Grace 24h 完成"; else log_error "Grace 24h 失败!"; fi

# SA3: CCA
log_info "[SA3/4] 完整个案 LASSO..."
if Rscript v26_9_sensitivity_cca.R; then log_ok "CCA 完成"; else log_error "CCA 失败!"; fi

# SA4: IPTW 截断
log_info "[SA4/4] IPTW 截断敏感性..."
if Rscript v26_9_sensitivity_trunc.R; then log_ok "Truncation 完成"; else log_error "Truncation 失败!"; fi

ELAPSED=$(( $(date +%s) - START_SUB ))
log_ok "敏感性分析完成 (耗时: $((ELAPSED/3600))h $((ELAPSED%3600/60))m)"

# ---- 步骤 5: 次要结局分析 ----
log_step "5/${TOTAL_STEPS}: 次要结局分析"
log_info "预计运行时间: 4-6 小时"

START_SUB=$(date +%s)
if Rscript v26_9_secondary_unified_seqtarget.R; then
    ELAPSED=$(( $(date +%s) - START_SUB ))
    log_ok "次要结局分析完成 (耗时: $((ELAPSED/3600))h $((ELAPSED%3600/60))m)"
else
    log_error "次要结局分析失败! 管线继续..."
fi

# ---- 步骤 6: CCW 诊断 + 发表级图表 ----
log_step "6/${TOTAL_STEPS}: CCW 诊断 + 发表级图表"
log_info "预计运行时间: ~1 小时"

START_SUB=$(date +%s)

# CCW 诊断
log_info "[1/2] CCW 诊断..."
if Rscript v26_9_ccw_diagnostics.R; then
    log_ok "CCW 诊断完成"
else
    log_error "CCW 诊断失败! 管线继续..."
fi

# 发表级图表
log_info "[2/2] 发表级图表..."
if Rscript v26_9_publication_figures.R; then
    log_ok "发表级图表完成"
else
    log_error "发表级图表失败! 管线继续..."
fi

ELAPSED=$(( $(date +%s) - START_SUB ))
log_ok "诊断+图表完成 (耗时: $((ELAPSED/60))m)"

# =============================================================================
# 完成汇总
# =============================================================================
END_TIME=$(date +%s)
TOTAL_ELAPSED=$(( END_TIME - START_TIME ))

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║          V26.9 完整 R 分析管线 — 完成 ✓                             ║"
echo "╠══════════════════════════════════════════════════════════════════════╣"
echo "║  完成时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "║  总耗时:   $((TOTAL_ELAPSED/3600))h $((TOTAL_ELAPSED%3600/60))m $((TOTAL_ELAPSED%60))s"
echo "║  日志文件: ${PIPELINE_LOG}"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""
echo "--- 管线各步骤汇总 ---"
echo "  步骤0: imputed parquet 重建"
echo "  步骤1: MICE 插补 + LASSO 筛选"
echo "  步骤2: MICE → Parquet 转换"
echo "  步骤3: SEQTaRget ITT ×20 Rubin 合并"
echo "  步骤4: 敏感性分析 (DAG13 / Grace24h / CCA / Trunc)"
echo "  步骤5: 次要结局分析"
echo "  步骤6: CCW 诊断 + 发表级图表"
echo ""
echo "--- 关键输出位置 ---"
echo "  主分析结果: TTE/v26/output/v26_9/tables/"
echo "  图表输出:   TTE/v26/output/v26_9/figures/"
echo "  RDS 对象:   TTE/v26/output/v26_9/rds/"
echo "=============================================================================="
