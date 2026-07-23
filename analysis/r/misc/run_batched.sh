#!/bin/bash
# =============================================================================
# V26.9 分批并行管线 — 受限于 16 Gi 内存
#
# Batch 1: SA1 (DAG13, boot=200, ~4Gi) + SA4 (Trunc, <0.5Gi) — 并行 ~2h
# Batch 2: SA2 (Grace, boot=50, ~2Gi) + SA3 (CCA, boot=50, ~2Gi) — 并行 ~6h
# Phase 2: 次要结局并行 (mclapply 3核, ~7h)
# Phase 3: CCW 诊断 + 发表级图表 (~1h)
# =============================================================================

set -o pipefail

PROJECT_ROOT="/data/小米云盘/文档/论文/血磷"
CODE_DIR="${PROJECT_ROOT}/TTE/v26/code/r"
LOG_DIR="${PROJECT_ROOT}/TTE/v26/output/v26_9/logs"
START_TIME=$(date +%s)

mkdir -p "${LOG_DIR}"

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
log_sec() {
    echo ""
    echo "══════════════════════════════════════════════════════════════"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "══════════════════════════════════════════════════════════════"
    echo ""
}
ok()   { log "[OK] $*"; }
err()  { log "[ERR] $*"; }
wait_proc() {
    local pid=$1 name=$2
    while kill -0 "$pid" 2>/dev/null; do sleep 30; done
    wait "$pid"
    local rc=$?
    if [ $rc -eq 0 ]; then ok "$name 完成"; else err "$name 失败 (exit=$rc)"; fi
    return $rc
}

cd "${CODE_DIR}" || { log "无法切换到 ${CODE_DIR}"; exit 1; }

# =============================================================================
# Batch 1: SA1 + SA4 并行
# =============================================================================
log_sec "Batch 1: SA1 (DAG13, boot=200) + SA4 (Trunc) 并行"

log "启动 SA1 (DAG13, 预计 ~2h)..."
Rscript --vanilla v26_9_sensitivity_dag13.R > "${LOG_DIR}/sa1_dag13_batch1.log" 2>&1 &
PID_SA1=$!
log "  PID=$PID_SA1"

log "启动 SA4 (Trunc, 预计 ~0.5h)..."
Rscript --vanilla v26_9_sensitivity_trunc.R > "${LOG_DIR}/sa4_trunc_batch1.log" 2>&1 &
PID_SA4=$!
log "  PID=$PID_SA4"

# 等待两者完成
wait_proc "$PID_SA1" "SA1 (DAG13)"
wait_proc "$PID_SA4" "SA4 (Trunc)"

BATCH1_END=$(date +%s)
log_sec "Batch 1 完成! 耗时: $(( (BATCH1_END - START_TIME) / 60 )) 分钟"

# =============================================================================
# Batch 2: SA2 + SA3 并行
# =============================================================================
log_sec "Batch 2: SA2 (Grace 24h, boot=50) + SA3 (CCA, boot=50) 并行"

BATCH2_START=$(date +%s)

log "启动 SA2 (Grace 24h, 预计 ~6h)..."
Rscript --vanilla v26_9_sensitivity_grace.R > "${LOG_DIR}/sa2_grace_batch2.log" 2>&1 &
PID_SA2=$!
log "  PID=$PID_SA2"

log "启动 SA3 (CCA, 预计 ~6h)..."
Rscript --vanilla v26_9_sensitivity_cca.R > "${LOG_DIR}/sa3_cca_batch2.log" 2>&1 &
PID_SA3=$!
log "  PID=$PID_SA3"

wait_proc "$PID_SA2" "SA2 (Grace)"
wait_proc "$PID_SA3" "SA3 (CCA)"

BATCH2_END=$(date +%s)
log_sec "Batch 2 完成! 耗时: $(( (BATCH2_END - BATCH2_START) / 60 )) 分钟"
SENS_ELAPSED=$(( (BATCH2_END - START_TIME) / 60 ))
log_sec "敏感性分析全部完成! 总耗时: ${SENS_ELAPSED} 分钟"

# =============================================================================
# Phase 2: 次要结局并行分析
# =============================================================================
log_sec "Phase 2: 次要结局并行分析 (mclapply 3核, 预计 ~7h)"

PH2_START=$(date +%s)
if Rscript --vanilla v26_9_secondary_parallel.R; then
    PH2_END=$(date +%s)
    ok "次要结局分析完成! 耗时: $(( (PH2_END - PH2_START) / 60 )) 分钟"
else
    err "次要结局分析失败! 管线继续..."
fi

# =============================================================================
# Phase 3: CCW 诊断 + 发表级图表
# =============================================================================
log_sec "Phase 3: CCW 诊断 + 发表级图表 (预计 ~1h)"

PH3_START=$(date +%s)

log "CCW 诊断..."
if Rscript --vanilla v26_9_ccw_diagnostics.R; then ok "CCW 诊断完成"; else err "CCW 诊断失败"; fi

log "发表级图表..."
if Rscript --vanilla v26_9_publication_figures.R; then ok "发表级图表完成"; else err "发表级图表失败"; fi

# =============================================================================
# 汇总
# =============================================================================
END_TIME=$(date +%s)
TOTAL=$(( END_TIME - START_TIME ))

echo ""
log_sec "══════════ V26.9 分批并行管线 — 全部完成 ══════════"
echo ""
echo "  总耗时:     $((TOTAL/3600))h $((TOTAL%3600/60))m"
echo "  Batch 1:    ~2h (SA1 DAG13 + SA4 Trunc)"
echo "  Batch 2:    ~6h (SA2 Grace + SA3 CCA)"
echo "  Phase 2:    ~7h (次要结局 3核并行)"
echo "  Phase 3:    ~1h (诊断+图表)"
echo ""
echo "  关键输出:"
echo "    主分析:   ${PROJECT_ROOT}/TTE/v26/output/v26_9/tables/v26_9_seqtarget_rubin_results.csv"
echo "              HR=0.9872 (0.8927-1.0916), p=0.8018 (DAG14, ×20 MICE)"
echo "    敏感性:   ${PROJECT_ROOT}/TTE/v26/output/v26_9/tables/v26_9_sensitivity_*.csv"
echo "    次要:     ${PROJECT_ROOT}/TTE/v26/output/v26_9/tables/v26_9_secondary_*.csv"
echo "    图表:     ${PROJECT_ROOT}/TTE/v26/output/v26_9/figures/"
echo "══════════════════════════════════════════════════════════════════════"
