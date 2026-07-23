#!/bin/bash
# =============================================================================
# V26.9 并行敏感性分析启动器 — 每个 SA 独立 nohup
# 用法: bash run_v26_9_parallel_sas.sh [monitor]
#   不传参数 = 启动所有 4 个 SA 后台运行
#   monitor   = 监控已启动的进程（通过 PID 文件）
# =============================================================================

set -o pipefail

PROJECT_ROOT="/data/小米云盘/文档/论文/血磷"
CODE_DIR="${PROJECT_ROOT}/TTE/v26/code/r"
LOG_DIR="${PROJECT_ROOT}/TTE/v26/output/v26_9/logs"
PID_FILE="${LOG_DIR}/sa_pids.txt"

mkdir -p "${LOG_DIR}"

cd "${CODE_DIR}" || exit 1

if [ "$1" = "monitor" ]; then
    # ===== 监控模式 =====
    if [ ! -f "${PID_FILE}" ]; then
        echo "错误: PID 文件 ${PID_FILE} 不存在"
        exit 1
    fi
    
    echo "开始监控 Phase 1 (4 项敏感性分析)..."
    echo "PID 文件: ${PID_FILE}"
    echo "日志目录: ${LOG_DIR}"
    echo ""
    
    # 读取 PID
    PID_SA1=$(grep '^SA1=' "${PID_FILE}" | cut -d= -f2)
    PID_SA2=$(grep '^SA2=' "${PID_FILE}" | cut -d= -f2)
    PID_SA3=$(grep '^SA3=' "${PID_FILE}" | cut -d= -f2)
    PID_SA4=$(grep '^SA4=' "${PID_FILE}" | cut -d= -f2)
    
    START_TIME=$(date +%s)
    
    while true; do
        sleep 30
        NOW=$(date '+%H:%M:%S')
        
        kill -0 $PID_SA1 2>/dev/null; S1=$?
        kill -0 $PID_SA2 2>/dev/null; S2=$?
        kill -0 $PID_SA3 2>/dev/null; S3=$?
        kill -0 $PID_SA4 2>/dev/null; S4=$?
        
        ELAPSED=$(( $(date +%s) - START_TIME ))
        ELAPSED_STR="$((ELAPSED/3600))h $((ELAPSED%3600/60))m"
        
        SA1_STAT=$([ $S1 -eq 0 ] && echo "运行中" || echo "已完成")
        SA2_STAT=$([ $S2 -eq 0 ] && echo "运行中" || echo "已完成")
        SA3_STAT=$([ $S3 -eq 0 ] && echo "运行中" || echo "已完成")
        SA4_STAT=$([ $S4 -eq 0 ] && echo "运行中" || echo "已完成")
        
        echo "[$NOW] 已运行 ${ELAPSED_STR} | SA1=${SA1_STAT} | SA2=${SA2_STAT} | SA3=${SA3_STAT} | SA4=${SA4_STAT}"
        
        # 全部完成则退出
        [ $S1 -ne 0 ] && [ $S2 -ne 0 ] && [ $S3 -ne 0 ] && [ $S4 -ne 0 ] && break
    done
    
    END_TIME=$(date +%s)
    TOTAL=$(( END_TIME - START_TIME ))
    echo ""
    echo "========================================"
    echo "Phase 1 全部完成! 总耗时: $((TOTAL/3600))h $((TOTAL%3600/60))m"
    echo "========================================"
    echo ""
    echo "日志文件:"
    echo "  SA1 (DAG13):  ${LOG_DIR}/sa1_dag13_*.log"
    echo "  SA2 (Grace):  ${LOG_DIR}/sa2_grace_*.log"
    echo "  SA3 (CCA):    ${LOG_DIR}/sa3_cca_*.log"
    echo "  SA4 (Trunc):  ${LOG_DIR}/sa4_trunc_*.log"
    exit 0
fi

# ===== 启动模式 =====
echo "启动 4 项敏感性分析并行运行..."
echo ""

# 清理旧日志
rm -f "${LOG_DIR}"/sa{1,2,3,4}_*.log 2>/dev/null

# SA1: DAG13 (boot=200, ~2h)
LOG1="${LOG_DIR}/sa1_dag13_$$.log"
nohup Rscript --vanilla v26_9_sensitivity_dag13.R > "${LOG1}" 2>&1 &
PID_SA1=$!
echo "SA1 (DAG13): PID=${PID_SA1}  log=${LOG1}"

# SA2: Grace 24h (boot=50, ~6h)
LOG2="${LOG_DIR}/sa2_grace_$$.log"
nohup Rscript --vanilla v26_9_sensitivity_grace.R > "${LOG2}" 2>&1 &
PID_SA2=$!
echo "SA2 (Grace): PID=${PID_SA2}  log=${LOG2}"

# SA3: CCA (boot=50, ~6h)
LOG3="${LOG_DIR}/sa3_cca_$$.log"
nohup Rscript --vanilla v26_9_sensitivity_cca.R > "${LOG3}" 2>&1 &
PID_SA3=$!
echo "SA3 (CCA):   PID=${PID_SA3}  log=${LOG3}"

# SA4: Trunc (~0.5h)
LOG4="${LOG_DIR}/sa4_trunc_$$.log"
nohup Rscript --vanilla v26_9_sensitivity_trunc.R > "${LOG4}" 2>&1 &
PID_SA4=$!
echo "SA4 (Trunc): PID=${PID_SA4}  log=${LOG4}"

# 保存 PID
cat > "${PID_FILE}" << EOF
SA1=${PID_SA1}
SA2=${PID_SA2}
SA3=${PID_SA3}
SA4=${PID_SA4}
EOF

echo ""
echo "PID 已保存至: ${PID_FILE}"
echo ""
echo "监控命令:"
echo "  bash run_v26_9_parallel_sas.sh monitor"
echo "  tail -f ${LOG1} ${LOG2} ${LOG3} ${LOG4}"
