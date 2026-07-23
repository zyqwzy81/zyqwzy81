# ============================================================================
# SEQTaRget 兼容分批引擎 V2 — 统一配置文件
# 模块路径: TTE/v26/code/r/seqtarget_batch_v2/
# 
# 设计目标: 忠实复现 SEQTaRget 的 CCW 扩展逻辑（SEQexpand），
#           通过分批处理实现内存可控，逐期验证结果一致性。
#
# 与 V1 (seqtarget_batch_ccw) 的关键改进:
#   1. 扩展深度: patient-period → 日粒度展开（匹配 SEQexpand）
#   2. PS 模型: 全量数据 pooled PS（匹配 weight.preexpansion=TRUE）
#   3. 结局模型: pooled logistic + followup² + trial²（匹配 SEQTaRget ITT）
#   4. 验证机制: 逐期累加验证 HR 收敛性
# ============================================================================

# ---- 加载项目级配置 ----
project_config <- file.path("..", "config.R")
if (file.exists(project_config)) {
  source(project_config, local = TRUE)
  USE_PROJECT_CONFIG <- TRUE
} else {
  USE_PROJECT_CONFIG <- FALSE
  # 手动检测项目根目录
  get_project_root <- function() {
    env_root <- Sys.getenv("PHOSPHATE_PROJECT_ROOT")
    if (env_root != "") return(env_root)
    current_dir <- getwd()
    for (i in 1:8) {
      if (file.exists(file.path(current_dir, "CLAUDE.md")) ||
          file.exists(file.path(current_dir, "framework"))) {
        return(current_dir)
      }
      parent <- dirname(current_dir)
      if (parent == current_dir) break
      current_dir <- parent
    }
    return(getwd())
  }
  PROJECT_ROOT <- get_project_root()
  VERSION_DIR <- file.path(PROJECT_ROOT, "TTE", "v26")
  CODE_DIR <- file.path(VERSION_DIR, "code", "r")
  DATA_DIR <- file.path(VERSION_DIR, "output")
  FIG_DIR <- file.path(DATA_DIR, "figures")
  TBL_DIR <- file.path(DATA_DIR, "tables")
  RDS_DIR <- file.path(DATA_DIR, "rds")
  LOG_DIR <- file.path(DATA_DIR, "logs")
  for (d in c(FIG_DIR, TBL_DIR, RDS_DIR, LOG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 模块目录 ----
MODULE_DIR <- file.path(CODE_DIR, "seqtarget_batch_v2")
OUTPUT_DIR <- file.path(DATA_DIR, "seqtarget_batch_v2")
WEIGHT_DIR <- file.path(OUTPUT_DIR, "weighted")       # Phase 1: pooled PS 权重输出
EXPAND_DIR <- file.path(OUTPUT_DIR, "expanded")        # Phase 2: 日粒度展开输出
RESULT_DIR <- file.path(OUTPUT_DIR, "results")         # Phase 3-4: 分析结果 + 质量报告
TEMP_DIR   <- file.path(PROJECT_ROOT, "framework", "temp")

for (d in c(WEIGHT_DIR, EXPAND_DIR, RESULT_DIR, TEMP_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 输入数据 ----
INPUT_PARQUET          <- file.path(DATA_DIR, "sequential_trials.parquet")
INPUT_PARQUET_AUGMENTED <- file.path(DATA_DIR, "sequential_trials_augmented.parquet")

# ---- 运行标识 ----
MODULE_PREFIX <- "seqtarget_batch_v2"
RUN_TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE <- file.path(LOG_DIR, sprintf("%s_%s.log", MODULE_PREFIX, RUN_TIMESTAMP))

# ---- 核心变量名映射（与 config.R 对齐） ----
ID_COL         <- "stay_id"
PERIOD_COL     <- "trial_num"           # 1-indexed
OUTCOME_COL    <- "death_28day_num"
SURV_TIME_COL  <- "surv_time"
TREATMENT_COL  <- "itt_group"
TREAT_LABEL    <- "treatment"
CTRL_LABEL     <- "conservative"

# ---- PS 模型协变量集（DAG 12 核心变量） ----
# 与 V1 保持一致，以便对比
PS_COVARIATES <- c(
  "age", "baseline_sofa", "charlson", "sepsis3", "egfr",
  "ever_crrt", "invasive_vent_at_t0", "vasoactive_at_t0",
  "refeeding_risk", "lactate_t0", "liver_disease", "aki_kdigo2_plus"
)

# ---- 字符型协变量（PS 模型中需转为 factor）----
CHAR_COVARIATES <- c()  # DAG12 无字符型，保留供扩展

# ---- IPTW 参数（与 config.R 对齐） ----
IPTW_TRUNC_UPPER <- 0.99   # 上侧截断
IPTW_TRUNC_LOWER <- 0.01   # 下侧截断（双侧截断）
IPTW_STABILIZE   <- TRUE   # 稳定化权重
PS_CLIP          <- c(1e-6, 1 - 1e-6)  # PS 裁剪

# ---- SEQexpand 参数 ----
# SEQTaRget 的 SEQexpand: 对每 patient-period，从 day=0 到 floor(surv_time) 展开
# 死亡日 outcome=1，其余日 outcome=0
# 最长展开到 28 天（28天死亡结局）
SEQEXPAND_MAX_DAY <- 28L     # 最大展开天数
SEQEXPAND_ALLOW_ZERO_FU <- TRUE  # 允许 surv_time=0 的期（展开为 1 条 day=0 记录）

# ---- 处理组/对照组水平 ----
TX_LEVEL   <- 1L   # 治疗组编码
CTRL_LEVEL <- 0L   # 对照组编码

# ---- 分批策略 ----
# 分批方式: 按 trial_period 分批
# 注意: period 7 (n=3, events=0) 和 period 8 (n=2, events=0) 会被自动跳过
# 每批次的展开规模 ≈ n_patients × avg_surv_time
# Period 1 展开约 22685 × 9.65 ≈ 218,910 行（全量最重批次）
PERIOD_FILTER_MIN_EVENTS <- 5L   # 事件数 < 5 的 period 自动跳过
PERIOD_FILTER_MIN_N      <- 30L  # 样本数 < 30 的 period 自动跳过

# ---- Pooled logistic 参数（与 SEQTaRget ITT 对齐） ----
# logit(P) = β₀ + β₁*tx_init + β₂*followup + β₃*followup² + β₄*trial + β₅*trial² + covs
POOL_LOGISTIC_INCLUDE_FOLLOWUP_SQ <- TRUE  # followup²
POOL_LOGISTIC_INCLUDE_TRIAL_SQ    <- TRUE  # trial²
POOL_LOGISTIC_CLUSTER_BY          <- ID_COL  # cluster SE

# ---- Bootstrap 参数（可选） ----
BOOTSTRAP_N <- 0L  # 0 = 不执行（仅使用 robust SE），>0 时使用非参数 bootstrap

# ---- 随机种子 ----
RANDOM_SEED <- 20260619L

# ---- 结构化日志（与 V1 一致） ----
if (exists("log_info") && is.function(log_info)) {
  module_log_info    <- log_info
  module_log_warn    <- if (exists("log_warn")) log_warn else function(msg) message("[WARN] ", msg)
  module_log_section <- if (exists("log_section")) log_section else function(msg) message("\n===== ", msg, " =====")
} else {
  module_log_info    <- function(msg) cat(sprintf("[INFO] %s\n", msg))
  module_log_warn    <- function(msg) cat(sprintf("[WARN] %s\n", msg))
  module_log_section <- function(msg) cat(sprintf("\n===== %s =====\n", msg))
}

# ---- 峰值内存监控 ----
get_peak_mem_mb <- function() {
  if (file.exists("/proc/self/status")) {
    lines <- readLines("/proc/self/status", warn = FALSE)
    peak_line <- grep("^VmPeak:", lines, value = TRUE)
    if (length(peak_line) > 0) {
      kb <- as.numeric(gsub("[^0-9]", "", peak_line))
      return(round(kb / 1024, 2))
    }
  }
  NA_real_
}

# ---- 启动日志 ----
module_log_section("SEQTaRget 兼容分批引擎 V2 配置加载完成")
module_log_info(sprintf("项目根目录: %s", PROJECT_ROOT))
module_log_info(sprintf("输入文件: %s", INPUT_PARQUET))
module_log_info(sprintf("输出目录: %s", OUTPUT_DIR))
module_log_info(sprintf("PS 协变量 (%d): %s", length(PS_COVARIATES),
                         paste(PS_COVARIATES, collapse = ", ")))
module_log_info(sprintf("SEQexpand 最大展开天数: %d", SEQEXPAND_MAX_DAY))
module_log_info(sprintf("IPTW 截断: [%.2f, %.2f]", IPTW_TRUNC_LOWER, IPTW_TRUNC_UPPER))
