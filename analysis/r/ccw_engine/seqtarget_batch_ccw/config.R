# ============================================================================
# SEQTaRget 分批扩展-加权-Cox 引擎配置文件
# 模块路径: TTE/v26/code/r/seqtarget_batch_ccw/
# 说明: 集中管理路径、变量、截断分位数、随机种子与批次参数
# ============================================================================

# ---- 加载项目级配置（路径、配色、日志函数） ----
# 若项目 config.R 不可用，则使用模块内默认路径
project_config <- file.path("..", "config.R")
if (file.exists(project_config)) {
  source(project_config, local = TRUE)
  USE_PROJECT_CONFIG <- TRUE
} else {
  USE_PROJECT_CONFIG <- FALSE
  # 自动检测项目根目录
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

# ---- 模块根目录与输出目录 ----
MODULE_DIR <- file.path(CODE_DIR, "seqtarget_batch_ccw")
OUTPUT_DIR <- file.path(DATA_DIR, "seqtarget_batch_ccw")
EXPANDED_DIR <- file.path(OUTPUT_DIR, "expanded")
WEIGHTED_DIR <- file.path(OUTPUT_DIR, "weighted_batches")
RESULTS_DIR <- file.path(OUTPUT_DIR, "results")
TEMP_DIR <- file.path(PROJECT_ROOT, "framework", "temp")

for (d in c(EXPANDED_DIR, WEIGHTED_DIR, RESULTS_DIR, TEMP_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 输入数据路径 ----
# 主分析使用 sequential_trials.parquet（19核心列）
INPUT_PARQUET <- file.path(DATA_DIR, "sequential_trials.parquet")
INPUT_PARQUET_AUGMENTED <- file.path(DATA_DIR, "sequential_trials_augmented.parquet")

# ---- 分析标识 ----
MODULE_PREFIX <- "seqtarget_batch_ccw"
RUN_TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE <- file.path(LOG_DIR, sprintf("%s_%s.log", MODULE_PREFIX, RUN_TIMESTAMP))

# ---- 核心变量名映射 ----
ID_COL <- "stay_id"
PERIOD_COL <- "trial_num"           # 数据库层 1-indexed 期号
OUTCOME_COL <- "death_28day_num"
SURV_TIME_COL <- "surv_time"
TREATMENT_SRC_COL <- "itt_group"
TREATMENT_LABEL <- "treatment"
CONTROL_LABEL <- "conservative"

# ---- PS 模型协变量集（DAG 12 核心变量，用于对比验证） ----
# 基于 DAG 后门路径分析确定的 12 个核心混杂变量
# 详见 TTE/v26/docs/v26_DAG.md
PS_COVARIATES <- c(
  "age",
  "baseline_sofa",
  "charlson",
  "sepsis3",
  "egfr",
  "ever_crrt",
  "invasive_vent_at_t0",
  "vasoactive_at_t0",
  "refeeding_risk",
  "lactate_t0",
  "liver_disease",
  "aki_kdigo2_plus"
)

# ---- IPTW 权重参数 ----
IPTW_TRUNC_UPPER <- 0.99   # 上侧截断分位数
IPTW_TRUNC_LOWER <- 0.01   # 下侧截断分位数（双侧截断）
IPTW_STABILIZE <- TRUE     # 稳定化权重
PS_CLIP <- c(1e-6, 1 - 1e-6)  # PS 裁剪区间，避免极端值

# ---- 随机种子 ----
RANDOM_SEED <- 20260619L

# ---- 批次参数 ----
# 本模块按 trial_period 分批，每个 period 为一个 batch
# BATCH_MAX_ROWS 可用于限制单个 batch 的内存占用（NULL 表示不限制）
BATCH_MAX_ROWS <- NULL

# ---- HR 汇总方法 ----
# "period_specific": 仅返回每期 HR 表格
# "fixed": 逆方差固定效应 Meta 分析
# "random": DerSimonian-Laird 随机效应 Meta 分析（需 meta 包）
POOL_METHOD <- "fixed"

# ---- SEQTaRget 对比参数 ----
COMPARE_SAMPLE_SIZES <- c(5000L)          # 验证使用的患者数（仅 n=5000 初步验证）
COMPARE_FULL_DATA <- TRUE                # 是否包含全量数据对比（约 22685 例）
COMPARE_BOOTSTRAP <- 100L               # SEQTaRget 对比时的 bootstrap 次数
COMPARE_NTHREADS <- 2L                  # SEQTaRget 并行线程数
COMPARE_COVARIATES <- PS_COVARIATES     # 对比使用相同协变量集（DAG 12 变量）

# ---- SEQTaRget 模型参数（本次对比试验专用） ----
COMPARE_FASTGLM_METHOD <- 1L            # fastglm 分解方法：1L = unpivoted QR
COMPARE_FOLLOWUP_INCLUDE <- TRUE        # 结局模型中包含 followup + followup²
COMPARE_TRIAL_INCLUDE <- TRUE           # 结局模型中包含 trial + trial²
COMPARE_SKIP_SEQTARGET <- TRUE          # 跳过 SEQTaRget（仅运行 batch 方法，提速）

# ---- 结构化日志函数 ----
# 优先使用项目 helper；否则使用简易 message 日志
if (exists("log_info") && is.function(log_info)) {
  module_log_info <- log_info
  module_log_warn <- if (exists("log_warn")) log_warn else function(msg) message("[WARN] ", msg)
  module_log_section <- if (exists("log_section")) log_section else function(msg) message("\n===== ", msg, " =====")
} else {
  module_log_info <- function(msg) cat(sprintf("[INFO] %s\n", msg))
  module_log_warn <- function(msg) cat(sprintf("[WARN] %s\n", msg))
  module_log_section <- function(msg) cat(sprintf("\n===== %s =====\n", msg))
}

module_log_section("SEQTaRget 分批引擎配置加载完成")
module_log_info(sprintf("项目根目录: %s", PROJECT_ROOT))
module_log_info(sprintf("输入文件: %s", INPUT_PARQUET))
module_log_info(sprintf("输出目录: %s", OUTPUT_DIR))
module_log_info(sprintf("PS 协变量: %s", paste(PS_COVARIATES, collapse = ", ")))
module_log_info(sprintf("HR 汇总方法: %s", POOL_METHOD))
