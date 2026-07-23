# ============================================================================
# V26.4 统一配置文件
# 所有分析脚本通过 source("config.R") 获取路径和参数配置
# V26.4: 51-covariate pool + sofa_t0 + ecg_abnormal - crp_t0(不可用) + BOOT=1000
# ============================================================================

# ---- 项目根目录（自动检测） ----
# 从当前脚本所在位置向上4级找到项目根目录
# 当前: TTE/v26/code/r/config.R → 项目根目录
get_project_root <- function() {
  # 优先使用环境变量
  env_root <- Sys.getenv("PHOSPHATE_PROJECT_ROOT")
  if (env_root != "") return(env_root)
  # 自动检测脚本路径
  current_dir <- getwd()
  # 方法1：使用命令行 --file= 参数（Rscript 方式）
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    current_dir <- normalizePath(dirname(sub("^--file=", "", file_arg)))
  }
  # 方法2：使用 sys.frame 获取 source 调用的脚本路径（仅交互模式）
  script_path <- try(sys.frame(1)$ofile, silent = TRUE)
  if (!inherits(script_path, "try-error") && !is.null(script_path)) {
    current_dir <- normalizePath(dirname(script_path))
  }
  # 向上遍历查找项目根目录（含 CLAUDE.md 或 framework 目录）
  for (i in 1:8) {
    if (file.exists(file.path(current_dir, "CLAUDE.md")) ||
        file.exists(file.path(current_dir, "framework"))) {
      return(current_dir)
    }
    parent <- dirname(current_dir)
    if (parent == current_dir) break  # 到达文件系统根
    current_dir <- parent
  }
  return(getwd())
}

PROJECT_ROOT <- get_project_root()

# ---- 版本目录 ----
VERSION_DIR <- file.path(PROJECT_ROOT, "TTE", "v26")
CODE_DIR    <- file.path(VERSION_DIR, "code", "r")
DATA_DIR    <- file.path(VERSION_DIR, "output")
FIG_DIR     <- file.path(DATA_DIR, "figures")
TBL_DIR     <- file.path(DATA_DIR, "tables")
RDS_DIR     <- file.path(DATA_DIR, "rds")
LOG_DIR     <- file.path(DATA_DIR, "logs")

# 确保输出目录存在
for (d in c(FIG_DIR, TBL_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 数据文件路径（主输入） ----
# 注：这些是默认值，各脚本可通过命令行参数覆盖
DATA_PRIMARY   <- file.path(DATA_DIR, "sequential_trials.parquet")
DATA_AUGMENTED <- file.path(DATA_DIR, "sequential_trials_augmented.csv")  # V26.4: CSV格式(86列)
DATA_AUGMENTED_FULL <- file.path(DATA_DIR, "sequential_trials_augmented_full.parquet")
DATA_IMPUTED   <- file.path(RDS_DIR, "v26_mice_imputed.rds")  # MICE 插补对象（RDS 格式，保留 mids 结构）

# ---- 🆕 共享 MICE Parquet 路径（多版本通用） ----
# 保存至 sql/result/mice_imputed/，供 V26/V27 共享
MICE_PARQUET_DIR <- file.path(PROJECT_ROOT, "sql", "result", "mice_imputed")
MICE_PARQUET_ALL   <- file.path(MICE_PARQUET_DIR, "v26_mice_imputed_all.parquet")     # 全量（含 .imp 列）
MICE_PARQUET_SINGLE <- file.path(MICE_PARQUET_DIR, "v26_mice_imputed_imp%02d.parquet")  # 单插补集模板
MICE_PARQUET_META   <- file.path(MICE_PARQUET_DIR, "mice_metadata.json")               # 元数据

# ---- 🆕 共享核心数据路径（多版本通用） ----
# 由 build_sequential_trials.py 生成，复制至 sql/result/ 供 V26/V27 共享
SHARED_RESULT_DIR <- file.path(PROJECT_ROOT, "sql", "result")
SHARED_SEQ_TRIALS       <- file.path(SHARED_RESULT_DIR, "sequential_trials.parquet")
SHARED_SEQ_TRIALS_AUG   <- file.path(SHARED_RESULT_DIR, "sequential_trials_augmented.parquet")
SHARED_SEQ_TRIALS_META  <- file.path(SHARED_RESULT_DIR, "sequential_trials_metadata.json")

# ---- 协变量管理（V26.4 修订版，2026-06-18） ----
# 原则：协变量集在框架搭建阶段定义为可配置列表，
#       最终协变量数量在框架确定后根据分析需求填充。
#
# COVARIATE_POOL 是 V26 所有分析脚本中协变量的**唯一权威来源**。
# 最终协变量集 = COVARIATE_POOL ∩ 数据现有列 ∩ LASSO筛选结果(≥80%选择率)。

# V26.4 候选池状态:
#   COVARIATE_POOL: 50 变量
#   数据中可用: 49 变量 (98.0%)
#   数据中缺失: 1 变量 — crp_t0 (MIMIC-IV 匹配率 <7%, 不可用)
#   高缺失警告(>50%): bilirubin_t0(54.8%), albumin_min(51.4%) — 保留但标注
#   新增可用: troponin_elevated_pre (M0_36), tpn_at_t0, en_at_t0, malnutrition_at_t0
#   移除: crp_t0 (覆盖率<7%, 不可修复)

# 全协变量集（候选池）— 所有可用协变量的最大集合
COVARIATE_POOL <- c(
  # 人口学 (5) — 全部可用, 0%缺失
  "age", "gender", "admission_type", "bmi", "race_group",
  # 疾病严重度 (9) — V26.4: 新增 sofa_t0
  "baseline_sofa", "sofa_t0", "sapsii", "charlson", "sepsis3", "aki_kdigo2_plus",
  "ever_crrt", "cardiac_surgery", "liver_disease",
  # 实验室 (5) — V26.4: 移除 crp_t0
  "baseline_phosphate", "egfr", "ph_worst", "bicarbonate_worst",
  "corrected_calcium_mmol_min", "magnesium_t0",
  # 电解质 (2) — 全部可用, <1%缺失
  "potassium_t0", "sodium_t0",
  # 肝功能 (2) — 高缺失警告: bilirubin(54.8%), albumin(51.4%)
  "bilirubin_t0", "albumin_min",
  # 肾功能 (1)
  "bun_t0",
  # 血常规 (3) — 全部可用, <3%缺失
  "hemoglobin_min", "platelets_min", "wbc_max",
  # 凝血 (2) — 可用, ~14%缺失
  "inr_max", "ptt_max",
  # 炎症/灌注 (1) — 仅 lactate(31.3%), crp_t0 不可用
  "lactate_t0",
  # 神经系统 (1) — 29.3%缺失
  "gcs",
  # 心肌状态 (2) — ecg_abnormal 已可用(0%), troponin 仍不可用
  "troponin_elevated_pre", "ecg_abnormal",
  # 呼吸 (2) — 全部可用
  "invasive_vent_at_t0", "sf_ratio",
  # 治疗/用药 (8) — 全部可用, 0%缺失
  "vasoactive_at_t0", "steroid_at_t0", "insulin_at_t0",
  "nahco3_infusion_at_t0", "diuretic_at_t0", "beta2_agonist_at_t0",
  "transfusion_at_t0", "ppi_use",
  # 营养/代谢 (5) — refeeding_risk 已重新定义(NICE NG21), tpn_at_t0/en_at_t0/malnutrition 不可用
  "tpn_at_t0", "en_at_t0", "malnutrition_at_t0", "refeeding_risk",
  "net_fluid_balance_24h",
  # ICU类型 (1)
  "icu_type"
)

# 数据中不可用的 COVARIATE_POOL 变量（自动排除）
COVARIATE_UNAVAILABLE <- c("crp_t0")

# 高缺失率警告变量 (>50%, 保留但标注)
COVARIATE_HIGH_MISSING <- c("bilirubin_t0", "albumin_min")
# 注: ppi_use、nutrition_support_type 等字符型变量在 MICE 时需要特殊处理
#      通过 v26_mice_lasso_variable_selection.R 的 exclude_for_mice 自动排除

# LASSO 筛选排除的字符型变量（自动保留但不参与 LASSO 惩罚选择）
LASSO_CHAR_EXCLUDE <- c("gender", "admission_type", "race_group", "icu_type")

# 时变协变量（time-varying covariates）
# 这些变量在 CCCW 框架中随 trial period 变化
COVARIATE_TIME_VARYING <- c(
  "baseline_phosphate"
)

# 默认主分析协变量集（全量，通过 intersect 过滤实际存在的列）
# 使用时: covariates <- get_covariates(data, COVARIATE_POOL)
get_covariates <- function(data, pool = COVARIATE_POOL) {
  intersect(pool, names(data))
}

# ---- 异常记录排除阈值配置 ----
ANOMALOUS_EXCLUSION_RULES <- list(
  egfr_threshold      = 15,    # eGFR < 15 视为极端异常
  short_stay_hours    = 24,    # ICU 住院 < 24h 视为 short stay
  extreme_phosphate_max = 5.0, # 血磷 > 5.0 mmol/L 视为极端异常
  extreme_age_max     = 120    # 年龄 > 120 视为异常
)

# ---- 协变量预处理配置 ----
MISSING_EXCLUDE_PCT <- 80  # 变量缺失率超过该阈值则排除（%）
LASSO_DISPLAY_WARN_PCT <- 50  # LASSO 选择率显示警告阈值（%）
CONFIDENCE_LEVEL <- 0.975            # 95% CI 置信水平（双侧）
Z_ALPHA_95 <- qnorm(CONFIDENCE_LEVEL) # 95% CI 标准正态分位数
ALPHA_LEVEL <- 0.05                  # 全局显著性水平

# 置信区间分位数（Rubin 合并等小样本场景使用 t 分布）
get_ci_quantile <- function(df = Inf) {
  qt(CONFIDENCE_LEVEL, df)
}

# ---- CCW vs IPTW 一致性评估阈值 ----
CONSISTENCY_HR_DISAGREE_THRESHOLD <- 0.15  # HR 相对差异 > 15% 视为方向性分歧
CONSISTENCY_NULL_TOLERANCE <- 0.05         # HR 与 1 的距离小于该值视为 "无效应"
CONSISTENCY_PS_OVERLAP_MIN <- 0.10         # PS 重叠区域下限
CONSISTENCY_PS_OVERLAP_MAX <- 0.90         # PS 重叠区域上限

# ---- MaxCombo 非比例风险检验配置 ----
MAXCOMBO_ENABLED <- TRUE
MAXCOMBO_NPERM   <- 2000L  # 置换次数
MAXCOMBO_SCORES  <- c("logrank", "gehan", "peto", "fh01")

# ---- 偏差量化配置 ----
BIAS_QUANTIFICATION_ENABLED <- TRUE
BIAS_ENROLLMENT_CORRELATION <- TRUE  # 入组次数-死亡率相关性
BIAS_MAX_ENROLLMENT_LEVELS  <- c(3L, 5L)  # 限制入组次数进行敏感性分析
BIAS_GRACE_HR_DIFF_QUANTIFY <- TRUE  # 量化宽限期 HR 差异

# ---- 阴性对照分析配置 ----
NEGATIVE_CONTROL_ENABLED <- TRUE
NEGATIVE_CONTROL_OUTCOMES <- c("non_cv_death_28day_num")  # 需数据中存在

# ---- Delta adjustment (MNAR) 配置 ----
DELTA_ADJUSTMENT_ENABLED <- TRUE
DELTA_ADJUSTMENT_VALUES  <- c(0.5, 1.0, 1.5, 2.0)

# ---- Weibull AFT 配置 ----
AFT_ENABLED <- TRUE
AFT_DIST    <- "weibull"  # 可选: weibull | loglogistic | lognormal

# ---- 输出文件命名配置 ----
OUTPUT_PREFIX <- "v26_4"

# ---- MICE 插补参数（统一配置） ----
MICE_M     <- 20L    # 插补数据集数量
MICE_MAXIT <- 20L    # 最大迭代次数
MICE_SEED  <- 12345L    # 实际MICE运行使用的种子（论文已同步声明为此值）

# ---- IPTW 权重参数（统一配置） ----
IPTW_TRUNC_PCT      <- 0.99   # 上侧截断分位数
IPTW_TRUNC_LOWER    <- 0.01   # 🆕 下侧截断分位数（双侧截断，Cole & Hernán 2008）
IPTW_STABILIZE      <- TRUE   # 启用稳定化权重
IPTW_BILATERAL      <- TRUE   # 🆕 启用双侧截断（否则仅上侧）

# ---- logistf Firth 逻辑回归配置 ----
LOGISTF_MAXIT <- 100L  # Firth 回归最大迭代次数

# ---- Bootstrap 参数 ----
BOOT_MAIN        <- 200L   # 主分析 bootstrap 次数 (SEQTaRget推荐≥200)
BOOT_SENSITIVITY <- 50L    # 敏感性分析 bootstrap 次数 (SEQTaRget推荐≥50)
BOOT_KM          <- 100L   # KM 曲线 bootstrap 次数

# 测试子集模式（NULL = 全量数据；整数 = 按 stay_id 抽取指定数量）
TEST_SUBSET_N    <- 5000L  # [测试专用] 全流程子集测试样本量

# ---- 🆕 LASSO 变量筛选配置（V26.4 修订） ----
# 方法: 在 MICE 插补后对每个插补集独立运行 LASSO (λ.1se, 10-fold CV)
#       PS模型(二分类) + 结局模型(Cox) 双模型取并集
#       选择率阈值: ≥80% (在 m=20 个插补集中被选中的比例)
#       字符型变量(gender/admission_type/race_group/icu_type) 不参与LASSO但自动保留
LASSO_FILL_METHOD <- "colmean"  # 插补前填充方法
LASSO_SELECTION_THRESHOLD <- 0.80  # 选择率阈值 (V26.4: 0.50→0.80, 更严格)
LASSO_NFOLDS <- 10L  # 交叉验证折数
LASSO_LAMBDA_RULE <- "lambda.1se"  # 惩罚参数选择规则
LASSO_DUAL_MODEL <- "union"  # PS+结局模型合并方式: "union" | "intersect"

# ---- 🆕 敏感性分析完整配置（V26.4 修订） ----
# 敏感性分析矩阵: 方法 × 参数组合

# 1. 宽限期敏感性
GRACE_PERIODS <- c(12, 24, 48)  # 小时 (12h=主分析)

# 2. IPTW 截断敏感性
TRUNCATION_LEVELS <- c(0.95, 0.975, 0.99, 0.995)  # 上侧截断分位数

# 3. 协变量集敏感性（数值阈值，用于 LASSO 百分比筛选）
COVARIATE_SETS <- list(
  C1 = 0.80,   # 主分析: LASSO筛选 (≥80%)
  C2 = 0.50,   # 敏感性: 更多变量 (≥50%)
  C3 = 0.90,   # 敏感性: 更少变量 (≥90%)
  C4 = 0.90    # 敏感性: 双稳健变量近似 (≥90% 代理交集)
)

# 4. 权重方法敏感性
WEIGHT_METHODS <- c("iptw", "overlap", "matching")
WEIGHT_PRIMARY <- "iptw"  # 主权重方法

# 5. 缺失数据处理敏感性
MISSING_DATA_METHODS <- c("mice_rubin", "mice_single", "cca")

# 6. 首次入组/入组次数限制敏感性
COPRIMARY_FIRST_PERIOD <- TRUE   # 仅使用 trial_num=1 作为 co-primary
MAX_ENROLLMENT_LEVELS  <- c(3L, 5L, 9999L)  # 限制最大入组次数的敏感性水平（9999 = 不限制）

# 7. 竞争风险敏感性
COMPETING_RISK_ENABLED <- TRUE

# 8. E-value 敏感性分析
EVALUE_ENABLED <- TRUE
EVALUE_DELTA <- c(1, 1.5, 2)  # HR 偏移量
EVALUE_STRONG_THRESHOLD    <- 2.0  # E-value 强稳健性阈值
EVALUE_MODERATE_THRESHOLD  <- 1.5  # E-value 中等稳健性阈值

# ---- 🆕 亚组分析完整配置（V26.4 修订） ----
# 预设亚组定义（基于先验临床假设，协议预注册级别）
# 每项亚组包含：变量名、分组条件、标签、预期方向
SUBGROUP_DEFINITIONS <- list(
  # === 血磷水平亚组（核心） ===
  list(var = "baseline_phosphate", condition = "< 0.50",  label = "中重度低磷 (<0.50)", expected = "benefit"),
  list(var = "baseline_phosphate", condition = ">= 0.50", label = "轻度低磷 (≥0.50)",  expected = "null"),
  
  # === 疾病严重度亚组 ===
  list(var = "sepsis3",     condition = "== 1",      label = "脓毒症",         expected = "harm"),
  list(var = "sepsis3",     condition = "== 0",      label = "非脓毒症",       expected = "null"),
  list(var = "sofa_t0",     condition = ">= 8",      label = "SOFA ≥ 8",      expected = "harm"),
  list(var = "sofa_t0",     condition = "< 8",       label = "SOFA < 8",      expected = "null"),
  
  # === 人口学亚组 ===
  list(var = "age",         condition = ">= 65",     label = "年龄 ≥ 65",      expected = "harm"),
  list(var = "age",         condition = "< 65",      label = "年龄 < 65",      expected = "null"),
  list(var = "gender_raw",  condition = "== 1",      label = "男性",           expected = "null"),
  list(var = "gender_raw",  condition = "== 0",      label = "女性",           expected = "null"),
  
  # === 治疗相关亚组 ===
  list(var = "invasive_vent_at_t0", condition = "== 1", label = "有创通气",    expected = "benefit"),
  list(var = "vasoactive_at_t0",    condition = "== 1", label = "血管活性药物", expected = "benefit"),
  
  # === 肾功能亚组 ===
  list(var = "egfr",        condition = "< 60",      label = "eGFR < 60",     expected = "harm"),
  list(var = "egfr",        condition = ">= 60",     label = "eGFR ≥ 60",     expected = "null"),
  
  # === 再喂养风险亚组 ===
  list(var = "refeeding_risk", condition = "== 1",   label = "再喂养高风险",   expected = "benefit"),
  list(var = "refeeding_risk", condition = "== 0",   label = "再喂养低风险",   expected = "null"),
  
  # === 营养支持亚组 ===
  list(var = "nutrition_support_type", condition = "!= 'none'", label = "营养支持", expected = "benefit"),
  list(var = "nutrition_support_type", condition = "== 'none'", label = "无营养支持", expected = "null")
)

# 亚组分析参数
SUBGROUP_MIN_SAMPLE  <- 50L    # 最小样本量
SUBGROUP_BOOT        <- 50L    # Bootstrap 次数
SUBGROUP_ALPHA       <- 0.05   # 显著性水平
SUBGROUP_MULTIPLE_CORRECTION <- "holm"  # 多重比较校正

# ---- 🆕 异质性检验配置 ----
HETEROGENEITY_MIN_STRATA <- 2L
HETEROGENEITY_METHOD     <- "DL"  # DerSimonian-Laird
HETEROGENEITY_I2_WARN    <- 50    # I² 警告阈值

# ---- 🆕 交互效应分析配置 ----
INTERACTION_COVARIATES <- c("baseline_phosphate", "sepsis3", "age", "egfr",
                            "invasive_vent_at_t0", "vasoactive_at_t0", "refeeding_risk")

# ---- 🆕 剂量-反应关系分析配置 ----
DOSE_RESPONSE_KNOTS <- 3L       # RCS 节点数
DOSE_RESPONSE_REF   <- 0.65    # 参考血磷值 (mmol/L)
DOSE_RESPONSE_BINS  <- c(0.32, 0.50, 0.65, 0.80)  # 血磷分箱

# ---- SEQTaRget 包配置 ----
# API 入口: SEQuential()（非 SEQdata()，后者为示例数据集）
# 正确用法: SEQuential(data, id.col, time.col, ..., method, options, verbose)
# 选项配置: SEQopts(bootstrap, weighted, weight.p99, ...)
# 结果提取: hazard_ratio(result), km_curve(result), diagnostics(result)
# 数据要求: period 必须从0开始（zero-indexed），ID/period 为整数型
# 变量数限制: 经测试验证，34个协变量可正常运行（tmp_seqtarget_varlimit_test.R v4）
# V26.4 方法学调整：文献标准以 ITT 为主分析，PP（censoring）为补充分析
SEQTAARGET_METHOD_MAIN <- "ITT"       # 主分析：意向性治疗 (ITT)
SEQTAARGET_METHOD_PP   <- "censoring" # 补充分析：符合方案效应 (PP / CCCW)
SEQTAARGET_METHOD      <- SEQTAARGET_METHOD_MAIN  # 默认向后兼容
SEQTAARGET_VERBOSE     <- FALSE
SEQTAARGET_NTHREADS    <- 2L          # V26.3: 4→2，控制内存避免OOM

# SEQTaRget 治疗变量水平映射（字符型治疗变量需映射为 0/1）
TREATMENT_LEVEL <- "treatment"     # 治疗组标签
CONTROL_LEVEL   <- "conservative"  # 对照组标签

# ---- 处理组/结局变量名（统一映射） ----
# 原始列名 → 分析用列名的映射
COL_TREATMENT_SRC <- "itt_group"
COL_OUTCOME_SRC   <- "death_28day_num"
COL_ID_SRC        <- "stay_id"
COL_PERIOD_SRC    <- "trial_num"

# ---- 发表级配色常量（规范八: JAMA Internal Medicine） ----
JOURNAL_COLORS <- list(
  treatment   = "#BC3C29",  # JAMA 红
  control     = "#0072B5",  # JAMA 蓝
  mortality   = "#D55E00",  # 橙色
  los         = "#0072B2",  # 深蓝
  respiratory = "#56B4E9",  # 浅蓝
  clinical    = "#009E73",  # 绿色
  renal       = "#CC79A7",  # 紫色
  e_weak      = "#F39C12",  # E-value 弱
  e_moderate  = "#E67E22",  # E-value 中等
  e_strong    = "#27AE60",  # E-value 强
  smd_line    = "#E74C3C",  # SMD 阈值线
  ref_line    = "#BDC3C7"   # 参考线
)

# ---- DAG14 协变量集（V26.9 新增 — 审稿 M-1.1 响应） ----
# DAG12 + 校正钙(corrected_calcium_mmol_min) + 营养支持类型(nutrition_support_type)
# 校正钙: 基线混杂，钙磷代谢耦合，0.7%缺失
# 营养支持类型: 医源性低磷重要驱动因素，0%缺失
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# ---- 🆕 Co-primary 分析配置（V26.4） ----
# 三重 Co-Primary: CCW + IPTW-Cox + Weibull AFT
COPRIMARY_RMST         <- TRUE    # RMST 作为 co-primary 估计量（PH 违反补充）
COPRIMARY_AFT          <- TRUE    # Weibull AFT 作为参数化验证

# ---- 🆕 TMLE 配置（V26.4） ----
TMLE_TAU     <- 28    # TMLE 分析时间窗口（天）

# ---- 🆕 竞争风险增强配置（V26.4） ----
COMPETING_RISK_CIF   <- TRUE  # 是否生成 CIF 曲线
COMPETING_RISK_RMTL  <- TRUE  # 是否计算 RMTL
COMPETING_RISK_CSH   <- TRUE  # 是否计算 Cause-specific Hazard

# ---- 🆕 血磷分层混杂诊断配置（V26.4） ----
PHOSPHATE_CONFOUNDING_PS_CUTOFF <- 0.20  # 分层 PS 诊断中 SMD 阈值
PHOSPHATE_CONFOUNDING_EVALUE_ENABLED <- TRUE  # 是否计算 E-value 评估未测量混杂

# ---- 硬盘缓存CCW配置（V26 新增） ----
# 缓存根目录
CCW_CACHE_DIR <- file.path(PROJECT_ROOT, "framework", "temp", "ccw_cache")

# 缓存有效期（小时）
CACHE_TTL_HOURS <- 24L

# 缓存强制刷新开关（通过环境变量 CCW_FORCE_REFRESH=TRUE 覆盖）
CCW_FORCE_REFRESH <- as.logical(Sys.getenv("CCW_FORCE_REFRESH", unset = "FALSE"))

# 缓存权重截断参数（默认值，可通过命令行参数覆盖）
CCW_TRUNC_LOWER <- 0.01
CCW_TRUNC_UPPER <- 0.99

# 缓存 Bootstrap 默认次数
CCW_BOOTSTRAP_N <- 100L

# 缓存随机种子
CCW_SEED <- 12345L

# 缓存并行线程数
CCW_NTHREADS <- 1L

# ---- 🆕 DAG 强制纳入变量（V26.4 修订，2026-06-21） ----
# 基于后门路径分析，以下变量是**治疗分配和死亡率的共同原因**，
# 必须无条件纳入 PS 模型，无论 LASSO 选择频率如何。
# 参考文献：Hernán & Robins (2020) §7.3; Greenland et al. (1999)
# 详细 DAG 文档见 TTE/v26/docs/v26_DAG.md
#
# 文献依据（新增变量）:
#   lactate_t0:
#     Suzuki et al. J Crit Care 2013 (PMID: 23265292): 低磷是疾病严重度标记物，
#       乳酸是休克灌注不足的直接标志，同时驱动细胞内磷转移(治疗指征)和死亡率
#     Padelli et al. Aust Crit Care 2021 (PMID: 32732023): 菌血症+低磷患者
#       校正严重度后 HR 2.10，乳酸是核心严重度混杂
#   liver_disease:
#     肝脏是磷代谢核心器官(糖原分解/合成)；肝病→磷代谢紊乱→补磷需求↑
#     肝病→凝血障碍、门脉高压→死亡率↑
#     Charlson 合并症指数中"中重度肝病"权重最高(3分)
#   aki_kdigo2_plus:
#     AKI→急性肾功能损伤→磷排泄/重吸收紊乱→补磷决策↑
#     AKI→死亡率独立相关(多篇ICU文献证实)
#     与 ever_crrt 互补：AKI=急性肾功能损伤，CRRT=治疗强度
DAG_FORCE_INCLUDE <- c(
  "age",                 # DAG 路径 1: 年龄→治疗, 年龄→死亡
  "baseline_sofa",       # DAG 路径 2: 疾病严重度→治疗, 疾病严重度→死亡
  "charlson",            # DAG 路径 3: 合并症负担→治疗, 合并症→死亡
  "sepsis3",             # DAG 路径 4: 脓毒症→治疗, 脓毒症→死亡
  "egfr",                # DAG 路径 5: 肾功能→治疗, 肾功能→死亡
  "ever_crrt",           # DAG 路径 6: CRRT→低磷(清除过滤)→治疗; CRRT→死亡(严重度标志)
  "invasive_vent_at_t0", # DAG 路径 7: 机械通气→治疗, 通气→死亡
  "vasoactive_at_t0",    # DAG 路径 8: 血管活性药→治疗, 药物→死亡
  "refeeding_risk",      # DAG 路径 9: 再喂养风险→治疗(指征), 风险→死亡
  "lactate_t0",          # DAG 路径10: 乳酸(灌注不足)→细胞内磷转移→治疗; 乳酸→死亡率(强预测)
  "liver_disease",       # DAG 路径11: 肝病→磷代谢紊乱→治疗; 肝病→凝血障碍/死亡率↑
  "aki_kdigo2_plus"      # DAG 路径12: AKI→肾性磷紊乱→治疗; AKI→死亡率↑
)

# ---- 🆕 分批 CCW 引擎配置（V26.4 - 替代 SEQTaRget） ----
# 分批引擎替代 SEQTaRget 作为主分析：
#   扩展 → 逐期 PS/IPTW → 逐期 Cox → Meta 汇总
#   优势: 内存可控（单期 O(n)），无 OOM 风险，支持任意协变量数
#   限制: 无 Bootstrap CI（使用聚类稳健 SE + Meta 合并）
CCW_BATCH_ENABLED     <- TRUE         # 启用分批引擎
CCW_BATCH_POOL_METHOD <- "fixed"      # Meta 汇总方法: "fixed" | "random"
CCW_BATCH_STABILIZE   <- TRUE         # 启用稳定化 IPTW
CCW_BATCH_COV_SET     <- "C1"          # 协变量集编号（同 SEQTaRget）
CCW_BATCH_OUTPUT_DIR  <- file.path(DATA_DIR, "batch_ccw")

# ---- 缓存协变量集默认配置
CCW_COVARIATE_SET <- "C1"

# ---- 🆕 renv 环境管理配置（P2-7修复） ----
# renv 初始化状态文件路径
RENV_LOCKFILE <- file.path(PROJECT_ROOT, "renv.lock")
RENV_ACTIVATE  <- file.path(PROJECT_ROOT, "renv", "activate.R")

# 检查 renv 是否已初始化
renv_check <- function() {
  if (file.exists(RENV_LOCKFILE)) {
    cat(sprintf("[config.R] renv 环境已初始化: %s\n", RENV_LOCKFILE))
    return(TRUE)
  } else {
    cat("[config.R] 警告: renv 未初始化，建议运行 renv::init() 锁定依赖版本\n")
    cat("[config.R] 最小依赖: data.table, survival, arrow, logistf, mice, glmnet, sandwich, rms, metafor\n")
    return(FALSE)
  }
}

# 自动检查 renv 状态（启动时静默检查）
RENV_AVAILABLE <- renv_check()

# ---- 日志启动信息 ----
log_startup <- function(script_name) {
  log_section(paste("V26", script_name, "- 分析启动"))
  log_info(sprintf("项目根目录: %s", PROJECT_ROOT))
  log_info(sprintf("输出目录: %s", DATA_DIR))
  log_info(sprintf("启动时间: %s", Sys.time()))
}

cat(sprintf("[config.R] 项目根目录: %s\n", PROJECT_ROOT))
cat(sprintf("[config.R] 协变量候选池: %d 个变量\n", length(COVARIATE_POOL)))