# ============================================================================
# V26.9 多重比较校正脚本 — FDR (交互检验) + Bonferroni (次要结局)
#
# 审稿要求 M4: 多重比较校正
#   Part 1: 10变量交互检验 → Benjamini-Hochberg FDR
#     - 对每个交互变量运行: coxph(Surv(...) ~ tx_int * iv + DAG14_excl, cluster = stay_id)
#     - 提取交互项 p 值, p.adjust(..., method = "BH")
#   Part 2: 5项次要结局 → Bonferroni 校正 (α=0.01)
#     - 对每项次要结局运行: glm(y ~ tx_int + DAG14, family = binomial)
#     - Bonferroni 阈值 = 0.05/5 = 0.01
#
# 数据: sequential_trials_imputed_v26.parquet (MIMIC-IV, MICE插补后)
# 输出: TTE/v26/output/v26_9/
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
})

# ============================================================================
# 1. 参数配置
# ============================================================================

# ---- 输出目录 ----
OUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
TBL_DIR <- file.path(OUT_DIR, "tables")
LOG_DIR <- file.path(OUT_DIR, "logs")
for (d in c(TBL_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 数据路径 ----
DATA_FILE <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"

# ---- 10个交互变量（同 V26.8 Phase 3.2） ----
INTERACTION_VARS <- c(
  "baseline_sofa", "age", "egfr", "charlson",
  "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "invasive_vent_at_t0", "vasoactive_at_t0", "sepsis3"
)

# ---- DAG14 协变量集（来自 config.R 第376-381行） ----
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# ---- 5项次要结局变量映射（审稿 M4 指定） ----
SECONDARY_OUTCOMES <- list(
  list(name = "AKI (KDIGO ≥2期)",          var = "aki_kdigo2_plus"),
  list(name = "CV复合结局",                 var = "cardiovascular_event_comprehensive"),
  list(name = "新发心律失常",               var = "new_arrhythmia_rhythm"),
  list(name = "心肌损伤",                   var = "myocardial_injury"),
  list(name = "输血",                       var = "transfusion_at_t0")
)
N_SECONDARY <- length(SECONDARY_OUTCOMES)               # 5
BONFERRONI_ALPHA <- 0.05 / N_SECONDARY                  # 0.01

# ---- 分析参数 ----
SEED <- 20260629L
set.seed(SEED)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_multiple_correction_%s.log",
                                        format(Sys.time(), "%Y%m%d_%H%M%S")))
log_con  <- file(log_file, open = "wt")

log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
      file = log_con, append = TRUE)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

log_section <- function(title) {
  cat(sprintf("\n========== %s ==========\n", title),
      file = log_con, append = TRUE)
  cat(sprintf("\n========== %s ==========\n", title))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

log_info("V26.9 多重比较校正脚本启动")
log_info(sprintf("种子: %d", SEED))
log_info(sprintf("输出目录: %s", OUT_DIR))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")

dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int = as.integer(itt_group == "treatment"),
  event  = death_28day_num
)]
setorder(dt_full, stay_id, trial_num)

log_info(sprintf("Parquet 数据: %d 行, %d 事件, %d 列",
                 nrow(dt_full), sum(dt_full$event), ncol(dt_full)))

# ---- 检查 DAG14 变量可用性 ----
dag14_avail <- intersect(DAG14, names(dt_full))
dag14_missing <- setdiff(DAG14, names(dt_full))
log_info(sprintf("DAG14 可用协变量: %d/%d", length(dag14_avail), length(DAG14)))
if (length(dag14_missing) > 0) {
  log_info(sprintf("DAG14 缺失变量（从 DAG14 排除）: %s",
                   paste(dag14_missing, collapse = ", ")))
}

# ---- 检查交互变量可用性 ----
iv_avail <- intersect(INTERACTION_VARS, names(dt_full))
iv_missing <- setdiff(INTERACTION_VARS, names(dt_full))
log_info(sprintf("交互变量可用: %d/%d", length(iv_avail), length(INTERACTION_VARS)))
if (length(iv_missing) > 0) {
  log_info(sprintf("交互变量缺失: %s", paste(iv_missing, collapse = ", ")))
}

# ---- 检查次要结局变量可用性 ----
so_vars <- sapply(SECONDARY_OUTCOMES, `[[`, "var")
so_avail <- intersect(so_vars, names(dt_full))
so_missing <- setdiff(so_vars, names(dt_full))
log_info(sprintf("次要结局变量可用: %d/%d", length(so_avail), length(so_vars)))
if (length(so_missing) > 0) {
  log_info(sprintf("次要结局缺失: %s", paste(so_missing, collapse = ", ")))
}

# ---- 字符型协变量转 factor ----
for (v in dag14_avail) {
  if (is.character(dt_full[[v]])) {
    dt_full[[v]] <- as.factor(dt_full[[v]])
  }
}

# ============================================================================
# 4. Part 1: 交互检验 — Benjamini-Hochberg FDR 校正
# ============================================================================
log_section("Part 1: 10变量交互检验 — Benjamini-Hochberg FDR 校正")

interaction_vars_active <- iv_avail
log_info(sprintf("将检验 %d 个交互变量: %s",
                 length(interaction_vars_active),
                 paste(interaction_vars_active, collapse = ", ")))

interaction_raw <- list()

for (iv in interaction_vars_active) {
  # ---- 4a. 变量预处理 ----
  # 字符型转 factor
  if (is.character(dt_full[[iv]])) {
    dt_full[[iv]] <- as.factor(dt_full[[iv]])
  }

  # 连续变量中心化（避免与主效应共线性）
  if (is.numeric(dt_full[[iv]])) {
    cen_var <- paste0(iv, "_c")
    dt_full[[cen_var]] <- scale(dt_full[[iv]], scale = FALSE)[, 1]
    iv_col <- cen_var
  } else {
    iv_col <- iv
  }

  # ---- 4b. DAG14 扣除当前交互变量 ----
  # 排除当前交互变量本身，避免完全共线性
  dag14_excl <- setdiff(dag14_avail, iv)
  # 如果用了中心化变量，也排除原始变量
  if (iv_col != iv) {
    dag14_excl <- setdiff(dag14_excl, iv)
  }

  # 再扣除字符/因子型变量中仅有单水平的
  dag14_excl <- dag14_excl[sapply(dag14_excl, function(v) {
    x <- dt_full[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x))  return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  # ---- 4c. Cox 交互检验 ----
  dag14_str <- paste(dag14_excl, collapse = " + ")
  formula_str <- sprintf("Surv(surv_time, event) ~ tx_int * %s + %s", iv_col, dag14_str)

  fit <- tryCatch(
    coxph(as.formula(formula_str), data = dt_full, cluster = stay_id),
    error = function(e) {
      log_info(sprintf("[交互] %s: Cox 失败 — %s", iv, e$message))
      return(NULL)
    }
  )

  if (!is.null(fit)) {
    int_names <- grep(":", names(coef(fit)), value = TRUE)
    if (length(int_names) > 0) {
      s <- summary(fit)
      int_hr  <- exp(coef(fit)[int_names[1]])
      int_ci  <- exp(confint(fit)[int_names[1], ])
      int_p   <- s$coefficients[int_names[1], "Pr(>|z|)"]

      log_info(sprintf("[交互] tx × %s: HR=%.4f (%.4f–%.4f), p=%.6f",
                       iv, int_hr, int_ci[1], int_ci[2], int_p))

      interaction_raw[[iv]] <- data.table(
        var       = iv,
        hr        = int_hr,
        ci_lower  = int_ci[1],
        ci_upper  = int_ci[2],
        raw_p     = int_p
      )
    } else {
      log_info(sprintf("[交互] %s: 模型成功但未检测到交互项", iv))
    }
  }

  # ---- 4d. 清理临时变量 ----
  if (iv_col != iv && iv_col %in% names(dt_full)) {
    dt_full[, (iv_col) := NULL]
  }
}

n_int <- length(interaction_raw)
log_info(sprintf("交互检验完成: %d/%d 成功获取 p 值",
                 n_int, length(interaction_vars_active)))

# ---- 4e. Benjamini-Hochberg FDR 校正 ----
if (n_int >= 2) {
  dt_int_fdr <- rbindlist(interaction_raw)

  # BH 校正
  dt_int_fdr[, `:=`(
    q_value   = p.adjust(raw_p, method = "BH"),
    rank_by_p = frank(raw_p)
  )]
  setorder(dt_int_fdr, raw_p)

  # 显著性标注
  dt_int_fdr[, significance_label := fcase(
    q_value < 0.10, "FDR-significant",
    raw_p   < 0.05, "Nominal",
    default =         "NS"
  )]

  # ---- 4f. 控制台打印汇总 ----
  log_info(sprintf("\n=== FDR 校正结果 (BH, %d 检验) ===", n_int))
  cat(sprintf("\n%-25s %8s %10s %10s %10s %10s %s\n",
              "变量", "HR", "CI_lower", "CI_upper", "raw_p", "q_value", "标注"),
      file = log_con, append = TRUE)
  cat(sprintf("\n%-25s %8s %10s %10s %10s %10s %s\n",
              "变量", "HR", "CI_lower", "CI_upper", "raw_p", "q_value", "标注"))
  for (i in seq_len(nrow(dt_int_fdr))) {
    r <- dt_int_fdr[i]
    line <- sprintf("  %-25s %8.4f %10.4f %10.4f %10.6f %10.6f %s",
                    r$var, r$hr, r$ci_lower, r$ci_upper,
                    r$raw_p, r$q_value, r$significance_label)
    cat(line, "\n", file = log_con, append = TRUE)
    cat(line, "\n")
  }

  n_fdr <- sum(dt_int_fdr$q_value < 0.10)
  n_nom <- sum(dt_int_fdr$raw_p < 0.05)
  log_info(sprintf("\nFDR 汇总: %d/%d FDR-significant (q<0.10), %d/%d Nominal (raw_p<0.05)",
                   n_fdr, n_int, n_nom, n_int))

  # ---- 4g. 保存 CSV ----
  fwrite(dt_int_fdr, file.path(TBL_DIR, "v26_9_interaction_fdr.csv"))
  log_info(sprintf("FDR 结果已保存: %s",
                   file.path(TBL_DIR, "v26_9_interaction_fdr.csv")))

} else {
  log_info("FDR 校正跳过: 不足 2 个成功交互检验（无法进行有意义的校正）")
  dt_int_fdr <- data.table()
}

# ============================================================================
# 5. Part 2: 次要结局 — Bonferroni 校正
# ============================================================================
log_section(sprintf("Part 2: %d项次要结局 — Bonferroni 校正 (α=%.4f)",
                    N_SECONDARY, BONFERRONI_ALPHA))

secondary_raw <- list()

for (so in SECONDARY_OUTCOMES) {
  outcome_var  <- so$var
  outcome_name <- so$name

  # ---- 5a. 变量存在性检查 ----
  if (!outcome_var %in% names(dt_full)) {
    log_info(sprintf("[次要] %s (%s): 数据中不存在，跳过", outcome_name, outcome_var))
    next
  }

  # ---- 5b. 事件数检查 ----
  n_event <- sum(dt_full[[outcome_var]] == 1, na.rm = TRUE)
  n_total <- sum(!is.na(dt_full[[outcome_var]]))
  if (n_event < 10) {
    log_info(sprintf("[次要] %s: 事件数不足 (events=%d), 跳过", outcome_name, n_event))
    next
  }

  # ---- 5c. DAG14 协变量准备（扣除当前结局变量） ----
  dag14_for_so <- setdiff(dag14_avail, outcome_var)

  # 排除单水平因子
  dag14_for_so <- dag14_for_so[sapply(dag14_for_so, function(v) {
    x <- dt_full[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x))  return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  dag14_so_str <- paste(dag14_for_so, collapse = " + ")
  formula_str <- sprintf("%s ~ tx_int + %s", outcome_var, dag14_so_str)

  # ---- 5d. Logistic 回归 ----
  fit <- tryCatch(
    glm(as.formula(formula_str), data = dt_full, family = binomial),
    error = function(e) {
      log_info(sprintf("[次要] %s: GLM 失败 — %s", outcome_name, e$message))
      return(NULL)
    }
  )

  if (!is.null(fit)) {
    s <- summary(fit)

    # 确认 tx_int 系数存在
    if (!"tx_int" %in% rownames(s$coefficients)) {
      log_info(sprintf("[次要] %s: 模型中无 tx_int 系数", outcome_name))
      next
    }

    coef_row <- s$coefficients["tx_int", ]
    or_val   <- exp(coef_row["Estimate"])
    or_ci    <- tryCatch(
      exp(confint(fit, "tx_int")),
      error = function(e) c(NA_real_, NA_real_)
    )
    raw_p    <- coef_row["Pr(>|z|)"]

    log_info(sprintf("[次要] %s: OR=%.4f (%.4f–%.4f), p=%.6f, events=%d/%d",
                     outcome_name, or_val, or_ci[1], or_ci[2],
                     raw_p, n_event, n_total))

    secondary_raw[[outcome_var]] <- data.table(
      outcome      = outcome_name,
      outcome_var  = outcome_var,
      or           = or_val,
      ci_lower     = or_ci[1],
      ci_upper     = or_ci[2],
      raw_p        = raw_p,
      n_events     = n_event,
      n_total      = n_total
    )
  }
}

n_so <- length(secondary_raw)
log_info(sprintf("次要结局分析完成: %d/%d 成功", n_so, N_SECONDARY))

# ---- 5e. Bonferroni 标注 ----
if (n_so >= 1) {
  dt_sec_bonf <- rbindlist(secondary_raw)
  dt_sec_bonf[, `:=`(
    bonferroni_threshold = BONFERRONI_ALPHA,
    significant_yn       = fifelse(raw_p < BONFERRONI_ALPHA, "Yes", "No")
  )]
  setorder(dt_sec_bonf, raw_p)

  # 控制台打印
  log_info(sprintf("\n=== Bonferroni 校正结果 (阈值=%.4f) ===", BONFERRONI_ALPHA))
  for (i in seq_len(nrow(dt_sec_bonf))) {
    r <- dt_sec_bonf[i]
    sig_mark <- if (r$significant_yn == "Yes") "★ Bonferroni 显著" else "NS"
    line <- sprintf("  %-25s OR=%.4f (%.4f–%.4f)  raw_p=%.6f  %s",
                    r$outcome, r$or, r$ci_lower, r$ci_upper,
                    r$raw_p, sig_mark)
    cat(line, "\n", file = log_con, append = TRUE)
    cat(line, "\n")
  }

  n_bonf <- sum(dt_sec_bonf$significant_yn == "Yes")
  log_info(sprintf("\nBonferroni 汇总: %d/%d 显著 (raw_p < %.4f)",
                   n_bonf, n_so, BONFERRONI_ALPHA))

  # 保存 CSV
  fwrite(dt_sec_bonf, file.path(TBL_DIR, "v26_9_secondary_bonferroni.csv"))
  log_info(sprintf("Bonferroni 结果已保存: %s",
                   file.path(TBL_DIR, "v26_9_secondary_bonferroni.csv")))

} else {
  log_info("Bonferroni 校正跳过: 无成功次要结局分析")
  dt_sec_bonf <- data.table()
}

# ============================================================================
# 6. 保存完整 RDS 结果对象
# ============================================================================
log_section("结果保存")

rds_out <- list(
  version              = "V26.9",
  script               = "v26_9_multiple_correction.R",
  timestamp            = as.character(Sys.time()),
  seed                 = SEED,
  dag14                = DAG14,
  dag14_available      = dag14_avail,
  dag14_missing        = dag14_missing,

  # Part 1: 交互 FDR
  interaction_vars          = INTERACTION_VARS,
  interaction_vars_active   = interaction_vars_active,
  interaction_n_success     = n_int,
  interaction_fdr           = dt_int_fdr,

  # Part 2: 次要结局 Bonferroni
  secondary_outcomes_config = SECONDARY_OUTCOMES,
  secondary_n_success       = n_so,
  secondary_bonferroni      = dt_sec_bonf,
  bonferroni_alpha          = BONFERRONI_ALPHA
)

saveRDS(rds_out, file.path(TBL_DIR, "v26_9_multiple_correction.rds"))
log_info(sprintf("完整 RDS 已保存: %s",
                 file.path(TBL_DIR, "v26_9_multiple_correction.rds")))

# ============================================================================
# 7. 流程自审核
# ============================================================================
log_section("流程自审核")

audit <- function(check, pass, detail = "") {
  mark <- if (pass) "PASS" else "FAIL"
  log_info(sprintf("[%s] %s: %s", mark, check, detail))
  return(pass)
}

all_pass <- TRUE
all_pass <- all_pass && audit("数据加载",
  nrow(dt_full) > 10000, sprintf("n=%d", nrow(dt_full)))
all_pass <- all_pass && audit("DAG14 可用性",
  length(dag14_avail) >= 12, sprintf("%d/14", length(dag14_avail)))
all_pass <- all_pass && audit("交互 FDR 校正",
  n_int >= 2, sprintf("%d/10 成功", n_int))
all_pass <- all_pass && audit("Bonferroni 校正",
  n_so >= 1, sprintf("%d/5 成功", n_so))
all_pass <- all_pass && audit("FDR 结果 CSV",
  file.exists(file.path(TBL_DIR, "v26_9_interaction_fdr.csv")), "")
all_pass <- all_pass && audit("Bonferroni 结果 CSV",
  file.exists(file.path(TBL_DIR, "v26_9_secondary_bonferroni.csv")), "")
all_pass <- all_pass && audit("完整 RDS",
  file.exists(file.path(TBL_DIR, "v26_9_multiple_correction.rds")), "")
all_pass <- all_pass && audit("日志文件",
  file.exists(log_file), "")

log_info(sprintf("\n审核结果: %s",
                 if (all_pass) "全部通过 ✓" else "存在未通过项 ✗"))

# ============================================================================
# 8. 收尾
# ============================================================================
close(log_con)

cat("\n")
cat("################################################################################\n")
cat("      V26.9 多重比较校正完成\n")
cat("      FDR (交互检验, BH) + Bonferroni (次要结局, α=0.01)\n")
cat("################################################################################\n\n")

cat(sprintf("Part 1 — 交互检验 FDR:\n"))
cat(sprintf("  成功/总: %d/%d\n", n_int, length(INTERACTION_VARS)))
if (n_int >= 2) {
  cat(sprintf("  FDR-significant (q<0.10): %d\n",
              sum(dt_int_fdr$q_value < 0.10)))
  cat(sprintf("  Nominal (raw_p<0.05): %d\n",
              sum(dt_int_fdr$raw_p < 0.05)))
}

cat(sprintf("\nPart 2 — 次要结局 Bonferroni (α=%.4f):\n", BONFERRONI_ALPHA))
cat(sprintf("  成功/总: %d/%d\n", n_so, N_SECONDARY))
if (n_so >= 1) {
  cat(sprintf("  显著: %d\n", sum(dt_sec_bonf$significant_yn == "Yes")))
}

cat(sprintf("\n日志: %s\n", log_file))
cat(sprintf("交互 FDR:  %s\n",
            file.path(TBL_DIR, "v26_9_interaction_fdr.csv")))
cat(sprintf("Bonferroni: %s\n",
            file.path(TBL_DIR, "v26_9_secondary_bonferroni.csv")))
cat(sprintf("RDS:        %s\n",
            file.path(TBL_DIR, "v26_9_multiple_correction.rds")))
cat("========== V26.9 多重比较校正完成 ==========\n")
