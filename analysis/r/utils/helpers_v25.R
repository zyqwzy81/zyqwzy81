# ============================================================================
# V25 统一辅助函数集
# 所有 V25 分析脚本通过 source("_helpers_v25.R") 获取以下函数
#
# 引用链：_helpers_v25.R → _helpers_v24.R → _helpers_v23.R → _helpers_v22.R → _helpers.R
# ============================================================================

# ---- 历史版本辅助函数（按需加载，通过 source() 引用） ----
# 注意：V25 框架中所有核心函数已在本文件中定义，无需加载历史版本
# 历史版本函数（如 V21 PS Variable Map）在各版本脚本中按需 source
# ... (V21/V22/V23 函数定义已整合到本文件)

# ============================================================================
# 日志工具函数（统一版本，所有脚本共用）
# ============================================================================

# 输出章节标题
log_section <- function(title) {
  cat(sprintf("\n============================================================\n=== %s ===\n============================================================\n", title))
}

# 输出信息日志
log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# 输出警告日志
log_warn <- function(msg) {
  cat(sprintf("[%s] WARNING: %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# 输出错误日志
log_error <- function(msg) {
  cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# ============================================================================
# SMD 计算函数（统一版本）
# ============================================================================

# 计算单个变量的标准化均值差异 (Standardized Mean Difference)
# 参数:
#   x     : 协变量向量
#   treat : 处理组指示向量 (0/1)
#   w     : 权重向量 (IPTW)，为 NULL 表示未加权
#   var_name: 变量名（可选，用于日志记录）— P2-8
# 返回:
#   SMD 值（数值）
calc_smd <- function(x, treat, w = NULL, var_name = NULL) {
  # 🆕 P2-8: 变量名日志记录（仅在首次调用时记录）
  if (!is.null(var_name)) {
    log_info(sprintf("SMD 计算: %s", var_name))
  }
  x0 <- x[treat == 0]
  x1 <- x[treat == 1]
  if (!is.null(w)) {
    w0 <- w[treat == 0]
    w1 <- w[treat == 1]
  }
  
  # 判断是否为二值/分类变量
  x_vals <- unique(na.omit(x))
  is_binary <- all(x_vals %in% c(0, 1, NA)) || length(x_vals) <= 5
  
  if (is_binary) {
    # 分类/二值变量：使用比例差
    if (is.null(w)) {
      p0 <- mean(x0, na.rm = TRUE)
      p1 <- mean(x1, na.rm = TRUE)
    } else {
      p0 <- weighted.mean(x0, w0, na.rm = TRUE)
      p1 <- weighted.mean(x1, w1, na.rm = TRUE)
    }
    if (is.na(p0) || is.na(p1) || is.nan(p0) || is.nan(p1)) return(0)
    pooled_sd <- sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / 2)  # P1-1修复: pooled SD除以2
    if (is.na(pooled_sd) || pooled_sd < 1e-10) return(0)
    return((p1 - p0) / pooled_sd)
  } else {
    # 连续变量：使用加权 Cohen's d
    if (is.null(w)) {
      m0 <- mean(x0, na.rm = TRUE)
      m1 <- mean(x1, na.rm = TRUE)
      s0 <- sd(x0, na.rm = TRUE)
      s1 <- sd(x1, na.rm = TRUE)
    } else {
      m0 <- weighted.mean(x0, w0, na.rm = TRUE)
      m1 <- weighted.mean(x1, w1, na.rm = TRUE)
      s0 <- sqrt(sum(w0 * (x0 - m0)^2, na.rm = TRUE) / sum(w0))
      s1 <- sqrt(sum(w1 * (x1 - m1)^2, na.rm = TRUE) / sum(w1))
    }
    pooled_sd <- sqrt((s0^2 + s1^2) / 2)
    if (is.na(pooled_sd) || pooled_sd < 1e-10) return(0)
    return((m1 - m0) / pooled_sd)
  }
}

# ============================================================================
# 变量名美化函数
# ============================================================================

# 将内部变量名映射为中文/可读标签
beautify_var <- function(v) {
  labels <- c(
    age = "年龄", gender = "性别", gender_female = "女性",
    baseline_phosphate = "基线血磷", baseline_sofa = "基线 SOFA",
    baseline_sofa_t0 = "基线 SOFA (T0)", admission_sofa = "入院 SOFA",
    apsiii = "APS III", sepsis3 = "Sepsis-3",
    egfr = "eGFR", invasive_vent_at_t0 = "有创通气",
    vasoactive_at_t0 = "血管活性药物", albumin_min = "白蛋白(min)",
    charlson = "Charlson合并症指数", charlson_comorbidity_index = "Charlson合并症指数",
    liver_disease = "肝病", insulin_at_t0 = "胰岛素",
    diuretic_at_t0 = "利尿剂", steroid_at_t0 = "类固醇",
    nahco3_infusion_at_t0 = "碳酸氢钠输注",
    beta2_agonist_at_t0 = "β2激动剂", ppi_use = "PPI使用",
    ph_worst = "最差pH", bicarbonate_worst = "最差HCO3",
    alkalosis = "碱中毒", resp_alkalosis = "呼吸性碱中毒",
    metabolic_alkalosis = "代谢性碱中毒",
    net_fluid_balance_24h = "24h净液体平衡", transfusion_at_t0 = "输血",
    cardiac_surgery = "心脏手术", aki_kdigo2_plus = "AKI KDIGO 2+",
    ever_crrt = "CRRT史", corrected_calcium_mmol_min = "校正钙(mmol/L)",
    admission_type = "入院类型",
    refeeding_risk = "再喂养风险",
    nutrition_support_type = "营养支持类型",
    nutrition_any = "任意营养支持",
    en_kcal_24h = "EN热卡(T0前24h)",
    tpn_kcal_24h = "TPN热卡(T0前24h)",
    sf_ratio = "S/F比值", sf_estimated_pf = "估算P/F",
    sf_fio2_defaulted = "FiO2默认值标记",
    # 🆕 V25.4: ICU类型 + 补充变量
    icu_micu = "MICU", icu_sicu = "SICU", icu_cvicu = "CVICU",
    icu_ccu = "CCU", icu_neuro = "Neuro ICU",
    magnesium_t0 = "血镁(T0)", phosphate_nadir = "血磷最低值",
    alcohol_abuse = "酒精滥用",
    peep_t0 = "PEEP(T0)", fio2_t0 = "FiO2(T0)",
    ventilation_duration_hours = "通气时长(h)",
    # 🆕 P2-5修复: 补齐 COVARIATE_POOL (30个) 中缺失的变量标签
    albumin_min = "白蛋白(min)", hemoglobin_min = "血红蛋白(min)",
    lactate_t0 = "乳酸(T0)", apsiii = "APS III",
    # 🆕 增强数据集扩展变量
    gnri = "GNRI", pni = "PNI", conut_score = "CONUT评分",
    pao2fio2ratio = "P/F比值", nutrition_any = "任意营养支持",
    has_oral_ng_phosphate_48h = "口服/鼻饲磷(48h)",
    cumulative_po_phosphate_mmol_48h = "口服磷累积剂量(48h)",
    iv_within_12h = "12h内IV补磷", iv_within_24h = "24h内IV补磷"
  )
  if (v %in% names(labels)) labels[v] else v
}

# ============================================================================
# IPTW 权重计算函数
# ============================================================================

# 计算稳定化 IPTW 权重
# 参数:
#   data       : data.table，需含 treatment 列和所有协变量
#   covariates : 协变量名称向量
#   trunc_p    : 权重截断分位数 (NULL = 不截断)
#   stabilize  : 是否使用稳定化权重
compute_iptw <- function(data, covariates, trunc_p = 0.99, lower_p = 0.01, stabilize = TRUE, bilateral = TRUE, treat_col = "tx_init") {
  covs <- intersect(covariates, names(data))
  if (length(covs) < 2) {
    log_warn("协变量不足 (需 >=2)，返回均匀权重")
    return(rep(1, nrow(data)))
  }
  
  # 🆕 P2-8: 记录候选变量名列表
  log_info(sprintf("IPTW 候选协变量 (%d): %s", length(covs), paste(covs, collapse = ", ")))
  
  # 确保处理列存在
  if (!(treat_col %in% names(data))) {
    if ("treatment" %in% names(data)) treat_col <- "treatment"
    else if ("tx_init" %in% names(data)) treat_col <- "tx_init"
    else stop("无法找到处理变量列")
  }
  
  treat_vec <- data[[treat_col]]
  
  ps_formula <- as.formula(paste(treat_col, "~", paste(covs, collapse = " + ")))
  ps_glm <- glm(ps_formula, data = data, family = binomial(link = "logit"))
  ps <- predict(ps_glm, type = "response")
  
  if (stabilize) {
    p_treat <- mean(treat_vec == 1, na.rm = TRUE)
    p_ctrl  <- 1 - p_treat
    w <- fifelse(treat_vec == 1, p_treat / ps, p_ctrl / (1 - ps))
  } else {
    w <- fifelse(treat_vec == 1, 1 / ps, 1 / (1 - ps))
  }
  
  # 🆕 双侧截断（Cole & Hernán, 2008）
  if (!is.null(trunc_p) && !is.na(trunc_p)) {
    w_upper <- quantile(w, trunc_p, na.rm = TRUE)
    w <- pmin(w, w_upper)
  }
  if (bilateral && !is.null(lower_p) && !is.na(lower_p)) {
    w_lower <- quantile(w, lower_p, na.rm = TRUE)
    w <- pmax(w, w_lower)
  }
  
  return(w)
}

# ============================================================================
# E-value 计算函数（使用标准公式）
# ============================================================================

# 参考: VanderWeele & Ding (2017), Annals of Internal Medicine
compute_evalue <- function(hr, ci_lower, ci_upper) {
  # 对于 HR < 1，取倒数
  hr_for_eval <- if (hr < 1) 1 / hr else hr
  
  # E-value (点估计): HR + sqrt(HR*(HR-1))
  e_point <- if (hr_for_eval > 1) hr_for_eval + sqrt(hr_for_eval * (hr_for_eval - 1)) else 1
  
  # E-value (CI 下界)
  if (hr > 1) {
    ci_bound <- ci_lower
  } else {
    ci_bound <- 1 / ci_upper
  }
  
  if (is.na(ci_bound) || ci_bound <= 1) {
    e_ci <- 1
  } else {
    e_ci <- ci_bound + sqrt(ci_bound * (ci_bound - 1))
  }
  
  # 分级
  grade <- if (e_point >= 2.0) {
    "强 (E>=2.0)"
  } else if (e_point >= 1.5) {
    "中等 (1.5<=E<2.0)"
  } else {
    "弱 (E<1.5)"
  }
  
  list(
    e_value_point = round(e_point, 2),
    e_value_ci    = round(e_ci, 2),
    grade         = grade
  )
}

# ============================================================================
# 数据加载函数（支持 Parquet 和 CSV）
# ============================================================================

load_data <- function(data_path) {
  if (grepl("\\.parquet$", data_path, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("需要 arrow 包来读取 Parquet 文件。请运行: install.packages('arrow')")
    }
    library(arrow)
    dt <- as.data.table(read_parquet(data_path))
  } else {
    dt <- fread(data_path)
  }
  log_info(sprintf("数据加载: %s 行 x %d 列", format(nrow(dt), big.mark = ","), ncol(dt)))
  return(dt)
}

# ============================================================================
# 变量分类函数
# ============================================================================

# 根据变量名后缀自动分类为时变变量(time_varying)和基线变量(fixed)
classify_variables <- function(vars) {
  tv_patterns <- c("_at_t0$", "_24hr$", "_48h$", "_any_48h$", "_dose$",
                   "_daily$", "_at_entry$", "_bl$", "_tight$")
  fixed_patterns <- c("^age$", "^gender", "^charlson", "^admission",
                      "^ethnicity", "^insurance", "^height",
                      "^weight", "^bmi$", "^ckd", "^comorbidity",
                      "^race", "^icu_type", "^icu_", "^apsiii$", "^oasis$",
                      "^sf_ratio", "^sf_", "^alcohol_", "^phosphate_nadir",
                      "^ventilation_duration")
  
  time_varying <- character(0)
  fixed <- character(0)
  
  for (v in vars) {
    is_tv <- FALSE
    for (pat in tv_patterns) {
      if (grepl(pat, v, ignore.case = TRUE)) { is_tv <- TRUE; break }
    }
    if (is_tv) { time_varying <- c(time_varying, v); next }
    for (pat in fixed_patterns) {
      if (grepl(pat, v, ignore.case = TRUE)) { fixed <- c(fixed, v); break }
    }
  }
  
  unclassified <- setdiff(vars, c(time_varying, fixed))
  if (length(unclassified) > 0) {
    time_varying <- c(time_varying, unclassified)
  }
  
  list(time_varying = time_varying, fixed = fixed)
}

# ===========================================================================
# 🆕 血磷分层分类函数（2026-06 重构版）
# ===========================================================================

# 将基线血磷值映射为3层分类（主分析）
# <0.50 (中重度), 0.50-0.65 (轻中度), 0.65-0.80 (轻度)
classify_phosphate_level_3 <- function(phosphate_val) {
  ifelse(phosphate_val < 0.50, "mod_severe",
    ifelse(phosphate_val < 0.65, "mild_mod", "mild"))
}

# 将基线血磷值映射为4层分类（保留，用于敏感性分析）
# <0.32 (严重), 0.32-0.50 (中度), 0.50-0.65 (轻中度), 0.65-0.80 (轻度)
classify_phosphate_level <- function(phosphate_val) {
  ifelse(phosphate_val < 0.32, "severe",
    ifelse(phosphate_val < 0.50, "moderate",
      ifelse(phosphate_val < 0.65, "mild_mod", "mild")))
}

# ===========================================================================
# 辅助函数加载完成（基础模块）
# ===========================================================================
log_info("_helpers_v25.R 基础模块加载完成（含日志、SMD、IPTW、E-value、数据加载、变量分类、血磷分层）")

# ===========================================================================
# 🆕 2026-06 重构：亚组分析、异质性检验、交互效应、剂量-反应关系模块
# ===========================================================================

# ---- 亚组筛选函数（安全的数据子集提取） ----
# 参数:
#   data         : data.table
#   var_name     : 筛选变量名
#   condition    : 条件表达式字符串（如 ">= 65"）
#   label        : 亚组标签
# 返回:
#   筛选后的 data.table
get_subgroup <- function(data, var_name, condition, label) {
  if (!(var_name %in% names(data))) {
    log_warn(sprintf("亚组 [%s]: 列 '%s' 不存在，跳过", label, var_name))
    return(data[0])
  }
  # 类型安全检查
  col_data <- data[[var_name]]
  if (is.factor(col_data) || is.character(col_data)) {
    data <- copy(data)
    data[, (var_name) := as.numeric(as.character(get(var_name)))]
  }
  # 构建筛选表达式
  expr_str <- paste0(var_name, " ", condition)
  result <- tryCatch({
    idx <- which(with(as.list(data), eval(parse(text = expr_str))))
    data[idx]
  }, error = function(e) {
    log_warn(sprintf("亚组筛选失败 [%s]: %s", label, e$message))
    return(data[0])
  })
  if (nrow(result) > 0) {
    log_info(sprintf("  亚组 [%s]: %d 观测入选 (条件: %s)", label, nrow(result), expr_str))
  } else {
    log_warn(sprintf("  亚组 [%s]: 0 观测入选", label))
  }
  return(result)
}

# ---- 统一亚组SEQTaRget CCCW分析函数 ----
# 参数:
#   data         : 完整分析数据集
#   subgroup_data: 亚组数据子集
#   label        : 亚组标签
#   tv_cols      : 时变协变量
#   fx_cols      : 基线协变量
#   boot_n       : bootstrap 次数
#   fallback     : 是否启用降级策略（简化协变量/无bootstrap）
# 返回:
#   data.table 含 HR, CI, n, events, converged 状态
run_subgroup_seqtarget <- function(data, subgroup_data, label,
                                    tv_cols, fx_cols, boot_n = 50,
                                    fallback = TRUE) {
  n_sub <- nrow(subgroup_data)
  n_events <- sum(subgroup_data$outcome, na.rm = TRUE)
  
  if (n_sub < SUBGROUP_MIN_SAMPLE) {
    log_warn(sprintf("  [%s]: 样本量 %d < %d，跳过", label, n_sub, SUBGROUP_MIN_SAMPLE))
    return(data.table(analysis = label, hr = NA, ci_lower = NA, ci_upper = NA,
      n = n_sub, events = n_events, ok = FALSE, converged = FALSE,
      note = sprintf("样本量不足 (%d < %d)", n_sub, SUBGROUP_MIN_SAMPLE)))
  }
  
  opts <- SEQopts(hazard = TRUE, weighted = TRUE, weight.preexpansion = TRUE,
    weight.p99 = TRUE, bootstrap = TRUE, bootstrap.nboot = boot_n, nthreads = SEQTAARGET_NTHREADS)
  
  # 主尝试
  result <- tryCatch({
    m <- SEQuential(data = subgroup_data,
      id.col = "ID", time.col = "period",
      eligible.col = "eligible", treatment.col = "tx_init",
      outcome.col = "outcome",
      time_varying.cols = tv_cols, fixed.cols = fx_cols,
      method = "censoring", options = opts)
    hr <- hazard_ratio(m)
    ci_width <- if (!is.na(hr[3]) && !is.na(hr[2]) && hr[2] > 0) hr[3] / hr[2] else NA
    converged <- !is.na(ci_width) && ci_width < 50
    data.table(analysis = label, hr = hr[1], ci_lower = hr[2], ci_upper = hr[3],
      n = n_sub, events = n_events, ok = TRUE, converged = converged,
      note = if (!converged) sprintf("CI宽度比=%.1f", ci_width) else "ok")
  }, error = function(e) {
    log_warn(sprintf("  [%s]: SEQTaRget 失败 (%s)", label, e$message))
    NULL
  })
  
  if (!is.null(result) && result$ok && result$converged) return(result)
  
  # 降级策略1：简化协变量
  if (fallback && (is.null(result) || !result$converged)) {
    log_info(sprintf("  [%s]: 降级策略1 — 简化协变量集", label))
    reduced_tv <- intersect(c("baseline_phosphate"), tv_cols)
    reduced_fx <- intersect(c("age", "baseline_sofa", "sepsis3"), fx_cols)
    result <- tryCatch({
      m <- SEQuential(data = subgroup_data,
        id.col = "ID", time.col = "period",
        eligible.col = "eligible", treatment.col = "tx_init",
        outcome.col = "outcome",
        time_varying.cols = reduced_tv, fixed.cols = reduced_fx,
        method = "censoring", options = opts)
      hr <- hazard_ratio(m)
      ci_width <- if (!is.na(hr[3]) && !is.na(hr[2]) && hr[2] > 0) hr[3] / hr[2] else NA
      converged <- !is.na(ci_width) && ci_width < 50
      data.table(analysis = label, hr = hr[1], ci_lower = hr[2], ci_upper = hr[3],
        n = n_sub, events = n_events, ok = TRUE, converged = converged,
        note = sprintf("简化协变量(%dT+%dF)", length(reduced_tv), length(reduced_fx)))
    }, error = function(e2) {
      log_warn(sprintf("  [%s]: 降级策略1也失败 (%s)", label, e2$message))
      NULL
    })
  }
  
  # 降级策略2：关闭 bootstrap
  if (fallback && (is.null(result) || !result$converged)) {
    log_info(sprintf("  [%s]: 降级策略2 — 关闭bootstrap", label))
    opts_noboot <- modifyList(opts, list(bootstrap = FALSE))
    result <- tryCatch({
      m <- SEQuential(data = subgroup_data,
        id.col = "ID", time.col = "period",
        eligible.col = "eligible", treatment.col = "tx_init",
        outcome.col = "outcome",
        time_varying.cols = tv_cols, fixed.cols = fx_cols,
        method = "censoring", options = do.call(SEQopts, opts_noboot))
      hr <- hazard_ratio(m)
      data.table(analysis = label, hr = hr[1], ci_lower = hr[2], ci_upper = hr[3],
        n = n_sub, events = n_events, ok = TRUE, converged = TRUE,
        note = "无bootstrap(点估计)")
    }, error = function(e3) {
      NULL
    })
  }
  
  if (is.null(result)) {
    return(data.table(analysis = label, hr = NA, ci_lower = NA, ci_upper = NA,
      n = n_sub, events = n_events, ok = FALSE, converged = FALSE,
      note = "所有策略失败"))
  }
  return(result)
}

# ---- 正式交互效应检验（个体层面模型） ----
# 在合并数据中运行 Cox 模型，加入 treatment × covariate 交互项
# 参数:
#   data         : 分析数据集
#   interact_var : 交互变量名
#   covariate    : 协变量集
#   weights      : IPTW 权重向量
# 返回:
#   list(hr_interaction, p_interaction, model_summary)
test_interaction <- function(data, interact_var, covariates, weights = NULL) {
  if (!(interact_var %in% names(data))) {
    return(list(hr = NA, p_value = NA, note = "变量不存在"))
  }
  
  covs <- intersect(covariates, names(data))
  covs <- setdiff(covs, interact_var)
  
  # 构建交互项公式
  formula_str <- paste0("Surv(surv_time, outcome) ~ tx_init * ", interact_var)
  if (length(covs) > 0) {
    formula_str <- paste0(formula_str, " + ", paste(covs, collapse = " + "))
  }
  formula_str <- paste0(formula_str, " + cluster(ID)")
  
  fit <- tryCatch({
    if (is.null(weights)) {
      coxph(as.formula(formula_str), data = data)
    } else {
      coxph(as.formula(formula_str), data = data, weights = weights)
    }
  }, error = function(e) {
    log_warn(sprintf("交互检验失败 [%s]: %s", interact_var, e$message))
    return(NULL)
  })
  
  if (is.null(fit)) {
    return(list(hr = NA, p_value = NA, note = "模型拟合失败"))
  }
  
  # 提取交互项系数
  coef_names <- names(coef(fit))
  interact_coef <- grep(paste0("tx_init:", interact_var, "|", interact_var, ":tx_init"),
                        coef_names, value = TRUE)
  if (length(interact_coef) == 0) {
    return(list(hr = NA, p_value = NA, note = "未找到交互项"))
  }
  
  coef_idx <- which(coef_names == interact_coef[1])
  hr_interact <- exp(coef(fit)[coef_idx])
  p_interact <- summary(fit)$coefficients[coef_idx, "Pr(>|z|)"]
  
  log_info(sprintf("交互检验 [%s]: HR=%.3f, p=%.4f", interact_var, hr_interact, p_interact))
  
  list(
    hr = hr_interact,
    p_value = p_interact,
    coef_name = interact_coef[1],
    model = fit,
    note = if (p_interact < 0.05) "交互效应显著" else "交互效应不显著"
  )
}

# ---- 异质性检验（meta-analysis 框架） ----
# 对分层分析结果进行 Cochran's Q 检验和 I² 计算
# 参数:
#   strata_results: data.frame 含 hr, ci_low, ci_up, label 列
#   method        : meta-analysis 方法
# 返回:
#   list(Q, Q_p, I2, pooled_hr, pooled_ci)
pool_heterogeneity <- function(strata_results, method = "DL") {
  valid <- strata_results[!is.na(strata_results$hr) & !is.na(strata_results$ci_low) &
                          !is.na(strata_results$ci_up) & strata_results$ci_low > 0, ]
  
  if (nrow(valid) < HETEROGENEITY_MIN_STRATA) {
    log_warn(sprintf("有效分层数 %d < %d，无法进行异质性检验",
      nrow(valid), HETEROGENEITY_MIN_STRATA))
    return(list(Q = NA, Q_p = NA, I2 = NA, pooled_hr = NA, pooled_ci_low = NA,
                pooled_ci_up = NA, n_strata = nrow(valid)))
  }
  
  valid$log_hr <- log(valid$hr)
  valid$se <- (log(valid$ci_up) - log(valid$ci_low)) / (2 * 1.96)
  
  meta_res <- rma(yi = log_hr, sei = se, data = valid, method = method)
  
  log_info(sprintf("异质性检验: Q=%.2f, df=%d, p=%.4f, I²=%.1f%%",
    meta_res$QE, meta_res$k - 1, meta_res$QEp, meta_res$I2))
  
  if (meta_res$I2 > HETEROGENEITY_I2_WARN) {
    log_info(sprintf("  注意: I² > %d%%，存在显著异质性", HETEROGENEITY_I2_WARN))
  }
  
  list(
    Q = meta_res$QE, Q_p = meta_res$QEp, I2 = meta_res$I2,
    pooled_hr = exp(meta_res$b),
    pooled_ci_low = exp(meta_res$ci.lb),
    pooled_ci_up = exp(meta_res$ci.ub),
    n_strata = nrow(valid),
    tau2 = meta_res$tau2,
    meta_obj = meta_res
  )
}

# ---- 剂量-反应关系分析（RCS 限制性立方样条） ----
# 评估血磷水平与治疗效应的连续剂量-反应关系
# 参数:
#   data       : 分析数据集
#   dose_var   : 剂量变量名（如 baseline_phosphate）
#   covariates : 协变量集
#   knots      : RCS 节点数
#   ref_value  : 参考值
# 返回:
#   list 含 RCS 拟合对象和预测值
run_dose_response <- function(data, dose_var = "baseline_phosphate",
                               covariates, knots = 3, ref_value = 0.65) {
  if (!requireNamespace("rms", quietly = TRUE)) {
    log_warn("rms 包未安装，无法运行剂量-反应分析")
    return(NULL)
  }
  
  covs <- intersect(covariates, names(data))
  covs <- setdiff(covs, dose_var)
  
  # 构建 RCS 公式
  rcs_formula <- paste0("Surv(surv_time, outcome) ~ rcs(", dose_var, ", ", knots, ") * tx_init")
  if (length(covs) > 0) {
    rcs_formula <- paste0(rcs_formula, " + ", paste(covs, collapse = " + "))
  }
  
  fit <- tryCatch({
    cph(as.formula(rcs_formula), data = data, x = TRUE, y = TRUE)
  }, error = function(e) {
    log_warn(sprintf("剂量-反应模型拟合失败: %s", e$message))
    return(NULL)
  })
  
  if (is.null(fit)) return(NULL)
  
  # 预测血磷值范围内的 HR
  phos_range <- quantile(data[[dose_var]], c(0.05, 0.95), na.rm = TRUE)
  phos_seq <- seq(phos_range[1], phos_range[2], length.out = 100)
  
  # 对治疗组和对照组分别预测
  pred_data <- data.frame(
    baseline_phosphate = rep(phos_seq, 2),
    tx_init = rep(c(0, 1), each = length(phos_seq))
  )
  # 填充协变量中位数
  for (cv in covs) {
    if (cv %in% names(data)) {
      pred_data[[cv]] <- median(data[[cv]], na.rm = TRUE)
    }
  }
  
  pred <- tryCatch({
    predict(fit, newdata = pred_data, type = "lp")
  }, error = function(e) {
    log_warn(sprintf("剂量-反应预测失败: %s", e$message))
    return(NULL)
  })
  
  list(
    model = fit,
    dose_var = dose_var,
    dose_range = phos_range,
    dose_seq = phos_seq,
    predictions = pred,
    ref_value = ref_value,
    knots = knots
  )
}

# ---- 多重比较校正 ----
# 参数:
#   p_values : p 值向量
#   method   : 校正方法
# 返回:
#   校正后的 p 值向量
adjust_pvalues <- function(p_values, method = "holm") {
  if (method == "none") return(p_values)
  p.adjust(p_values, method = method)
}

# ---- 亚组分析森林图生成函数 ----
# 参数:
#   results   : data.table 含 analysis, hr, ci_lower, ci_upper, p_value, n, events
#   title     : 图标题
#   filename  : 输出文件名（不含扩展名）
# 返回:
#   无（生成 PDF + PNG 双格式）
generate_forest_plot <- function(results, title = "亚组分析森林图",
                                  filename = "v25_forest_subgroup") {
  # 过滤有效结果
  plot_data <- results[ok == TRUE & !is.na(hr) & !is.na(ci_lower) & !is.na(ci_upper)]
  if (nrow(plot_data) == 0) {
    log_warn("无有效结果，跳过森林图生成")
    return(invisible(NULL))
  }
  
  plot_data[, idx := .I]
  # 颜色：HR>1 红色，HR<1 蓝色
  plot_data[, col := fifelse(hr > 1, JOURNAL_COLORS$treatment, JOURNAL_COLORS$control)]
  
  x_limits <- range(c(plot_data$ci_lower, plot_data$ci_upper), na.rm = TRUE)
  x_limits <- c(min(x_limits, 0.3), max(x_limits, 3.0))
  
  # PDF 输出
  pdf(file.path(FIG_DIR, paste0(filename, ".pdf")), width = 12, height = 2 + nrow(plot_data) * 0.4)
  par(mar = c(4, 18, 3, 2), family = "sans")
  
  plot(NA, NA, xlim = x_limits, ylim = c(0.5, nrow(plot_data) + 0.5),
       xlab = "Hazard Ratio (95% CI)", ylab = "", yaxt = "n",
       main = "", log = "x")
  title(title, adj = 0, cex.main = 1.1)
  
  abline(v = 1, lty = 2, col = JOURNAL_COLORS$ref_line, lwd = 1.5)
  
  for (i in 1:nrow(plot_data)) {
    y_pos <- nrow(plot_data) - i + 1
    hr_val <- plot_data$hr[i]
    ci_l <- plot_data$ci_lower[i]
    ci_u <- plot_data$ci_upper[i]
    
    # CI 线
    lines(c(ci_l, ci_u), c(y_pos, y_pos), lwd = 2, col = plot_data$col[i])
    # 点估计
    points(hr_val, y_pos, pch = 22, bg = plot_data$col[i], cex = 1.5)
    # 标签
    mtext(plot_data$analysis[i], side = 2, at = y_pos, las = 1, cex = 0.8,
          adj = 1, line = 0.5)
  }
  
  dev.off()
  log_info(sprintf("森林图 PDF 已保存: %s", file.path(FIG_DIR, paste0(filename, ".pdf"))))
  
  # PNG 输出
  png(file.path(FIG_DIR, paste0(filename, ".png")), width = 12, height = 2 + nrow(plot_data) * 0.4,
      units = "in", res = 300)
  par(mar = c(4, 18, 3, 2), family = "sans")
  plot(NA, NA, xlim = x_limits, ylim = c(0.5, nrow(plot_data) + 0.5),
       xlab = "Hazard Ratio (95% CI)", ylab = "", yaxt = "n",
       main = "", log = "x")
  title(title, adj = 0, cex.main = 1.1)
  abline(v = 1, lty = 2, col = JOURNAL_COLORS$ref_line, lwd = 1.5)
  for (i in 1:nrow(plot_data)) {
    y_pos <- nrow(plot_data) - i + 1
    hr_val <- plot_data$hr[i]
    ci_l <- plot_data$ci_lower[i]
    ci_u <- plot_data$ci_upper[i]
    lines(c(ci_l, ci_u), c(y_pos, y_pos), lwd = 2, col = plot_data$col[i])
    points(hr_val, y_pos, pch = 22, bg = plot_data$col[i], cex = 1.5)
    mtext(plot_data$analysis[i], side = 2, at = y_pos, las = 1, cex = 0.8,
          adj = 1, line = 0.5)
  }
  dev.off()
  log_info(sprintf("森林图 PNG 已保存: %s", file.path(FIG_DIR, paste0(filename, ".png"))))
  
  invisible(NULL)
}

# ===========================================================================
# 2026-06 重构模块加载完成
# ===========================================================================
log_info("_helpers_v25.R 重构模块加载完成（含亚组分析、异质性检验、交互效应、剂量-反应、森林图）")

# ===========================================================================
# 🆕 P0-CP01: RMST 分析（限制性平均生存时间）— Co-primary 估计量
# 用于 PH 假设违反时的互补估计
# 参考: Royston & Parmar (2013), BMC Med Res Methodol
#       Uno et al. (2014), Stat Med
# ===========================================================================

# 计算 RMST 差值
# 参数:
#   data    : data.table 含 surv_time, outcome (0/1), tx_init (0/1)
#   tau     : 截断时间（天），默认 28
#   weights : IPTW 权重向量（可选）
# 返回:
#   list 含 RMST 差值、SE、CI、RMST 比率
compute_rmst <- function(data, tau = 28, weights = NULL) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    log_error("survival 包未安装，无法计算 RMST")
    return(NULL)
  }
  
  # 构建 Surv 对象
  surv_obj <- Surv(data$surv_time, data$outcome)
  
  # 拟合 KM（按处理组）
  fit <- survfit(surv_obj ~ data$tx_init, weights = weights)
  
  # 提取 RMST
  rmst_result <- summary(fit, rmean = tau)
  rmst_tbl <- rmst_result$table
  
  if (nrow(rmst_tbl) < 2) {
    log_error("RMST 计算失败：需要两个处理组")
    return(NULL)
  }
  
  # 动态获取列名（survival 包版本可能不同）
  col_names <- colnames(rmst_tbl)
  rmean_col <- grep("rmean", col_names, value = TRUE)[1]
  se_col <- grep("se\\(rmean\\)|std\\.err|rmean\\.se", col_names, value = TRUE, ignore.case = TRUE)[1]
  records_col <- grep("records", col_names, value = TRUE)[1]
  events_col <- grep("events", col_names, value = TRUE)[1]
  
  # 如果找不到 se 列，尝试从 rmean CI 反算
  if (is.na(se_col)) {
    ci_lower_col <- grep("0\\.95LCL|lower", col_names, value = TRUE)[1]
    ci_upper_col <- grep("0\\.95UCL|upper", col_names, value = TRUE)[1]
    if (!is.na(ci_lower_col) && !is.na(ci_upper_col)) {
      rmst_tbl <- cbind(rmst_tbl, 
        se_est = (rmst_tbl[, ci_upper_col] - rmst_tbl[, ci_lower_col]) / (2 * 1.96))
      se_col <- "se_est"
    }
  }
  
  if (is.na(se_col) || is.na(rmean_col)) {
    log_error(sprintf("RMST 表列名: %s", paste(col_names, collapse = ", ")))
    return(NULL)
  }
  
  rmst_treat <- rmst_tbl[2, rmean_col]
  rmst_ctrl  <- rmst_tbl[1, rmean_col]
  se_treat   <- rmst_tbl[2, se_col]
  se_ctrl    <- rmst_tbl[1, se_col]
  
  rmst_diff <- rmst_treat - rmst_ctrl
  rmst_se <- sqrt(se_treat^2 + se_ctrl^2)
  
  list(
    rmst_treatment = rmst_treat,
    rmst_control   = rmst_ctrl,
    rmst_diff      = rmst_diff,
    rmst_se        = rmst_se,
    rmst_ci_lower  = rmst_diff - 1.96 * rmst_se,
    rmst_ci_upper  = rmst_diff + 1.96 * rmst_se,
    rmst_ratio     = rmst_treat / rmst_ctrl,
    tau            = tau,
    n_treat        = if (!is.na(records_col)) rmst_tbl[2, records_col] else NA_integer_,
    n_ctrl         = if (!is.na(records_col)) rmst_tbl[1, records_col] else NA_integer_,
    events_treat   = if (!is.na(events_col)) rmst_tbl[2, events_col] else NA_integer_,
    events_ctrl    = if (!is.na(events_col)) rmst_tbl[1, events_col] else NA_integer_
  )
}

# RMST 的 Bootstrap 置信区间
# 参数:
#   data      : data.table
#   tau       : 截断时间
#   weights_fn: 权重计算函数（返回权重向量）
#   n_boot    : bootstrap 次数
# 返回:
#   data.table 含 RMST 差值、Bootstrap CI
compute_rmst_bootstrap <- function(data, tau = 28, weights_fn = NULL, n_boot = 200) {
  n <- nrow(data)
  rmst_boot <- numeric(n_boot)
  
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    boot_data <- data[idx]
    w <- if (!is.null(weights_fn)) weights_fn(boot_data) else NULL
    rmst_result <- tryCatch({
      compute_rmst(boot_data, tau = tau, weights = w)
    }, error = function(e) NULL)
    rmst_boot[i] <- if (!is.null(rmst_result)) rmst_result$rmst_diff else NA
  }
  
  rmst_boot <- na.omit(rmst_boot)
  ci <- quantile(rmst_boot, c(0.025, 0.975))
  
  data.table(
    rmst_diff_boot_mean = mean(rmst_boot),
    rmst_diff_boot_se   = sd(rmst_boot),
    rmst_diff_boot_ci_l = ci[1],
    rmst_diff_boot_ci_u = ci[2],
    n_valid_boot        = length(rmst_boot)
  )
}

# ===========================================================================
# 🆕 P0-CP02: Weibull AFT 模型 — PH 违反时的参数化互补估计
# 参考: Wei (1992), Biometrika; Cleves et al. (2010), Stata
# ===========================================================================

# 运行 Weibull AFT 模型
# 参数:
#   data       : data.table
#   covariates : 协变量向量
#   weights    : IPTW 权重向量（可选）
# 返回:
#   list 含 HR、CI、scale 参数
compute_aft <- function(data, covariates, weights = NULL) {
  covs <- intersect(covariates, names(data))
  formula_str <- paste0("Surv(surv_time, outcome) ~ tx_init + ",
                        paste(covs, collapse = " + "))
  
  fit <- tryCatch({
    if (is.null(weights)) {
      survreg(as.formula(formula_str), data = data, dist = "weibull")
    } else {
      survreg(as.formula(formula_str), data = data,
              weights = weights, dist = "weibull")
    }
  }, error = function(e) {
    log_warn(sprintf("Weibull AFT 拟合失败: %s", e$message))
    return(NULL)
  })
  
  if (is.null(fit)) return(NULL)
  
  # AFT 系数转 HR: HR = exp(-coef / scale)
  scale_param <- fit$scale
  aft_coef <- -coef(fit)["tx_init"] / scale_param
  aft_hr <- exp(aft_coef)
  aft_se <- sqrt(vcov(fit)["tx_init", "tx_init"]) / scale_param
  
  list(
    hr = aft_hr,
    ci_lower = exp(aft_coef - 1.96 * aft_se),
    ci_upper = exp(aft_coef + 1.96 * aft_se),
    scale = scale_param,
    loglik = fit$loglik[2],
    converged = fit$convergence == 0,
    n_covariates = length(covs)
  )
}

# ===========================================================================
# 🆕 P1-ST01: 重叠权重（Overlap Weighting）函数
# 参考: Li, Morgan & Zaslavsky (2018), JASA
# 重叠权重将目标人群聚焦于"临床均衡区"（PS 接近 0.5），
# 对极端 PS 值的观测给予低权重，显著减小方差
# ===========================================================================

# 计算重叠权重
# 参数:
#   ps    : 倾向得分向量
#   treat : 处理指示向量 (0/1)
# 返回:
#   重叠权重向量
compute_overlap_weights <- function(ps, treat) {
  if (length(ps) != length(treat)) {
    stop("ps 和 treat 长度不匹配")
  }
  ow <- numeric(length(ps))
  ow[treat == 1] <- 1 - ps[treat == 1]
  ow[treat == 0] <- ps[treat == 0]
  # 避免极端权重
  ow <- pmax(ow, 1e-10)
  return(ow)
}

# ===========================================================================
# 🆕 P1-ST02: 匹配权重（Matching Weight）函数
# 参考: Li & Greene (2013), Stat Med
# 匹配权重模拟 1:1 匹配，将 PS 分布映射到两组重叠区域
# ===========================================================================

# 计算匹配权重
# 参数:
#   ps    : 倾向得分向量
#   treat : 处理指示向量 (0/1)
# 返回:
#   匹配权重向量
compute_matching_weights <- function(ps, treat) {
  if (length(ps) != length(treat)) {
    stop("ps 和 treat 长度不匹配")
  }
  mw <- numeric(length(ps))
  mw[treat == 1] <- pmin(ps[treat == 1], 1 - ps[treat == 1]) / ps[treat == 1]
  mw[treat == 0] <- pmin(ps[treat == 0], 1 - ps[treat == 0]) / (1 - ps[treat == 0])
  # 截断极端值
  mw <- pmin(mw, quantile(mw, 0.99, na.rm = TRUE))
  mw <- pmax(mw, 1e-10)
  return(mw)
}

# ===========================================================================
# 🆕 P1-ST03: 权重方法对比函数
# 对比三种权重方法的诊断指标
# ===========================================================================

# 权重诊断对比
# 参数:
#   data       : data.table 含 treatment 和协变量
#   covariates : 协变量向量
# 返回:
#   data.table 含三种权重方法的 ESS、max_weight、SMD_before/after
compare_weight_methods <- function(data, covariates) {
  covs <- intersect(covariates, names(data))
  if (length(covs) < 2) {
    log_warn("协变量不足 (需 >=2)，无法进行权重对比")
    return(NULL)
  }
  
  # 估计 PS
  ps_formula <- as.formula(paste("tx_init ~", paste(covs, collapse = " + ")))
  ps_glm <- glm(ps_formula, data = data, family = binomial(link = "logit"))
  ps <- predict(ps_glm, type = "response")
  treat <- data$tx_init
  
  # 计算三种权重
  w_iptw <- compute_iptw(data, covs, trunc_p = 0.99, lower_p = 0.01,
                         stabilize = TRUE, bilateral = TRUE)
  w_overlap <- compute_overlap_weights(ps, treat)
  w_matching <- compute_matching_weights(ps, treat)
  
  # 诊断对比
  calc_diagnostics <- function(w, label) {
    ess <- (sum(w)^2) / sum(w^2)
    ess_ratio <- ess / length(w)
    max_w <- max(w, na.rm = TRUE)
    mean_w <- mean(w, na.rm = TRUE)
    sd_w <- sd(w, na.rm = TRUE)
    
    # 平均 SMD（加权后）
    smd_values <- sapply(covs, function(v) {
      x <- data[[v]]
      if (!is.numeric(x)) return(NA_real_)
      tryCatch(calc_smd(x, treat, w = w), error = function(e) NA_real_)
    })
    mean_smd <- mean(abs(na.omit(smd_values)))
    n_smd_lt_0.1 <- sum(abs(na.omit(smd_values)) < 0.1, na.rm = TRUE)
    n_smd_lt_0.2 <- sum(abs(na.omit(smd_values)) < 0.2, na.rm = TRUE)
    
    data.table(
      method = label,
      n_original = length(w),
      ess = round(ess, 1),
      ess_ratio = round(ess_ratio * 100, 1),
      max_weight = round(max_w, 2),
      mean_weight = round(mean_w, 3),
      sd_weight = round(sd_w, 3),
      mean_smd = round(mean_smd, 4),
      n_smd_lt_0.1 = n_smd_lt_0.1,
      n_smd_lt_0.2 = n_smd_lt_0.2,
      n_vars = length(smd_values)
    )
  }
  
  diag_iptw <- calc_diagnostics(w_iptw, "IPTW (稳定化)")
  diag_overlap <- calc_diagnostics(w_overlap, "Overlap Weighting")
  diag_matching <- calc_diagnostics(w_matching, "Matching Weight")
  
  result <- rbindlist(list(diag_iptw, diag_overlap, diag_matching))
  log_info("权重方法对比完成")
  return(result)
}

# ===========================================================================
# 2026-06 V25.5 RMST/AFT/权重对比模块加载完成
# ===========================================================================
log_info("_helpers_v25.R 增强模块加载完成（含RMST、AFT、重叠权重、匹配权重、权重对比）")

# ===========================================================================
# 🆕 P0-M03: Cluster-Robust Standard Error (Sandwich Estimator)
# 处理重复入组（同一 stay_id 多个 trial period）的 within-person correlation
# 参考: Williams (2000), Biometrika; Lee, Wei, & Amato (1992), Stat Med
# ===========================================================================

# 计算 cluster-robust 方差-协方差矩阵
# 参数:
#   cox_fit   : coxph() 拟合对象
#   cluster_id: 聚类标识向量（如 stay_id）
# 返回:
#   cluster-robust 协方差矩阵
compute_cluster_robust_vcov <- function(cox_fit, cluster_id) {
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    log_warn("sandwich 包未安装，跳过 cluster-robust SE")
    return(NULL)
  }
  
  # 计算 score residuals（按个体聚类）
  cluster_id <- as.factor(cluster_id)
  
  # 使用 sandwich::vcovCL 计算 cluster-robust 协方差
  vcov_cl <- tryCatch({
    sandwich::vcovCL(cox_fit, cluster = cluster_id, type = "HC0")
  }, error = function(e) {
    log_warn(sprintf("vcovCL 失败: %s，尝试 vcovCL(type='HC1')", e$message))
    tryCatch({
      sandwich::vcovCL(cox_fit, cluster = cluster_id, type = "HC1")
    }, error = function(e2) {
      log_warn(sprintf("vcovCL(HC1) 也失败: %s", e2$message))
      return(NULL)
    })
  })
  
  return(vcov_cl)
}

# 比较 cluster-robust SE 与 naive SE
# 参数:
#   cox_fit   : coxph() 拟合对象
#   cluster_id: 聚类标识向量
# 返回:
#   data.table 含 SE 对比和膨胀因子
compare_cluster_se <- function(cox_fit, cluster_id) {
  vcov_cl <- compute_cluster_robust_vcov(cox_fit, cluster_id)
  
  if (is.null(vcov_cl)) {
    return(data.table(
      coef = names(coef(cox_fit)),
      naive_se = sqrt(diag(vcov(cox_fit))),
      cluster_se = NA_real_,
      inflation_factor = NA_real_,
      ci_width_change_pct = NA_real_
    ))
  }
  
  naive_se <- sqrt(diag(vcov(cox_fit)))
  cluster_se <- sqrt(diag(vcov_cl))
  
  result <- data.table(
    coef = names(coef(cox_fit)),
    naive_se = round(naive_se, 6),
    cluster_se = round(cluster_se, 6),
    inflation_factor = round(cluster_se / naive_se, 3),
    ci_width_change_pct = round((cluster_se / naive_se - 1) * 100, 1)
  )
  
  # 标记处理组变量
  result[coef == "tx_init", coef := "tx_init (treatment)"]
  
  return(result)
}

# 首次 trial period 敏感性分析（仅使用 enrollment_seq=1 的观测）
# 用于验证重复入组是否影响结果
run_first_period_sensitivity <- function(data, tv_cols, fx_cols, boot_n = 50) {
  # 仅保留首次入组
  if ("period" %in% names(data)) {
    first_period_data <- data[period == 0]
  } else if ("enrollment_seq" %in% names(data)) {
    first_period_data <- data[enrollment_seq == 1]
  } else {
    log_warn("无法识别首次入组（缺少 period/enrollment_seq），跳过首次期敏感性分析")
    return(NULL)
  }
  
  n_unique <- uniqueN(first_period_data$ID)
  n_events <- sum(first_period_data$outcome)
  log_info(sprintf("首次期敏感性分析: N=%d (唯一患者), events=%d", n_unique, n_events))
  
  # 运行 PP IPTW 分析
  opts <- SEQopts(hazard = TRUE, weighted = TRUE, weight.preexpansion = TRUE,
    weight.p99 = TRUE, bootstrap = TRUE, bootstrap.nboot = boot_n, nthreads = SEQTAARGET_NTHREADS)
  
  result <- tryCatch({
    m <- SEQuential(data = first_period_data,
      id.col = "ID", time.col = "period",
      eligible.col = "eligible", treatment.col = "tx_init",
      outcome.col = "outcome",
      time_varying.cols = tv_cols, fixed.cols = fx_cols,
      method = "censoring", options = opts)
    hr <- hazard_ratio(m)
    data.table(
      analysis = "First_Period_Only",
      hr = hr[1], ci_lower = hr[2], ci_upper = hr[3],
      n_patients = n_unique, n_events = n_events,
      note = "仅首次入组，消除重复入组相关性"
    )
  }, error = function(e) {
    data.table(
      analysis = "First_Period_Only",
      hr = NA, ci_lower = NA, ci_upper = NA,
      n_patients = n_unique, n_events = n_events,
      note = sprintf("失败: %s", e$message)
    )
  })
  
  return(result)
}

# ===========================================================================
# 🆕 P1-S01: PS模型诊断函数（C-statistic + 校准度 + 重叠度）
# 参考: Austin & Steyerberg (2015), Stat Med; Dehejia & Wahba (2002), REStat
# ===========================================================================

# 计算 PS 模型的 C-statistic（判别度）
compute_ps_cstat <- function(ps_glm) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    log_warn("pROC 包未安装，无法计算 C-statistic")
    return(NA_real_)
  }
  predicted <- predict(ps_glm, type = "response")
  actual <- ps_glm$y
  tryCatch({
    roc_obj <- pROC::roc(actual, predicted, quiet = TRUE)
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) {
    log_warn(sprintf("C-statistic 计算失败: %s", e$message))
    NA_real_
  })
}

# 计算 PS 模型的 Hosmer-Lemeshow 校准度检验
compute_hl_calibration <- function(ps_glm, g = 10) {
  if (!requireNamespace("ResourceSelection", quietly = TRUE)) {
    log_warn("ResourceSelection 包未安装，无法计算 HL 校准度")
    return(data.table(hl_stat = NA_real_, hl_df = NA_integer_, hl_p = NA_real_))
  }
  predicted <- predict(ps_glm, type = "response")
  actual <- ps_glm$y
  tryCatch({
    hl <- ResourceSelection::hoslem.test(actual, predicted, g = g)
    data.table(hl_stat = hl$statistic, hl_df = hl$parameter, hl_p = hl$p.value)
  }, error = function(e) {
    log_warn(sprintf("HL 检验失败: %s", e$message))
    data.table(hl_stat = NA_real_, hl_df = NA_integer_, hl_p = NA_real_)
  })
}

# PS 分布重叠度评估
assess_ps_overlap <- function(ps, treat) {
  ps_treat <- ps[treat == 1]
  ps_ctrl  <- ps[treat == 0]
  
  min_treat <- min(ps_treat, na.rm = TRUE)
  max_treat <- max(ps_treat, na.rm = TRUE)
  min_ctrl  <- min(ps_ctrl, na.rm = TRUE)
  max_ctrl  <- max(ps_ctrl, na.rm = TRUE)
  
  csr_lower <- max(min_treat, min_ctrl)
  csr_upper <- min(max_treat, max_ctrl)
  
  treat_in_csr <- mean(ps_treat >= csr_lower & ps_treat <= csr_upper, na.rm = TRUE)
  ctrl_in_csr  <- mean(ps_ctrl >= csr_lower & ps_ctrl <= csr_upper, na.rm = TRUE)
  
  if (treat_in_csr < 0.90 || ctrl_in_csr < 0.90) {
    diagnosis <- "警告：PS重叠度不足，建议考虑截断或子集分析"
  } else if (treat_in_csr < 0.95 || ctrl_in_csr < 0.95) {
    diagnosis <- "PS重叠度可接受，但边缘区域需关注"
  } else {
    diagnosis <- "PS重叠度良好"
  }
  
  list(
    csr = c(csr_lower, csr_upper),
    treat_overlap_pct = round(treat_in_csr * 100, 1),
    ctrl_overlap_pct  = round(ctrl_in_csr * 100, 1),
    diagnosis = diagnosis,
    ps_range_treat = c(min_treat, max_treat),
    ps_range_ctrl  = c(min_ctrl, max_ctrl)
  )
}

# 完整 PS 模型诊断报告
run_ps_diagnostics <- function(data, covariates) {
  covs <- intersect(covariates, names(data))
  if (length(covs) < 2) {
    log_warn("协变量不足 (需 >=2)，无法运行 PS 诊断")
    return(list(cstat = NA, hl = NULL, overlap = NULL, vif = NULL))
  }
  
  ps_formula <- as.formula(paste("treatment ~", paste(covs, collapse = " + ")))
  ps_glm <- glm(ps_formula, data = data, family = binomial(link = "logit"))
  ps <- predict(ps_glm, type = "response")
  
  cstat <- compute_ps_cstat(ps_glm)
  log_info(sprintf("PS C-statistic: %.3f", cstat))
  if (!is.na(cstat) && cstat > 0.80) {
    log_info("  C-statistic > 0.80: PS模型判别度良好")
  } else if (!is.na(cstat) && cstat < 0.70) {
    log_info("  注意：C-statistic < 0.70，PS模型判别度偏低")
  }
  
  hl <- compute_hl_calibration(ps_glm)
  if (!is.na(hl$hl_p)) {
    log_info(sprintf("HL χ²=%.2f, df=%d, p=%.4f", hl$hl_stat, hl$hl_df, hl$hl_p))
    if (hl$hl_p < 0.05) {
      log_info("  HL p < 0.05: PS模型校准度不佳（建议使用GBM或npCBPS替代）")
    }
  }
  
  overlap <- assess_ps_overlap(ps, data$treatment)
  log_info(sprintf("PS 重叠度: 治疗组 %.1f%%, 对照组 %.1f%%",
    overlap$treat_overlap_pct, overlap$ctrl_overlap_pct))
  log_info(sprintf("  %s", overlap$diagnosis))
  
  vif <- tryCatch({
    if (requireNamespace("car", quietly = TRUE)) {
      vif_vals <- car::vif(ps_glm)
      high_vif <- names(vif_vals[vif_vals > 5])
      if (length(high_vif) > 0) {
        log_info(sprintf("VIF>5 变量 (%d): %s", length(high_vif),
          paste(high_vif, collapse = ", ")))
      }
      vif_vals
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  list(
    cstat = cstat,
    hl = hl,
    overlap = overlap,
    vif = vif,
    ps = ps,
    ps_glm = ps_glm
  )
}

# ===========================================================================
# 🆕 P1-S02: PH假设检验（Schoenfeld残差）
# 参考: Grambsch & Therneau (1994), Biometrika
# ===========================================================================

test_ph_assumption <- function(cox_fit) {
  schoenfeld_test <- cox.zph(cox_fit)
  
  global_chisq <- schoenfeld_test$table["GLOBAL", "chisq"]
  global_df <- schoenfeld_test$table["GLOBAL", "df"]
  global_p <- schoenfeld_test$table["GLOBAL", "p"]
  
  log_info(sprintf("PH 全局检验: χ²=%.2f, df=%d, p=%.4f", global_chisq, global_df, global_p))
  if (global_p < 0.05) {
    log_info("  警告：全局 PH 假设被拒绝（p<0.05），建议补充 RMST 分析")
  }
  
  var_names <- rownames(schoenfeld_test$table)
  var_names <- var_names[var_names != "GLOBAL"]
  
  ph_results <- data.table(
    variable = var_names,
    chisq = schoenfeld_test$table[var_names, "chisq"],
    df = schoenfeld_test$table[var_names, "df"],
    p_value = schoenfeld_test$table[var_names, "p"],
    ph_violated = ifelse(schoenfeld_test$table[var_names, "p"] < 0.05, "YES", "No")
  )
  
  ph_results[variable == "tx_init", variable := "tx_init (treatment)"]
  
  n_violated <- sum(ph_results$ph_violated == "YES")
  if (n_violated > 0) {
    violated_vars <- ph_results[ph_violated == "YES"]$variable
    log_info(sprintf("  PH违反变量 (%d): %s", n_violated, paste(violated_vars, collapse = ", ")))
  }
  
  list(
    global = data.table(chisq = global_chisq, df = global_df, p = global_p),
    by_variable = ph_results,
    schoenfeld_obj = schoenfeld_test
  )
}

# ===========================================================================
# 🆕 P0-M06: CCW删失报告函数
# ===========================================================================

report_ccw_censoring <- function(model_pp, data) {
  log_section("CCW 删失报告")
  
  rd <- tryCatch(risk_data(model_pp), error = function(e) NULL)
  
  if (is.null(rd)) {
    log_warn("无法从 SEQTaRget 模型提取风险数据，跳过 CCW 删失报告")
    return(NULL)
  }
  
  n_total <- nrow(data)
  n_censored <- sum(rd$status == 0, na.rm = TRUE)
  n_events <- sum(rd$status == 1, na.rm = TRUE)
  
  log_info(sprintf("CCCW 删失统计: 总观测=%d, 事件=%d, 删失=%d (%.1f%%)",
    n_total, n_events, n_censored, n_censored / n_total * 100))
  
  if ("period" %in% names(rd) || "time" %in% names(rd)) {
    time_col <- if ("period" %in% names(rd)) "period" else "time"
    period_censoring <- rd[, .(
      n_total = .N,
      n_censored = sum(status == 0),
      n_events = sum(status == 1)
    ), by = time_col]
    period_censoring[, censor_rate := round(n_censored / n_total * 100, 1)]
    
    log_info("按 period 的删失率:")
    print(period_censoring)
    
    high_censor_periods <- period_censoring[censor_rate > 20]
    if (nrow(high_censor_periods) > 0) {
      log_info(sprintf("  警告：%d 个 period 删失率>20%%", nrow(high_censor_periods)))
    }
  }
  
  list(
    n_total = n_total, n_events = n_events, n_censored = n_censored,
    censor_rate = round(n_censored / n_total * 100, 1),
    period_censoring = if (exists("period_censoring")) period_censoring else NULL
  )
}

# ===========================================================================
# 🆕 P0-M02: time_zero 验证与补磷时间分布报告
# ===========================================================================

report_phosphate_timing <- function(data) {
  log_section("补磷时间分布与 time_zero 验证")
  
  results <- list()
  
  has_first_iv <- "first_iv_time" %in% names(data)
  has_time_zero <- "time_zero" %in% names(data)
  has_iv_12h <- "iv_within_12h" %in% names(data)
  
  if (!has_first_iv || !has_time_zero) {
    log_warn("缺少 first_iv_time 或 time_zero 列，跳过补磷时间分布报告")
    return(NULL)
  }
  
  treat_data <- data[tx_init == 1]
  if (nrow(treat_data) == 0) {
    log_warn("无治疗组观测")
    return(NULL)
  }
  
  treat_data[, time_to_iv_hours := as.numeric(difftime(
    as.POSIXct(first_iv_time), as.POSIXct(time_zero), units = "hours"))]
  
  treat_data <- treat_data[!is.na(time_to_iv_hours) & time_to_iv_hours >= 0]
  
  n_valid <- nrow(treat_data)
  log_info(sprintf("有效补磷时间记录: %d/%d (%.1f%%)",
    n_valid, sum(data$tx_init == 1), n_valid / sum(data$tx_init == 1) * 100))
  
  time_quantiles <- quantile(treat_data$time_to_iv_hours, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  log_info(sprintf("补磷时间: 中位数 %.1fh, IQR %.1f-%.1fh",
    time_quantiles[2], time_quantiles[1], time_quantiles[3]))
  
  windows <- c(6, 12, 24, 48)
  for (w in windows) {
    n_within <- sum(treat_data$time_to_iv_hours <= w, na.rm = TRUE)
    log_info(sprintf("  %dh内补磷: %d (%.1f%%)", w, n_within, n_within / n_valid * 100))
  }
  
  if (has_iv_12h) {
    iv_12h_match <- treat_data[, .(
      n_total = .N,
      n_12h_true = sum(iv_within_12h == TRUE | iv_within_12h == 1, na.rm = TRUE),
      n_12h_calc = sum(time_to_iv_hours <= 12, na.rm = TRUE)
    )]
    log_info(sprintf("iv_within_12h 一致性: 标记=%d, 计算=%d",
      iv_12h_match$n_12h_true, iv_12h_match$n_12h_calc))
    results$iv_12h_validation <- iv_12h_match
  }
  
  results$time_quantiles <- time_quantiles
  results$time_distribution <- treat_data[, .(
    n = .N,
    mean_h = mean(time_to_iv_hours),
    median_h = median(time_to_iv_hours),
    sd_h = sd(time_to_iv_hours),
    within_6h = sum(time_to_iv_hours <= 6),
    within_12h = sum(time_to_iv_hours <= 12),
    within_24h = sum(time_to_iv_hours <= 24),
    within_48h = sum(time_to_iv_hours <= 48)
  )]
  
  return(results)
}