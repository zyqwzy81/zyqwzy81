# ============================================================================
# common.R
# 主分析公共模块：把"试验行 + MICE 插补"构建成 SEQTaRget 所需的
# 日级 person-period 表（每个 stay 一条连续日历时间线，自然日锚定）。
#
# 关键口径（代码审计 C1/C2 修复）：
#   * ID = stay_id；时间轴按"自然日历日期"锚定，day0 = 该 stay 最早起点日期
#   * eligible=1 落在每个试验起点所在自然日（起点不冲突）
#   * tx_init 在起点日按 itt_group 定臂（0/1），纯 ITT 随访（不删失、不切换）
#   * outcome 仅在死亡日落 1；死亡日用 stay 层共识（origin_day+round(surv_time) 众数）
#   * 网格延伸到 max(origin)+28，出院不删失（B 口径）
#   * DAG14 全部作为"起点基线"写入起点日并前向填充（包内 *_bas 机制）
# 版本: V28.4   日期: 2026-10-09
# ============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

# ---- 全局配置 ----
V28_BASE   <- "/data/小米云盘/文档/论文/血磷"
AUG_FILE   <- file.path(V28_BASE, "TTE/v26/output/sequential_trials_augmented.parquet")
MICE_FMT   <- file.path(V28_BASE, "sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
SEED       <- 20260629L
TAU        <- 28L

# DAG14 变量来源
DAG_FROM_IMP <- c("age", "baseline_sofa", "sepsis3", "charlson", "egfr",
                  "invasive_vent_at_t0", "vasoactive_at_t0", "baseline_phosphate",
                  "lactate_t0", "refeeding_risk", "corrected_calcium_mmol_min")
DAG_FROM_AUG <- c("gender", "admission_type", "nutrition_support_type")
DAG14        <- c(DAG_FROM_IMP, DAG_FROM_AUG)

# 分类变量水平（保证 model.matrix 跨 stay/插补一致）
GENDER_LEV   <- c("F", "M")
ADMIT_LEV    <- c("Observation", "Emergency", "Surgical_Elective")
NUTR_LEV     <- c("none", "enteral", "parenteral", "mixed")

# 需转成 0/1 数值的二分类（MICE 中为 factor）
BIN_COLS     <- c("sepsis3", "invasive_vent_at_t0", "vasoactive_at_t0",
                  "refeeding_risk")

# ---- 小工具 ----
to01 <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  as.integer(as.character(x))
}
mode1 <- function(x) {
  tb <- sort(table(x), decreasing = TRUE)
  as.integer(names(tb)[1])
}

# ---- 1. 读取并合并某个插补数据集的试验行 ----
# imp_id: 1..20；返回每个试验一行，含 origin_day / arm / DAG14
load_trials <- function(imp_id) {
  aug <- as.data.table(read_parquet(AUG_FILE))
  impf <- sprintf(MICE_FMT, as.integer(imp_id))
  imp <- as.data.table(read_parquet(impf))
  setkey(aug, stay_id, trial_num)
  setkey(imp, stay_id, trial_num)
  if (!identical(aug[, .(stay_id, trial_num)], imp[, .(stay_id, trial_num)]))
    stop("MICE 与 augmented 的 (stay_id,trial_num) 不对齐: imp", imp_id)

  dt <- aug[, .(stay_id, trial_num, time_zero, itt_group,
                death_28day_num, surv_time)]
  # 起点自然日（自然日历锚定）
  dt[, origin_day := as.integer(as.Date(time_zero) - min(as.Date(time_zero))),
     by = stay_id]
  # 臂：MICE treatment factor 0/1（与 itt_group 完全一致）
  dt[, arm := to01(imp$treatment)]

  # DAG14 协变量
  for (v in DAG_FROM_IMP) dt[, (v) := imp[[v]]]
  for (v in DAG_FROM_AUG) dt[, (v) := aug[[v]]]
  for (v in BIN_COLS)    dt[, (v) := to01(dt[[v]])]
  dt[, gender := factor(gender, levels = GENDER_LEV)]
  dt[, admission_type := factor(admission_type, levels = ADMIT_LEV)]
  dt[, nutrition_support_type := factor(nutrition_support_type, levels = NUTR_LEV)]
  dt
}

# ---- 2. stay 层共识死亡日 ----
# 返回 data.table: stay_id, death_day（幸存者无此行）
stay_death_day <- function(trials) {
  dth <- trials[death_28day_num == 1L,
                .(cand = pmin(TAU, pmax(0L,
                    origin_day + as.integer(round(surv_time))))),
                by = stay_id]
  dth[, .(death_day = mode1(cand)), by = stay_id]
}

# ---- 3. 构建日级 person-period 表 ----
# 返回 data.table，每个 stay 的 day = 0..max(origin)+28 连续
build_daily_input <- function(trials) {
  dd <- stay_death_day(trials)
  dt <- merge(trials, dd, by = "stay_id", all.x = TRUE)

  # 排除"起点严格晚于死亡日"的不可能试验（死亡后起点）
  dt <- dt[is.na(death_day) | origin_day <= death_day]

  # 每个 stay 的网格长度
  grid <- dt[, .(max_origin = max(origin_day),
                 dday = fifelse(any(death_28day_num == 1L),
                               death_day[1], NA_integer_)),
             by = stay_id]
  grid[, last_day := max_origin + TAU]

  # 展开网格
  daily <- grid[, .(day = seq(0L, last_day[1])), by = stay_id]

  # outcome：仅死亡日落 1
  daily <- merge(daily, grid[, .(stay_id, dday)], by = "stay_id", all.x = TRUE)
  daily[, outcome := 0L]
  daily[!is.na(dday) & day == dday, outcome := 1L]
  daily[, dday := NULL]

  # eligible 与 tx_init：起点日
  orig <- dt[, .(stay_id, origin_day, arm)]
  daily[, eligible := 0L]
  daily[, tx_init := 0L]
  daily[orig, on = c(stay_id = "stay_id", day = "origin_day"),
        `:=`(eligible = 1L, tx_init = i.arm)]

  # DAG14：写在起点日，随后前向填充（包按 *_bas 只读起点值）
  covcols <- DAG14
  for (v in covcols) {
    daily[, (v) := na_like(dt[[v]])]
    daily[dt, on = c(stay_id = "stay_id", day = "origin_day"),
          (v) := get(paste0("i.", v))]
    setorder(daily, stay_id, day)
    daily[, (v) := zoo_locf(get(v)), by = stay_id]
  }
  setorder(daily, stay_id, day)
  daily
}

# 生成与 x 同类型的标量 NA（用于初始化列）
na_like <- function(x) {
  if (is.factor(x))  return(factor(NA, levels = levels(x)))
  if (is.integer(x)) return(NA_integer_)
  if (is.numeric(x)) return(NA_real_)
  if (is.character(x)) return(NA_character_)
  NA
}

# data.table 无内置 LOCF；用 base 实现前向填充
zoo_locf <- function(x) {
  if (is.factor(x)) {
    lev <- levels(x)
    ii <- zoo_locf(as.integer(x))
    return(factor(ii, levels = seq_along(lev), labels = lev))
  }
  nas <- which(is.na(x))
  if (length(nas) == 0L) return(x)
  # day0 必有值；对每个 NA 用最近的非 NA 填充
  last <- NA
  out <- x
  for (i in seq_along(x)) {
    if (!is.na(x[i])) last <- x[i] else if (!is.na(last[1])) out[i] <- last
  }
  out
}

# ---- 4. 断言检查 ----
check_daily_input <- function(daily, label = "") {
  res <- list()
  # (a) 每个 stay 的 day 从 0 连续
  cc <- daily[, .(ok = identical(day, seq(0L, .N - 1L))), by = stay_id]
  res$contiguous <- all(cc$ok)
  # (b) (stay_id,day) 唯一
  res$unique_key <- nrow(daily[, .N, by = .(stay_id, day)][N > 1]) == 0L
  # (c) 无 NA
  res$no_na <- nrow(daily[!complete.cases(daily)]) == 0L
  # (d) eligible 计数 = 试验数；tx_init 两水平都在
  res$n_eligible <- daily[, sum(eligible)]
  res$tx_levels   <- sort(unique(daily$tx_init))
  # (e) 每个起点试验都有至少 28 天网格可供随访（包内 pmin 截到 28）
  last_dt <- daily[, .(last = max(day)), by = stay_id]
  span <- daily[eligible == 1L, .(stay_id, origin_day = day)]
  span <- merge(span, last_dt, by = "stay_id")
  span[, fu_avail := last - origin_day]
  res$all_tau <- all(span$fu_avail >= TAU)
  # (f) 事件只发生一次且在随访窗内
  ev <- daily[, .(nevt = sum(outcome)), by = stay_id]
  res$event_once <- all(ev$nevt <= 1L)

  cat(sprintf("[check %s] contiguous=%s unique=%s noNA=%s allTau=%s evtOnce=%s | n_elig=%d tx=%s\n",
      label, res$contiguous, res$unique_key, res$no_na, res$all_tau,
      res$event_once, res$n_eligible, paste(res$tx_levels, collapse = "/")))
  invisible(res)
}

# ---- 5. 从 SEQTaRget hazard 结果提取 logHR / SE ----
# h0 = c("Hazard ratio","LCI","UCI")（HR 尺度）；bootstrap.CI_method="se"
extract_log <- function(h0, ci_level = 0.95) {
  z <- qnorm(1 - (1 - ci_level) / 2)
  logq <- log(as.numeric(h0["Hazard ratio"]))
  se   <- (log(as.numeric(h0["UCI"])) - log(as.numeric(h0["LCI"]))) / (2 * z)
  list(loghr = logq, se = se, hr = as.numeric(h0["Hazard ratio"]),
       lci = as.numeric(h0["LCI"]), uci = as.numeric(h0["UCI"]))
}

# ---- 6. Rubin 合并（logHR 尺度）----
rubin_pool <- function(loghr, se, m = length(loghr)) {
  qbar <- mean(loghr)
  ubar <- mean(se^2)
  bv   <- if (m > 1) var(loghr) else 0
  tv   <- ubar + bv * (1 + 1 / m)
  dfm  <- if (m > 1 && bv > 0)
    (m - 1) * (1 + ubar / (bv * (1 + 1 / m)))^2 else Inf
  zt <- if (is.finite(dfm)) qt(0.975, dfm) else qnorm(0.975)
  list(loghr = qbar, hr = exp(qbar),
       lci = exp(qbar - zt * sqrt(tv)),
       uci = exp(qbar + zt * sqrt(tv)),
       p = 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
                else pnorm(-abs(qbar / sqrt(tv)))),
       T = tv, B = bv, df = dfm)
}
