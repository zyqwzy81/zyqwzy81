# ============================================================================
# multi_common.R
# 多次入组（multi-enrollment）序贯目标试验公共库
#
# 思路：在 common.R 的"stay 单条连续自然日时间线"几何上，
#       把每个二分类结局的【绝对事件日】铺到同一条时间线（每结局一列），
#       包在 SEQuential() 内做 expand + 窗口化，自动把事件归因到
#       各试验起点的克隆。因此同一份宽日级表可驱动全部 28 天结局，
#       仅需切换 outcome.col。
#
# 事件日来源（全部逐起点 32,330 行，时间锚定 time_zero）：
#   * 28d 全因死亡  oc_death：stay 共识死亡日（沿用主分析）
#   * ICU 死亡      oc_icu  / 院内死亡 oc_hosp：死亡日（stay 层死亡场所一致）
#   * 新发心律失常  oc_arrh：M0_14a rhythm_event_time
#   * 心肌损伤      oc_myo ：M0_14c first_troponin_elevated_time
#   * 新发 AKI≥2    oc_aki ：kdigo_stages 起点后首次 stage>=2（排除起点已≥2）
#   * PRBC 输血     oc_prbc：inputevents 起点后首次 PRBC（225168/226368）
#   * CV 2/3 复合   oc_cv2/oc_cv3：组分事件日的最早者
# 另：血管加压升级是 0-24h 标志（vaso24），以 followup.max=1 单独估计。
#
# 单次入组敏感性：SEQuential 内 selection.first_trial=TRUE（包内置）。
# 版本 V28.4  日期 2026-10-10
# ============================================================================

CODE_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r"
source(file.path(CODE_DIR, "common.R"))
suppressPackageStartupMessages(library(SEQTaRget))

V28_BASE <- "/data/小米云盘/文档/论文/血磷"
RRES     <- file.path(V28_BASE, "sql/result")
EXTRACT  <- file.path(V28_BASE, "TTE/v28.4/output/extraction")
SEED     <- 20260629L

# 28 天二分类结局列（宽表列名）
BIN_OUTCOMES <- c("oc_death", "oc_icu", "oc_hosp",
                  "oc_arrh", "oc_myo", "oc_aki", "oc_prbc",
                  "oc_cv2", "oc_cv3")

# ---- 相对事件日（floor，天）；NA 表示窗内无事件 ----
rel_day <- function(event_t, tz) {
  d <- as.integer(floor(as.numeric(difftime(event_t, tz, units = "days"))))
  ifelse(is.na(d) | d < 0L | d > TAU, NA_integer_, d)
}

# ---- 1. 逐起点事件偏移表（不随 imp 变，只算一次）----
# 返回 data.table: stay_id, time_zero, origin_anchor 不需要；
#   e_arrh,e_myo,e_aki,e_prbc（相对日 0..28 或 NA），
#   vaso24（0/1），death_icu/death_hospital（stay 层场所）
build_event_offsets <- function() {
  # 心律失常
  a <- as.data.table(read_parquet(file.path(RRES, "M0_14a_arrhythmia.parquet")))
  a <- a[, .(stay_id, time_zero, new_arrhythmia_rhythm, rhythm_event_time)]
  a[, e_arrh := fifelse(new_arrhythmia_rhythm == 1L,
                        rel_day(rhythm_event_time, time_zero), NA_integer_)]

  # 心肌损伤
  m <- as.data.table(read_parquet(file.path(RRES, "M0_14c_myocardial.parquet")))
  m <- m[, .(stay_id, time_zero, myocardial_injury, first_troponin_elevated_time)]
  m[, e_myo := fifelse(myocardial_injury == 1L,
                       rel_day(first_troponin_elevated_time, time_zero), NA_integer_)]

  # KDIGO 时间序列（含各分期时点）
  ks <- as.data.table(read_parquet(file.path(EXTRACT, "kdigo_stages_ts.parquet")))
  # AKI 定义：肌酐分期或新发 RRT（creat-or-RRT），避免 UO 分量过度触发；
  # 与项目 M0_13 肌酐口径及多数目标试验文献一致
  ks[, aki_cr := as.integer(pmax(aki_stage_creat, aki_stage_crrt, na.rm = TRUE))]
  ks[is.infinite(aki_cr) | is.na(aki_cr), aki_cr := -1L]
  setorder(ks, stay_id, charttime)
  # 仅保留 aki_cr>=2 的记录，供起点后非等价连接
  ks2 <- ks[aki_cr >= 2L, .(stay_id, ctime = charttime)]

  # PRBC 时间序列
  pb <- as.data.table(read_parquet(file.path(EXTRACT, "prbc_events.parquet")))
  setorder(pb, stay_id, starttime)

  # 起点前最近一次 AKI（creat-or-RRT）分期（rolling：<=time_zero 最后一条），识别 prevalent
  kpr <- ks[, .(stay_id, charttime, stage_pre = aki_cr)]
  anc <- a[, .(stay_id, time_zero)]
  setnames(anc, "time_zero", "charttime")
  setorder(kpr, stay_id, charttime); setorder(anc, stay_id, charttime)
  rolled <- kpr[anc, on = c("stay_id", "charttime"), roll = Inf, mult = "last"]
  setnames(rolled, "charttime", "time_zero")
  anchor <- rolled[, .(stay_id, time_zero, stage_pre)]
  setorder(anchor, stay_id, time_zero)

  # AKI：对每起点，取严格在 tz 之后、28 天内首次 aki_cr>=2
  cat("[offsets] computing post-origin AKI (creat/RRT) events...\n")
  # 非等价连接：同 stay，ctime>time_zero
  ak <- ks2[anchor, on = "stay_id", allow.cartesian = TRUE]
  ak <- ak[ctime > time_zero]
  ak[, rd := as.integer(floor(as.numeric(difftime(ctime, time_zero, units = "days"))))]
  ak <- ak[rd >= 0L & rd <= TAU]
  aki_first <- ak[order(stay_id, time_zero, rd),
                  .SD[1], by = .(stay_id, time_zero)]
  aki_first <- aki_first[, .(stay_id, time_zero, rd_aki = rd)]

  # PRBC：起点后首次
  cat("[offsets] computing post-origin PRBC events...\n")
  pj <- pb[anchor, on = "stay_id", allow.cartesian = TRUE]
  pj <- pj[starttime > time_zero]
  pj[, rd := as.integer(floor(as.numeric(difftime(starttime, time_zero, units = "days"))))]
  pj <- pj[rd >= 0L & rd <= TAU]
  pb_first <- pj[order(stay_id, time_zero, rd),
                 .SD[1], by = .(stay_id, time_zero)]
  pb_first <- pb_first[, .(stay_id, time_zero, rd_prbc = rd)]

  # 血管加压升级（0-24h 标志，已逐起点）
  h <- as.data.table(read_parquet(file.path(RRES, "M0_14b_hemodynamic.parquet")))
  h <- h[, .(stay_id, time_zero,
             vaso24 = as.integer(vasopressor_escalation))]

  # 死亡场所（stay 层，核对一致性）
  o07 <- as.data.table(read_parquet(file.path(RRES, "M0_07_outcome_extraction.parquet")))
  setorder(o07, stay_id, time_zero)
  dplace <- o07[, .(icu_any = as.integer(any(death_icu == 1L)),
                    hosp_any = as.integer(any(death_hospital == 1L))), by = stay_id]

  # 合并逐起点表
  off <- anchor[, .(stay_id, time_zero, stage_pre)]
  off <- merge(off, a[, .(stay_id, time_zero, e_arrh)], by = c("stay_id", "time_zero"), all.x = TRUE)
  off <- merge(off, m[, .(stay_id, time_zero, e_myo)], by = c("stay_id", "time_zero"), all.x = TRUE)
  off <- merge(off, aki_first, by = c("stay_id", "time_zero"), all.x = TRUE)
  off <- merge(off, pb_first, by = c("stay_id", "time_zero"), all.x = TRUE)
  off <- merge(off, h, by = c("stay_id", "time_zero"), all.x = TRUE)

  # AKI 仅保留"新发"：起点前最近分期 >=2 者排除（persistent/prevalent）
  off[, e_aki := fifelse(!is.na(stage_pre) & stage_pre >= 2L, NA_integer_, rd_aki)]
  off[, e_prbc := rd_prbc]

  # 复合：vaso 视为第 1 天；取组分最早日
  off[, vaso_d := fifelse(vaso24 == 1L, 1L, NA_integer_)]
  off[, e_cv2 := pmin(e_arrh, vaso_d, na.rm = TRUE)]
  off[, e_cv3 := pmin(e_arrh, e_myo, vaso_d, na.rm = TRUE)]
  off[is.infinite(e_cv2), e_cv2 := NA_integer_]
  off[is.infinite(e_cv3), e_cv3 := NA_integer_]

  attr(off, "dplace") <- dplace
  cat(sprintf("[offsets] rows=%d stays=%d\n", nrow(off), uniqueN(off$stay_id)))
  off
}

# ---- 2. 多次入组宽日级输入（含全部 28 天结局列）----
# trials: load_trials(imp_id)；off: build_event_offsets()
build_multi_daily <- function(trials, off, H = TAU) {
  dd <- stay_death_day(trials)
  dt <- merge(trials, dd, by = "stay_id", all.x = TRUE)
  dt <- dt[is.na(death_day) | origin_day <= death_day]

  grid <- dt[, .(max_origin = max(origin_day),
                 dday = fifelse(any(death_28day_num == 1L), death_day[1], NA_integer_)),
             by = stay_id]
  grid[, last_day := max_origin + H]

  daily <- grid[, .(day = seq(0L, last_day[1])), by = stay_id]

  # eligible / tx_init
  orig <- dt[, .(stay_id, origin_day, arm)]
  daily[, `:=`(eligible = 0L, tx_init = 0L)]
  daily[orig, on = c(stay_id = "stay_id", day = "origin_day"),
        `:=`(eligible = 1L, tx_init = i.arm)]

  # DAG14：起点日写入 + 前向填充
  for (v in DAG14) {
    daily[, (v) := na_like(dt[[v]])]
    daily[dt, on = c(stay_id = "stay_id", day = "origin_day"),
          (v) := get(paste0("i.", v))]
    setorder(daily, stay_id, day)
    daily[, (v) := zoo_locf(get(v)), by = stay_id]
  }

  # ---- 结局列初始化 ----
  for (oc in BIN_OUTCOMES) daily[, (oc) := 0L]

  # oc_death：死亡日
  daily[grid[!is.na(dday)], on = c(stay_id = "stay_id", day = "dday"),
        oc_death := 1L]

  # oc_icu / oc_hosp：死亡日 × 死亡场所
  dplace <- as.data.table(attr(off, "dplace"))
  g_dplace <- merge(grid[!is.na(dday), .(stay_id, dday)], dplace, by = "stay_id")
  daily[g_dplace[icu_any == 1L], on = c(stay_id = "stay_id", day = "dday"),
        oc_icu := 1L]
  daily[g_dplace[hosp_any == 1L], on = c(stay_id = "stay_id", day = "dday"),
        oc_hosp := 1L]

  # ---- 其余结局：逐起点绝对事件日 → 去重铺到时间线 ----
  # off 与 dt 均逐 (stay_id,time_zero)，按该键对齐取 origin_day
  omap <- merge(dt[, .(stay_id, time_zero, origin_day)],
                off, by = c("stay_id", "time_zero"), all.x = TRUE)

  place_events <- function(col_offset, col_out) {
    ev <- omap[!is.na(get(col_offset)),
               .(abs_day = unique(origin_day + get(col_offset))), by = stay_id]
    daily[ev, on = c(stay_id = "stay_id", day = "abs_day"),
          (col_out) := 1L]
  }
  place_events("e_arrh", "oc_arrh")
  place_events("e_myo",  "oc_myo")
  place_events("e_aki",  "oc_aki")
  place_events("e_prbc", "oc_prbc")
  place_events("e_cv2",  "oc_cv2")
  place_events("e_cv3",  "oc_cv3")

  setorder(daily, stay_id, day)
  daily
}

# ---- 3. 宽表断言 ----
check_multi_daily <- function(daily, label = "") {
  cc <- daily[, .(ok = identical(day, seq(0L, .N - 1L))), by = stay_id]
  n_elig <- daily[, sum(eligible)]
  counts <- daily[, lapply(.SD, sum), .SDcols = BIN_OUTCOMES]
  cat(sprintf("[multi-check %s] contiguous=%s n_elig=%d\n",
              label, all(cc$ok), n_elig))
  print(counts)
  invisible(list(contiguous = all(cc$ok), n_elig = n_elig, counts = counts))
}
