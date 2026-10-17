#!/usr/bin/env Rscript
# =============================================================================
# v28.4_ste_rmst_m20.R
# STE g-computation RMST — m=20 MI / B=200 subject-cluster bootstrap / Rubin
#
# v28.4.1（2026-10-14）两处方法学修正：
#   1) 死亡日改用 stay 层绝对日历日共识（v28.4_common.R::stay_death_day），
#      不再用 ceiling(death_time_from_t0) 逐起点规则；与主分析完全同几何。
#   2) 离散时间 pooled logistic 风险集正确截断：每试验仅贡献 day<=dday
#      （事件日 Y=1，事件前 Y=0），不再保留死后 person-day（旧脚本死后
#      行仍以 Y=0 入模，系统性压低 hazard）。
# 协变量 DAG14（与主分析同一集合、起点基线值）；arm×day 交互（Limozin 2025）。
# 每 imp：点估计 + B=200 subject 级 cluster bootstrap（seed=SEED+imp）。
# 断点续跑：per-imp 结果逐 imp 落盘，重跑自动跳过已完成 imp。
# 用法: OMP_NUM_THREADS=2 LD_LIBRARY_PATH=/opt/miniforge3/lib \
#         Rscript v28.4_ste_rmst_m20.R [imps,默认1..20] [B,默认200]
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(arrow) })
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_common.R")

argv <- commandArgs(trailingOnly = TRUE)
IMPS <- if (length(argv) >= 1) as.integer(strsplit(argv[1], ",")[[1]]) else 1:20
B    <- if (length(argv) >= 2) as.integer(argv[2]) else 200L

OUT <- file.path(V28_BASE, "TTE/v28.4/output/tables/fig2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PER_F <- file.path(OUT, "v28.4_ste_rmst_m20_per_imp.csv")
BOOT_F <- file.path(OUT, "v28.4_ste_rmst_m20_boot_daily.csv")
CURVE_F <- file.path(OUT, "v28.4_ste_rmst_m20_curve_per_imp.csv")

# subject_id 映射（load_trials 不含）
aug <- as.data.table(read_parquet(AUG_FILE))
subj_map <- unique(aug[, .(stay_id, trial_num, subject_id)])
rm(aug); gc()

# 因子水平（与 common 一致）
LEV <- list(gender = GENDER_LEV, admission_type = ADMIT_LEV,
            nutrition_support_type = NUTR_LEV)

# ---- 逐 imp 构建 person 表（新死亡日几何）----
build_person <- function(imp_id) {
  trials <- load_trials(imp_id)
  dd <- stay_death_day(trials)
  dt <- merge(trials, dd, by = "stay_id", all.x = TRUE)
  dt <- dt[is.na(death_day) | origin_day <= death_day]
  dt[, rel := death_day - origin_day]
  # 事件归因以 stay 层共识绝对死亡日为准（与 build_multi_daily/主分析一致）：
  # 窗口 [origin, origin+28] 含日 28 边界（32 例）；不与行级 death_28day_num
  # 相交（38 个晚起点行级标签为 0 但共识日落窗内，主管线口径为事件）。
  dt[, event := as.integer(!is.na(rel) & rel >= 0L & rel <= TAU)]
  # 事件日索引 0..TAU（29 个离散区间，day28 边界事件入模），否则按 TAU 截尾
  dt[, dday := fifelse(event == 1L, as.integer(rel), TAU)]
  dt <- merge(dt, subj_map, by = c("stay_id", "trial_num"), all.x = TRUE)
  dt[, idx := seq_len(.N)]
  dt
}

# 离散时间风险集 long：day 0..TAU（29 区间），仅 day <= dday；事件日 Y=1。
# 注意 dday==TAU 既可能是 day28 边界事件（event=1）也可能是截尾，需 event 区分。
build_long <- function(pdat) {
  lng <- merge(CJ(idx = pdat$idx, day = 0:TAU),
               pdat[, c("idx", "stay_id", "subject_id", "arm", "dday",
                        "event", DAG14), with = FALSE],
               by = "idx")
  lng <- lng[day <= dday]
  lng[, Y := as.integer(day == dday & event == 1L)]
  lng
}

FML <- as.formula(paste("Y ~ arm + day + arm:day +",
                        paste(DAG14, collapse = " + ")))

gcomp <- function(lng, pdat, wcol = NULL) {
  if (is.null(wcol)) {
    lng[, .w := 1]; pdat[, .w := 1]
  } else {
    lng[, .w := get(wcol)]; pdat[, .w := get(wcol)]
  }
  fit <- suppressWarnings(glm(FML, data = lng, family = binomial(),
                              weights = .w, control = list(maxit = 200)))
  pg <- CJ(idx = pdat$idx, day = 0:TAU, arm = c(0L, 1L))
  pg <- merge(pg, pdat[, c("idx", ".w", DAG14), with = FALSE], by = "idx")
  pg[, P := predict(fit, newdata = pg, type = "response")]
  pg <- pg[order(idx, arm, day)]
  pg[, S := cumprod(1 - P), by = .(idx, arm)]
  ES <- pg[, .(S = weighted.mean(S, .w)), by = .(day, arm)]
  ESw <- dcast(ES, day ~ arm, value.var = "S")
  setnames(ESw, c("day", "S0", "S1"))
  # RMST28 = 积分 0..27（28 个整日）；日 28 行仅用于 28 天风险与边界事件
  ESm <- ESw[day < TAU]
  list(rmst0 = sum(ESm$S0), rmst1 = sum(ESm$S1),
       diff = sum(ESm$S1) - sum(ESm$S0),
       risk0 = 1 - ESw[day == TAU]$S0, risk1 = 1 - ESw[day == TAU]$S1,
       curve = ESw)
}

# 一次 subject 级 cluster 重抽：以 multiplicity 权重实现（等价行复制）
cluster_sample <- function(lng, pdat, subs) {
  pick <- data.table(subject_id = sample(subs, length(subs), replace = TRUE))
  wtab <- pick[, .(w = .N), by = subject_id]
  bl <- merge(lng, wtab, by = "subject_id")
  bp <- merge(pdat, wtab, by = "subject_id")
  list(bl = bl, bp = bp)
}

run_imp <- function(imp_id) {
  cat(sprintf("\n================ STE imp %d ================\n", imp_id))
  pdat <- build_person(imp_id)
  lng <- build_long(pdat)
  cat(sprintf("  TP=%d stays=%d subjects=%d events=%d person-days=%s\n",
              nrow(pdat), uniqueN(pdat$stay_id), uniqueN(pdat$subject_id),
              sum(pdat$event), formatC(nrow(lng), big.mark = ",")))
  g <- gcomp(lng, pdat)
  cat(sprintf("  point: RMST0=%.4f RMST1=%.4f diff=%+.4f | risk0=%.4f risk1=%.4f\n",
              g$rmst0, g$rmst1, g$diff, g$risk0, g$risk1))

  # subject 级 cluster bootstrap
  set.seed(SEED + imp_id)
  subs <- unique(lng$subject_id)
  bd <- matrix(NA_real_, B, 3); colnames(bd) <- c("rmst0", "rmst1", "diff")
  bc <- vector("list", B)
  nb <- 0L; t0 <- Sys.time()
  while (nb < B) {
    cs <- cluster_sample(lng, pdat, subs)
    gg <- tryCatch(gcomp(cs$bl, cs$bp, wcol = "w"),
                   error = function(e) NULL)
    if (is.null(gg) || !all(is.finite(c(gg$rmst0, gg$rmst1)))) next
    nb <- nb + 1L
    bd[nb, ] <- c(gg$rmst0, gg$rmst1, gg$diff)
    bc[[nb]] <- gg$curve
    if (nb %% 50 == 0)
      cat(sprintf("  boot %d/%d  %.0fs\n", nb, B,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  curve_boot <- rbindlist(lapply(seq_len(B), function(i)
    cbind(boot = i, bc[[i]])), use.names = TRUE)
  curve_boot[, imp := imp_id]
  fwrite(curve_boot, BOOT_F, append = file.exists(BOOT_F))
  list(per = data.table(imp = imp_id, rmst0 = g$rmst0, rmst1 = g$rmst1,
                        diff = g$diff, risk0 = g$risk0, risk1 = g$risk1,
                        se = sd(bd[, "diff"], na.rm = TRUE),
                        n_tp = nrow(pdat), n_event = sum(pdat$event),
                        n_boot = nb, curve = list(g$curve)),
       curve_boot = curve_boot)
}

# ---- 断点续跑 ----
done_imps <- if (file.exists(PER_F)) fread(PER_F)$imp else integer(0)
todo <- setdiff(IMPS, done_imps)
# 若有将被重跑的 imp（PER/BOOT/CURVE 三文件须同步），整体重置
if (length(todo)) {
  for (f in c(PER_F, BOOT_F, CURVE_F))
    if (file.exists(f)) file.remove(f)
  done_imps <- integer(0)
  todo <- IMPS
}
cat(sprintf("STE m20: todo=%s done=%s B=%d\n",
            paste(IMPS, collapse = ","),
            if (length(done_imps)) paste(done_imps, collapse = ",") else "none", B))

res_l <- list()
for (im in todo) {
  r <- run_imp(im)
  per_row <- copy(r$per); per_row$curve <- NULL
  fwrite(per_row, PER_F, append = file.exists(PER_F))
  cv <- copy(r$per$curve[[1]]); cv[, imp := im]
  fwrite(cv, CURVE_F, append = file.exists(CURVE_F))
  res_l[[length(res_l) + 1]] <- r
}

# ---- Rubin 合并（读全量 per 文件）----
per <- fread(PER_F); per <- per[imp %in% IMPS][order(imp)]
m <- nrow(per)
qbar <- mean(per$diff)
bv <- if (m > 1) var(per$diff) else 0
ubar <- mean(per$se^2, na.rm = TRUE)
tv <- ubar + (1 + 1 / m) * bv
dfm <- if (bv > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bv))^2 else Inf
zt <- if (is.finite(dfm)) qt(0.975, dfm) else qnorm(0.975)
pval <- 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
             else pnorm(-abs(qbar / sqrt(tv))))
pooled <- data.table(analysis = "STE_gcomp_RMST_m20", tau = TAU, m = m,
                     rmst0 = mean(per$rmst0), rmst1 = mean(per$rmst1),
                     diff = qbar, ci_lower = qbar - zt * sqrt(tv),
                     ci_upper = qbar + zt * sqrt(tv), p_value = pval,
                     risk0 = mean(per$risk0), risk1 = mean(per$risk1),
                     between_var = bv, within_var = ubar, total_var = tv,
                     n_boot = B, rubin_df = dfm)
fwrite(pooled, file.path(OUT, "v28.4_ste_rmst_m20_pooled.csv"))

# ---- 合并生存曲线：20 imp 点曲线均值 + 日级 Rubin 方差带 ----
cv_imp <- fread(CURVE_F); cv_imp <- cv_imp[imp %in% IMPS]
bt <- fread(BOOT_F);  bt <- bt[imp %in% IMPS]
# within：每日 bootstrap 方差按 imp 平均；between：20 imp 点曲线日间方差
cv_mean <- cv_imp[, .(S0 = mean(S0), S1 = mean(S1)), by = day]
w0 <- bt[, .(w = mean(var(S0))), by = .(imp, day)][, .(u0 = mean(w)), by = day]
w1 <- bt[, .(w = mean(var(S1))), by = .(imp, day)][, .(u1 = mean(w)), by = day]
b0 <- cv_imp[, .(b0 = var(S0)), by = day]
b1 <- cv_imp[, .(b1 = var(S1)), by = day]
cv_final <- Reduce(function(x, y) merge(x, y, by = "day"),
                   list(cv_mean, w0, w1, b0, b1))
cv_final[, `:=`(
  S0_lci = S0 - 1.96 * sqrt(u0 + (1 + 1 / m) * b0),
  S0_uci = S0 + 1.96 * sqrt(u0 + (1 + 1 / m) * b0),
  S1_lci = S1 - 1.96 * sqrt(u1 + (1 + 1 / m) * b1),
  S1_uci = S1 + 1.96 * sqrt(u1 + (1 + 1 / m) * b1))]
fwrite(cv_final, file.path(OUT, "v28.4_ste_rmst_m20_curve.csv"))

cat("\n===== RUBIN POOLED STE RMST (m=", m, ", B=", B, ") =====\n", sep = "")
print(pooled)
cat("\n输出:\n ", file.path(OUT, "v28.4_ste_rmst_m20_pooled.csv"),
    "\n ", file.path(OUT, "v28.4_ste_rmst_m20_curve.csv"), "\n")
