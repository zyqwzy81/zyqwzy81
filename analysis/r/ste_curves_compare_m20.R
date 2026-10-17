#!/usr/bin/env Rscript
# =============================================================================
# v28.4_ste_curves_compare_m20.R — 生存曲线三法评估（m20，新死亡日几何）
#
# A法 STE g-comp：v28.4_ste_rmst_m20_curve.csv（pooled-logistic，arm×day，
#                 与 RMST 同模型，新 stay 共识死亡日，含 Rubin+bootstrap 带）
# B法 包内 km.curves：SEQTaRget SEQuential(km.curves=TRUE) 20 imp 点曲线，
#                 与 9 次要结局 oc_death 同模型同数据（包内 g-formula 标准化）
# C法 IPCW 加权 KM：fig2_fullkm_curve_m20b200.csv（注意：旧 ceiling 死亡日、
#                 死亡日修复前产物，仅作形态/量级参照，不做最终并列）
#
# 输出 tables/fig2/v28.4_curves_methods_m20.csv（long: method,arm,day,S）
#      tables/fig2/v28.4_curves_methods_summary_m20.csv（RMST/28d 风险对比）
# 用法: Rscript v28.4_ste_curves_compare_m20.R [imps,默认1..20]
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(SEQTaRget) })
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")

argv <- commandArgs(trailingOnly = TRUE)
IMPS <- if (length(argv) >= 1) as.integer(strsplit(argv[1], ",")[[1]]) else 1:20

OUT <- file.path(V28_BASE, "TTE/v28.4/output/tables/fig2")
B_F <- file.path(OUT, "v28.4_curves_pkg_m20_per_imp.csv")

cat("build event offsets once...\n")
off <- build_event_offsets()

run_imp <- function(imp_id) {
  cat(sprintf("\n================ pkg-km curves imp %d ================\n", imp_id))
  trials <- load_trials(imp_id)
  daily <- build_multi_daily(trials, off)
  opts <- SEQopts(km.curves = TRUE, bootstrap = FALSE, followup.max = TAU,
                  risk.times = 0:TAU, selection.first_trial = FALSE)
  mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                    eligible.col = "eligible", treatment.col = "tx_init",
                    outcome.col = "oc_death", time_varying.cols = DAG14,
                    fixed.cols = list(), method = "ITT",
                    options = opts, verbose = FALSE)
  rd <- as.data.table(mod@risk.data[[1]])
  cat("  risk.data cols:", paste(names(rd), collapse = ","), "\n")
  rd <- rd[, .(day = as.integer(Followup), arm = A,
               risk = as.numeric(Risk))]
  rd[, S := 1 - risk]
  rd[, imp := imp_id]
  rd[]
}

# 断点续跑（per-imp 行独立累积；缺任一 imp 则整体重跑所请求的 IMPS）
done <- if (file.exists(B_F)) unique(fread(B_F)$imp) else integer(0)
todo <- setdiff(IMPS, done)
if (length(todo) && file.exists(B_F)) { file.remove(B_F); todo <- IMPS }
for (im in todo) {
  x <- tryCatch(run_imp(im), error = function(e) {
    cat("  ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(x)) fwrite(x, B_F, append = file.exists(B_F))
}

bp <- fread(B_F); bp <- bp[imp %in% IMPS]
cat("Followup grid:", paste(sort(unique(bp$day)), collapse = ","), "\n")
# 与 STE 口径对齐：仅保留 day 0..TAU
bp <- bp[day >= 0 & day <= TAU]
bmean <- bp[, .(S = mean(S)), by = .(arm, day)]
bmean[, `:=`(method = "B_pkg_gformula", arm = as.character(arm))]
setcolorder(bmean, c("method", "arm", "day", "S"))

# A 法
a <- fread(file.path(OUT, "v28.4_ste_rmst_m20_curve.csv"))
a0 <- a[, .(method = "A_ste_gcomp", arm = "0", day, S = S0)]
a1 <- a[, .(method = "A_ste_gcomp", arm = "1", day, S = S1)]
ac <- rbind(a0, a1)

# C 法（修复前，标记）
cfile <- file.path(OUT, "fig2_fullkm_curve_m20b200.csv")
cc <- if (file.exists(cfile)) {
  x <- fread(cfile)
  x[, .(method = "C_ipcw_km_predathfix", arm = as.character(arm),
        day = as.integer(day), S = as.numeric(S))]
} else NULL

allc <- rbindlist(list(ac, bmean, cc), fill = TRUE)
fwrite(allc, file.path(OUT, "v28.4_curves_methods_m20.csv"))

# RMST28（积分 day 0..27）与 28d 风险（day 28）
sm <- allc[, .(
  rmst = sum(S[day < TAU]),
  risk28 = 1 - S[day == TAU][1]
), by = .(method, arm)]
fwrite(sm, file.path(OUT, "v28.4_curves_methods_summary_m20.csv"))
cat("\n===== methods summary (m=20) =====\n")
print(sm[order(arm, method)])
cat("\n臂间差（治疗-对照）:\n")
print(dcast(sm[, .(method, arm, rmst, risk28)],
            method ~ arm, value.var = c("rmst", "risk28"))[
            , .(method, rmst_diff = rmst_1 - rmst_0,
                risk_diff = risk28_1 - risk28_0)])
cat("\n保存:", file.path(OUT, "v28.4_curves_methods_m20.csv"), "\n")
