# ============================================================================
# v28.4_subgroup_multi.R
# 多次入组（multi-enrollment）框架下：9 个基线变量亚组 HR + 层间异质性检验
#
# 亚组通过包内置 SEQopts(subgroup=...) 在同一次 g-formula 展开中分层；
# 血磷阈值 <0.50 作为 baseline phosphate 亚组（不再单列阈值分析）。
# 异质性检验：Δ=logHR_g1-logHR_g0，SE=sqrt(SE_g1^2+SE_g0^2)（层为
# 互斥患者群），z 检验。
#
# 开关同主脚本：--imps/--nboot/--bootimps/--firstonly/--tag
# 版本 V28.4  日期 2026-10-10
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")

argv <- commandArgs(trailingOnly = TRUE)
gv <- function(k, d) { i <- match(k, argv); if (is.na(i)) d else argv[i + 1] }
IMPS   <- as.integer(strsplit(gv("--imps", "1"), ",")[[1]])
NBOOT  <- as.integer(gv("--nboot", "200"))
DEF_BOOT <- if (NBOOT > 0) paste(IMPS, collapse = ",") else ""
BOOTIMPS <- as.integer(strsplit(gv("--bootimps", DEF_BOOT), ",")[[1]])
TAG    <- gv("--tag", "sub")
FIRST  <- "--firstonly" %in% argv

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/subgroup_multi")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

suffix <- paste0(if (FIRST) "firstonly_" else "", TAG)
PER_F <- file.path(OUT_T, paste0("v28.4_subgroup_", suffix, "_per_imp.csv"))

cat("build event offsets once...\n")
off <- build_event_offsets()

# ---- 亚组拆分定义 ----
# col=新增分层列；labs[1]=flag0 层标签,labs[2]=flag1 层标签
# flag 作用于起点日（eligible==1）已写入 daily 的基线变量
SPEC <- list(
  list(variable = "baseline_phosphate", col = "sub_phos",
       labs = c("Phosphate >= 0.50 mmol/L", "Phosphate < 0.50 mmol/L"),
       flag = function(o) as.integer(o$baseline_phosphate < 0.50)),
  list(variable = "sepsis3", col = "sub_sepsis",
       labs = c("Non-Sepsis-3", "Sepsis-3"),
       flag = function(o) as.integer(o$sepsis3 == 1L)),
  list(variable = "baseline_sofa", col = "sub_sofa",
       labs = c("SOFA < 8", "SOFA >= 8"),
       flag = function(o) as.integer(o$baseline_sofa >= 8)),
  list(variable = "age", col = "sub_age",
       labs = c("Age < 65 yr", "Age >= 65 yr"),
       flag = function(o) as.integer(o$age >= 65)),
  list(variable = "gender", col = "sub_sex",
       labs = c("Female", "Male"),
       flag = function(o) as.integer(o$gender == "M")),
  list(variable = "invasive_vent_at_t0", col = "sub_vent",
       labs = c("No invasive ventilation", "Invasive ventilation"),
       flag = function(o) as.integer(o$invasive_vent_at_t0 == 1L)),
  list(variable = "vasoactive_at_t0", col = "sub_vaso",
       labs = c("No vasoactive agents", "Vasoactive agents"),
       flag = function(o) as.integer(o$vasoactive_at_t0 == 1L)),
  list(variable = "egfr", col = "sub_egfr",
       labs = c("eGFR >= 60", "eGFR < 60"),
       flag = function(o) as.integer(o$egfr < 60)),
  list(variable = "refeeding_risk", col = "sub_refeed",
       labs = c("No refeeding risk", "Refeeding risk"),
       flag = function(o) as.integer(o$refeeding_risk == 1L))
)

# 在 daily 上为某个 spec 写入分层列（起点日写 flag，前向填充）
add_subcol <- function(daily, trials, sp) {
  col <- sp$col
  # gender 不在 DAG14，需先补到 daily（stay 层常量）
  if (!"gender" %in% names(daily)) {
    gmap <- unique(trials[, .(stay_id, gender)])
    daily[gmap, on = "stay_id", gender := i.gender]
  }
  daily[, (col) := NA_integer_]
  o <- daily[eligible == 1L]
  o <- o[, .(stay_id, day, fl = sp$flag(o))]
  daily[o, on = c("stay_id", "day"), (col) := i.fl]
  setorder(daily, stay_id, day)
  daily[, (col) := zoo_locf(get(col)), by = stay_id]
  # 前向填充后起点前仍可能 NA，用首个值回补（分层仅在展开的试验行使用）
  daily[, (col) := { x <- get(col); if (all(is.na(x))) x else
    zoo::na.locf(x, fromLast = FALSE, na.rm = FALSE) }, by = stay_id]
  daily[, (col) := zoo::na.locf(get(col), na.rm = FALSE, fromLast = TRUE),
        by = stay_id]
  invisible(daily)
}

run_imp <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  trials <- load_trials(imp_id)
  daily0 <- build_multi_daily(trials, off)
  do_boot <- NBOOT > 0 && imp_id %in% BOOTIMPS
  res <- list()
  for (sp in SPEC) {
    daily <- copy(daily0)
    add_subcol(daily, trials, sp)
    opts <- SEQopts(hazard = TRUE, subgroup = sp$col,
                    bootstrap = do_boot,
                    bootstrap.nboot = max(1L, NBOOT),
                    bootstrap.sample = 1,
                    followup.max = TAU,
                    seed = SEED + if (do_boot) imp_id else 0L,
                    bootstrap.CI_method = "se",
                    selection.first_trial = FIRST)
    mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                      eligible.col = "eligible", treatment.col = "tx_init",
                      outcome.col = "oc_death", time_varying.cols = DAG14,
                      fixed.cols = list(sp$col), method = "ITT",
                      options = opts, verbose = FALSE)
    hh <- hazard_ratio(mod)           # list: "<col>_0"->vec, "<col>_1"->vec
    one <- function(lvl) {
      h <- hh[[paste0(sp$col, "_", lvl)]]
      lhr <- unname(log(h["Hazard ratio"]))
      se <- NA_real_
      if (do_boot && is.finite(h["LCI"]) && is.finite(h["UCI"]))
        se <- (as.numeric(log(h["UCI"])) - as.numeric(log(h["LCI"]))) /
              (2 * qnorm(0.975))
      c(loghr = lhr, se = se)
    }
    res[[sp$variable]] <- list(g0 = one("0"), g1 = one("1"))
  }
  res
}

# ---- checkpoint ----
done_imps <- if (file.exists(PER_F)) {
  x <- fread(PER_F, colClasses = list(imp = "integer"))$imp
  sort(unique(x))
} else integer(0)
todo <- setdiff(IMPS, done_imps)
if (length(done_imps) > 0 && !setequal(done_imps, intersect(done_imps, IMPS))) {
  cat(sprintf("[checkpoint] PER_F imps=%s ⊄ IMPS → reset\n",
              paste(done_imps, collapse = ",")))
  file.remove(PER_F); done_imps <- integer(0); todo <- IMPS
}
cat(sprintf("[checkpoint] todo=%s done=%s\n",
            if (length(todo)) paste(todo, collapse = ",") else "(none)",
            if (length(done_imps)) paste(done_imps, collapse = ",") else "(none)"))

# ---- 跨 imp 执行（长表 per_imp：imp × variable × level）----
res_l <- vector("list", length(todo))
for (k in seq_along(todo)) {
  im <- todo[k]
  r <- run_imp(im)
  per_rows <- rbindlist(lapply(SPEC, function(sp) {
    g0 <- r[[sp$variable]]$g0; g1 <- r[[sp$variable]]$g1
    rbind(
      data.table(imp = im, variable = sp$variable, level = "g0",
                 loghr = g0["loghr"], se = g0["se"]),
      data.table(imp = im, variable = sp$variable, level = "g1",
                 loghr = g1["loghr"], se = g1["se"])
    )
  }))
  fwrite(per_rows, PER_F, append = file.exists(PER_F))
  res_l[[k]] <- r
  rm(r, per_rows); invisible(gc())
}

# ---- Rubin 合并（从 PER_F 读全量）----
per_dt <- fread(PER_F, colClasses = list(imp = "integer"))
stopifnot(setequal(sort(unique(per_dt$imp)), IMPS))

pool_stratum <- function(v, lv) {
  sub <- per_dt[variable == v & level == lv]
  q <- sub$loghr; se_all <- sub$se
  qbar <- mean(q); bv <- if (length(q) > 1) var(q) else 0
  sev <- se_all[is.finite(se_all)]
  if (length(sev)) {
    ubar <- mean(sev^2); tv <- ubar + (1 + 1 / length(q)) * bv
    dfm <- if (bv > 0) (length(q) - 1) * (1 + ubar / ((1 + 1 / length(q)) * bv))^2 else Inf
  } else { tv <- NA_real_; dfm <- Inf }
  zt <- if (is.finite(dfm)) qt(.975, dfm) else qnorm(.975)
  list(loghr = qbar, se = if (is.finite(tv)) sqrt(tv) else NA_real_,
       df = dfm,
       lci = if (is.finite(tv)) qbar - zt * sqrt(tv) else NA_real_,
       uci = if (is.finite(tv)) qbar + zt * sqrt(tv) else NA_real_)
}

out <- rbindlist(lapply(SPEC, function(sp) {
  v <- sp$variable
  p0 <- pool_stratum(v, "g0")
  p1 <- pool_stratum(v, "g1")
  p_het <- NA_real_
  if (is.finite(p0$se) && is.finite(p1$se)) {
    delta <- p1$loghr - p0$loghr
    sed <- sqrt(p0$se^2 + p1$se^2)
    df_het <- min(p0$df, p1$df)
    p_het <- 2 * (if (is.finite(df_het)) pt(-abs(delta / sed), df_het)
                  else pnorm(-abs(delta / sed)))
  }
  rbind(data.table(variable = v, level = sp$labs[1],
                   HR = exp(p0$loghr),
                   HR_lci = exp(p0$lci), HR_uci = exp(p0$uci),
                   p_heterogeneity = p_het),
        data.table(variable = v, level = sp$labs[2],
                   HR = exp(p1$loghr),
                   HR_lci = exp(p1$lci), HR_uci = exp(p1$uci),
                   p_heterogeneity = p_het))
}))

cat("\n========== SUBGROUP multi-enrollment ==========\n")
print(out)
POOL_F <- file.path(OUT_T, paste0("v28.4_subgroup_", suffix, ".csv"))
fwrite(out, POOL_F)
cat("saved:", POOL_F, "\n")
cat("per_imp already appended into:", PER_F, "\n")
