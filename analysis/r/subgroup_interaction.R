# ============================================================================
# subgroup_interaction.R
# stay 唯一队列：18 亚组（9 变量×2 层）+ 10 乘法交互
# 稳定化 IPTW + coxph 稳健 SE + Rubin m=20（探索性，与 v28.3 方法一致，无 bootstrap）
# 版本: V28.4  日期: 2026-10-10
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/stay_common.R")
OUTDIR <- file.path(V28_BASE, "TTE/v28.4/output/tables/nonprimary_redo")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("loading", M_IMP, "imputations...\n")
DS <- lapply(seq_len(M_IMP), load_stay)
WS <- lapply(DS, stab_weights)

# 剔除层内零方差或造成完全分离（该水平内零事件 或 该水平外零事件）的设计列
clean_X <- function(X, event) {
  keep <- rep(TRUE, ncol(X))
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    if (is.na(var(x)) || var(x) == 0) { keep[j] <- FALSE; next }
    if (all(x %in% c(0, 1))) {
      ein  <- sum(event[x == 1])
      eout <- sum(event[x == 0])
      if (ein == 0 || eout == 0) keep[j] <- FALSE
    }
  }
  X[, keep, drop = FALSE]
}

# 通用加权 cox：返回指定系数名的 loghr/robust se
cox_est <- function(d, w, term, covs) {
  X <- clean_X(dag_matrix(d, covs), d$death_28day_num)
  df <- as.data.frame(X)
  df$arm <- d$arm; df$mod <- d[[term]]
  df$.t <- d$surv_time; df$.e <- d$death_28day_num
  covn <- setdiff(colnames(df), c("arm", "mod", ".t", ".e"))
  rhs <- paste(c("arm*mod", covn), collapse = "+")
  cx <- coxph(as.formula(paste("Surv(.t,.e)~", rhs)), data = df,
              weights = w, ties = "efron", robust = TRUE,
              control = coxph.control(iter.max = 50))
  s <- summary(cx)$coefficients
  rn <- if (is.matrix(s)) s["arm:mod", ] else s
  list(loghr = unname(rn["coef"]), se = unname(rn["se(coef)"]))
}

rubin <- function(q, se) {
  m <- length(q); qbar <- mean(q); ubar <- mean(se^2)
  bv <- if (m > 1) var(q) else 0
  tv <- ubar + bv * (1 + 1 / m)
  dfm <- if (m > 1 && bv > 0) (m - 1) * (1 + ubar / (bv * (1 + 1 / m)))^2 else Inf
  zt <- if (is.finite(dfm)) qt(.975, dfm) else qnorm(.975)
  p <- 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
            else pnorm(-abs(qbar / sqrt(tv))))
  list(q = qbar, L = qbar - zt * sqrt(tv), U = qbar + zt * sqrt(tv), p = p)
}

# ---- 亚组定义 ----
# labs[1]=g==0 层标签, labs[2]=g==1 层标签；grp 返回 0/1
subg_def <- list(
  list(v = "baseline_phosphate", excl = "baseline_phosphate",
       labs = c("Phosphate >= 0.50 mmol/L", "Phosphate < 0.50 mmol/L"),
       grp = function(d) as.integer(d$baseline_phosphate < 0.50)),
  list(v = "sepsis3", excl = "sepsis3",
       labs = c("Non Sepsis-3", "Sepsis-3"),
       grp = function(d) as.integer(d$sepsis3 == 1)),
  list(v = "baseline_sofa", excl = "baseline_sofa",
       labs = c("SOFA < 8", "SOFA >= 8"),
       grp = function(d) as.integer(d$baseline_sofa >= 8)),
  list(v = "age", excl = "age",
       labs = c("Age < 65 yr", "Age >= 65 yr"),
       grp = function(d) as.integer(d$age >= 65)),
  list(v = "gender", excl = "gender",
       labs = c("Female", "Male"),
       grp = function(d) as.integer(d$gender == "M")),
  list(v = "invasive_vent_at_t0", excl = "invasive_vent_at_t0",
       labs = c("No invasive ventilation", "Invasive ventilation"),
       grp = function(d) as.integer(d$invasive_vent_at_t0 == 1)),
  list(v = "vasoactive_at_t0", excl = "vasoactive_at_t0",
       labs = c("No vasoactive agents", "Vasoactive agents"),
       grp = function(d) as.integer(d$vasoactive_at_t0 == 1)),
  list(v = "egfr", excl = "egfr",
       labs = c("eGFR >= 60", "eGFR < 60"),
       grp = function(d) as.integer(d$egfr < 60)),
  list(v = "refeeding_risk", excl = "refeeding_risk",
       labs = c("No refeeding risk", "Refeeding risk"),
       grp = function(d) as.integer(d$refeeding_risk == 1)))

# 单层内加权 cox（arm + 协变量），稳健 SE
strat_cox <- function(d, w, covs) {
  X <- clean_X(dag_matrix(d, covs), d$death_28day_num)
  df <- as.data.frame(X)
  df$arm <- d$arm
  df$.t <- d$surv_time; df$.e <- d$death_28day_num
  covn <- setdiff(colnames(df), c("arm", ".t", ".e"))
  rhs <- paste(c("arm", covn), collapse = "+")
  cx <- coxph(as.formula(paste("Surv(.t,.e)~", rhs)), data = df,
              weights = w, ties = "efron", robust = TRUE,
              control = coxph.control(iter.max = 50))
  s <- summary(cx)$coefficients
  rn <- if (is.matrix(s)) s["arm", ] else s
  list(loghr = unname(rn["coef"]), se = unname(rn["se(coef)"]))
}

res_sub <- rbindlist(lapply(subg_def, function(sg) {
  rbindlist(lapply(0:1, function(gv) {
    q <- numeric(M_IMP); se <- numeric(M_IMP); nn <- 0L; ne <- 0L
    for (k in seq_len(M_IMP)) {
      g <- sg$grp(DS[[k]])
      ix <- which(g == gv)
      covs <- setdiff(DAG14, sg$excl)
      r <- strat_cox(DS[[k]][ix], WS[[k]][ix], covs)
      q[k] <- r$loghr; se[k] <- r$se
      if (k == 1L) { nn <- length(ix); ne <- sum(DS[[k]]$death_28day_num[ix]) }
    }
    rp <- rubin(q, se)
    data.table(subgroup = sg$labs[gv + 1L], n = nn, n_events = ne,
               hr = exp(rp$q), ci_lower = exp(rp$L), ci_upper = exp(rp$U),
               p_value = rp$p,
               significant = ifelse(rp$p < 0.05, "p<0.05", "NS"))
  }))
}))

# ---- 10 交互 ----
int_def <- list(
  c("vasoactive_at_t0", "Vasoactive agents"),
  c("baseline_sofa", "Baseline SOFA (per unit)"),
  c("lactate_t0", "Lactate (per unit)"),
  c("invasive_vent_at_t0", "Invasive ventilation"),
  c("charlson", "Charlson (per unit)"),
  c("age", "Age (per year)"),
  c("sepsis3", "Sepsis-3"),
  c("baseline_phosphate", "Baseline phosphate (per unit)"),
  c("refeeding_risk", "Refeeding risk"),
  c("egfr", "eGFR (per unit)"))

res_int <- rbindlist(lapply(int_def, function(sp) {
  q <- numeric(M_IMP); se <- numeric(M_IMP)
  for (k in seq_len(M_IMP)) {
    covs <- setdiff(DAG14, sp[1])
    r <- cox_est(DS[[k]], WS[[k]], sp[1], covs)
    q[k] <- r$loghr; se[k] <- r$se
  }
  rp <- rubin(q, se)
  data.table(variable = sp[1], label = sp[2],
             hr = exp(rp$q), ci_lower = exp(rp$L), ci_upper = exp(rp$U),
             p_value = rp$p)
}))
res_int[, q_value := p.adjust(p_value, method = "BH")]
res_int[, significant := ifelse(q_value < 0.10, "q<0.10", "NS")]

fwrite(res_sub, file.path(OUTDIR, "subgroup.csv"))
fwrite(res_int, file.path(OUTDIR, "interaction.csv"))
cat("\n=== Subgroups ===\n"); print(res_sub)
cat("\n=== Interactions ===\n"); print(res_int)
