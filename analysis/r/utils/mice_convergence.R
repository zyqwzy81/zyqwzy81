# ============================================================================
# V26.9 MICE 收敛诊断图
# 生成 MICE 链均值/标准差 trace plot，作为 Supplementary Figure S4
# 直接回应 structured review R05 要求
# ============================================================================

source("config.R")
source("_helpers_v26.R")
log_section("V26.9 MICE Convergence Diagnostics")

# ---- 加载 MICE 对象 ----
mice_path <- DATA_IMPUTED  # config.R: TTE/v26/output/rds/v26_mice_imputed.rds
if (!file.exists(mice_path)) {
  # 回退到备用路径
  mice_path <- file.path(DATA_DIR, "mice", "rds", "v26_mice_imp.rds")
}
if (!file.exists(mice_path)) {
  log_error(sprintf("MICE 对象未找到: %s", mice_path))
  log_info("请检查 MICE 对象路径。跳过收敛诊断图生成。")
  quit(save = "no", status = 0)
}

log_info(sprintf("加载 MICE 对象: %s", mice_path))
imp <- readRDS(mice_path)
log_info(sprintf("MICE 对象类型: %s", class(imp)[1]))
log_info(sprintf("插补集数: %d", imp$m))
log_info(sprintf("迭代次数: %d", imp$iteration))
log_info(sprintf("插补变量数: %d", length(imp$method)))

# ---- 生成收敛诊断图 ----
fig_file <- file.path(FIG_DIR, "S4_MICE_Convergence.pdf")
log_info(sprintf("输出路径: %s", fig_file))

# 选择关键变量绘制 trace plot（仅选择有实际链数据的连续型变量）
all_vars <- dimnames(imp$chainMean)[[1]]
# 检查哪些变量有非NA链数据
has_chain <- sapply(all_vars, function(v) {
  cm <- imp$chainMean[v, , ]
  !all(is.na(cm))
})
chain_vars <- all_vars[has_chain]
log_info(sprintf("有收敛链数据的变量: %d/%d", sum(has_chain), length(has_chain)))

# 优先选择临床重要的连续型变量
priority_vars <- c(
  "egfr", "lactate_t0", "bmi", "bilirubin_t0", "albumin_min",
  "hemoglobin_min", "bicarbonate_worst", "ph_worst",
  "platelets_min", "wbc_max", "inr_max", "ptt_max",
  "potassium_t0", "sodium_t0", "magnesium_t0", "bun_t0",
  "baseline_sofa", "sofa_t0"
)
key_vars <- intersect(priority_vars, chain_vars)
# 如果关键变量不够，添加其他可用变量
if (length(key_vars) < 5) {
  key_vars <- head(chain_vars, 10)
}
log_info(sprintf("将绘制 %d 个关键变量的收敛轨迹", min(length(key_vars), 10)))
# 限制最多10个变量，避免图形过大
key_vars <- head(key_vars, 10)

pdf(fig_file, width = 12, height = min(8, 2 * length(key_vars)))
par(mfrow = c(length(key_vars), 2), mar = c(2, 3, 1.5, 0.5))

for (v in key_vars) {
  # 提取链均值
  chain_mean <- imp$chainMean[v, , ]
  chain_var  <- imp$chainVar[v, , ]
  
  # 绘制均值 trace
  matplot(t(chain_mean), type = "l", lty = 1,
          col = c("#BC3C29", "#0072B5"),
          xlab = "迭代次数", ylab = "链均值",
          main = paste(v, "- 均值"))
  abline(h = mean(chain_mean, na.rm = TRUE), lty = 2, col = "gray60")
  
  # 绘制标准差 trace
  matplot(t(sqrt(chain_var)), type = "l", lty = 1,
          col = c("#BC3C29", "#0072B5"),
          xlab = "迭代次数", ylab = "链标准差",
          main = paste(v, "- 标准差"))
}

dev.off()
log_info(sprintf("收敛诊断图已保存: %s", fig_file))

# 检查收敛性（仅针对有链数据的变量）
conv_issues <- c()
for (v in key_vars) {
  cm <- imp$chainMean[v, , ]
  if (is.null(dim(cm)) || all(is.na(cm))) next
  if (nrow(cm) >= 3) {
    last_vals <- cm[nrow(cm), ]
    second_last <- cm[nrow(cm) - 1, ]
    rel_diff <- abs(last_vals - second_last) / (abs(second_last) + 0.001)
    if (any(rel_diff > 0.05, na.rm = TRUE)) {
      conv_issues <- c(conv_issues, v)
    }
  }
}

if (length(conv_issues) > 0) {
  log_warn(sprintf("以下变量可能存在收敛问题: %s", 
                   paste(conv_issues, collapse = ", ")))
} else {
  log_info("所有检查变量均显示良好收敛性。")
}

log_info("V26.9 MICE 收敛诊断完成。")
