# ============================================================================
# V26.9 PS 重叠度评估
# 生成倾向性评分分布重叠图，作为 Supplementary Figure S3
# 直接回应 structured review R01 要求（Causal Assumptions — positivity）
# ============================================================================

source("config.R")
source("_helpers_v26.R")
log_section("V26.9 Propensity Score Overlap Assessment")

# ---- 加载分析数据 ----
# 使用 augmented 数据集（含 DAG14 全部协变量）
aug_path <- DATA_AUGMENTED_FULL
if (!file.exists(aug_path)) {
  aug_path <- DATA_AUGMENTED
}
if (!file.exists(aug_path)) {
  log_error("增强数据集未找到，尝试主数据集。")
  aug_path <- file.path(DATA_DIR, "sequential_trials.parquet")
}
if (!file.exists(aug_path)) {
  log_error("所有数据路径均未找到，跳过。")
  quit(save = "no", status = 0)
}

log_info(sprintf("加载数据: %s", aug_path))
if (grepl("\\.parquet$", aug_path)) {
  suppressPackageStartupMessages(library(arrow))
  df <- arrow::read_parquet(aug_path)
} else if (grepl("\\.csv$", aug_path)) {
  suppressPackageStartupMessages(library(data.table))
  df <- as.data.frame(fread(aug_path))
} else {
  log_error(sprintf("不支持的文件格式: %s", aug_path))
  quit(save = "no", status = 0)
}

log_info(sprintf("数据加载完成: %d 行, %d 列", nrow(df), ncol(df)))

# ---- 构建 PS 模型（DAG14 协变量） ----
# 准备模型公式
treatment_var <- "itt_group"
df$treat_binary <- ifelse(df[[treatment_var]] == "treatment", 1L, 0L)

# 提取 DAG14 中存在于数据中的变量
dag14_available <- intersect(DAG14, names(df))
log_info(sprintf("DAG14 可用变量: %d/14", length(dag14_available)))

# 构建公式
formula_ps <- as.formula(
  paste("treat_binary ~", paste(dag14_available, collapse = " + "))
)

# 拟合 PS 模型（使用 na.action = na.exclude 确保预测长度为 nrow(df)）
log_info("拟合 PS 模型（DAG14 logistic 回归）...")
ps_model <- glm(formula_ps, data = df, family = binomial(link = "logit"),
                na.action = na.exclude)
log_info("PS 模型拟合完成。")

# 提取倾向性评分
df$propensity_score <- predict(ps_model, type = "response", na.action = na.exclude)

# 过滤 NA 值（因缺失协变量导致的预测缺失）
ps_valid <- df[!is.na(df$propensity_score), ]
log_info(sprintf("PS 有效行数: %d/%d (%.1f%%)", nrow(ps_valid), nrow(df),
                 100 * nrow(ps_valid) / nrow(df)))

# ---- 检查积极性（positivity） ----
# 按治疗组拆分 PS 分布
ps_treated <- na.omit(ps_valid$propensity_score[ps_valid$treat_binary == 1])
ps_control <- na.omit(ps_valid$propensity_score[ps_valid$treat_binary == 0])

log_info(sprintf("治疗组 PS: 均值=%.3f, SD=%.3f, 范围=[%.3f, %.3f]",
                 mean(ps_treated), sd(ps_treated), min(ps_treated), max(ps_treated)))
log_info(sprintf("对照组 PS: 均值=%.3f, SD=%.3f, 范围=[%.3f, %.3f]",
                 mean(ps_control), sd(ps_control), min(ps_control), max(ps_control)))

# 重叠度量化：计算两组 PS 的 Kolmogorov-Smirnov 统计量
ks_test <- ks.test(ps_treated, ps_control)
log_info(sprintf("KS 检验: D=%.3f, p=%s", ks_test$statistic, format(ks_test$p.value, digits=3)))

# 计算重叠区域比例（使用 min-max 方法）
overlap_min <- max(min(ps_treated), min(ps_control))
overlap_max <- min(max(ps_treated), max(ps_control))
overlap_range <- overlap_max - overlap_min
total_range <- max(max(ps_treated), max(ps_control)) - min(min(ps_treated), min(ps_control))
overlap_prop <- overlap_range / total_range
log_info(sprintf("PS 重叠区域比例: %.1f%%", overlap_prop * 100))

# ---- 生成重叠度图 ----
fig_file <- file.path(FIG_DIR, "S3_PS_Overlap.pdf")
pdf(fig_file, width = 10, height = 6)

# 密度图
par(mar = c(4.5, 4.5, 2, 1))
plot(density(ps_control), 
     col = "#0072B5", lwd = 2,
     xlim = c(0, 1), 
     xlab = "倾向性评分 (Propensity Score)",
     ylab = "密度",
     main = "Propensity Score Overlap: Treated vs. Control",
     cex.lab = 1.1, cex.main = 1.2)
lines(density(ps_treated), col = "#BC3C29", lwd = 2, lty = 2)
# 填充重叠区域
d_control <- density(ps_control)
d_treated <- density(ps_treated)
x_overlap <- seq(max(min(d_control$x), min(d_treated$x)),
                 min(max(d_control$x), max(d_treated$x)), length.out = 100)
y_control <- approx(d_control$x, d_control$y, x_overlap)$y
y_treated <- approx(d_treated$x, d_treated$y, x_overlap)$y
y_min <- pmin(y_control, y_treated, na.rm = TRUE)
polygon(c(x_overlap, rev(x_overlap)),
        c(y_min, rep(0, length(y_min))),
        col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
legend("topright",
       legend = c("对照组 (Control)", "治疗组 (Treated)"),
       col = c("#0072B5", "#BC3C29"),
       lwd = 2, lty = c(1, 2), bty = "n")
# 添加统计量注释
text(0.05, max(density(ps_control)$y) * 0.9,
     paste0("重叠区域: ", round(overlap_prop * 100), "%\n",
            "KS D = ", round(ks_test$statistic, 3)),
     adj = 0, cex = 0.9)

dev.off()
log_info(sprintf("PS 重叠度图已保存: %s", fig_file))

# ---- 生成重叠度统计表 ----
tbl_file <- file.path(TBL_DIR, "S3_PS_Overlap_Summary.csv")
overlap_summary <- data.frame(
  Metric = c("治疗组 PS 均值", "治疗组 PS 标准差", "治疗组 PS 最小值", "治疗组 PS 最大值",
             "对照组 PS 均值", "对照组 PS 标准差", "对照组 PS 最小值", "对照组 PS 最大值",
             "重叠区域比例", "KS D 统计量", "KS p 值"),
  Value = c(mean(ps_treated), sd(ps_treated), min(ps_treated), max(ps_treated),
            mean(ps_control), sd(ps_control), min(ps_control), max(ps_control),
            overlap_prop, ks_test$statistic, ks_test$p.value)
)
overlap_summary$Value <- round(overlap_summary$Value, 4)
write.csv(overlap_summary, tbl_file, row.names = FALSE)
log_info(sprintf("PS 重叠度统计表已保存: %s", tbl_file))

# ---- 极值检查 ----
near_zero <- sum(ps_valid$propensity_score < 0.05, na.rm = TRUE) / nrow(ps_valid) * 100
near_one  <- sum(ps_valid$propensity_score > 0.95, na.rm = TRUE) / nrow(ps_valid) * 100
log_info(sprintf("PS < 0.05: %.1f%% | PS > 0.95: %.1f%%", near_zero, near_one))

if (near_zero < 1 && near_one < 1) {
  log_info("积极性检查通过：无极端的 PS 值（<1% 在边缘区域）。")
} else {
  log_warn(sprintf("积极性警告：%.1f%% 的 PS 值在 [0, 0.05] 或 [0.95, 1] 区间。", 
                   near_zero + near_one))
}

log_info("V26.9 PS 重叠度评估完成。")
