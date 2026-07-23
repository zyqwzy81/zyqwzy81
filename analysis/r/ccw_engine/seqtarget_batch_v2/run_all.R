#!/usr/bin/env Rscript
# ============================================================================
# run_all.R — SEQTaRget 兼容分批引擎 V2 完整管线
#
# 运行顺序:
#   01_pooled_ps.R        Phase 1: Pooled PS + IPTW（全量一次性拟合）
#   02_batched_seqexpand.R  Phase 2: 按 period 分批日粒度展开
#   03_batched_logistic.R   Phase 3: 按 period 累加式 pooled logistic
#   04_quality_check.R      Phase 4: 质量检查 + 对比验证
#
# 每个 phase 可独立运行（读取上一 phase 的输出），支持断点续跑。
#
# 使用方法:
#   cd TTE/v26/code/r/seqtarget_batch_v2/
#   Rscript run_all.R
# ============================================================================

message("")
message("================================================================")
message("  SEQTaRget 兼容分批引擎 V2 — 完整管线")
message("================================================================")
message("")
message(sprintf("  启动时间: %s", Sys.time()))
message("")

script_dir <- if (interactive()) getwd() else dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
setwd(script_dir)
message(sprintf("  工作目录: %s", getwd()))

# ---- Phase 1: Pooled PS + IPTW ----
message("\n>>> Phase 1: Pooled PS + IPTW <<<")
t1 <- Sys.time()
source("01_pooled_ps.R")
t1_elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
message(sprintf("Phase 1 完成: %.1f 秒", t1_elapsed))

# ---- Phase 2: Batched SEQexpand ----
message("\n>>> Phase 2: Batched SEQexpand <<<")
t2 <- Sys.time()
source("02_batched_seqexpand.R")
t2_elapsed <- as.numeric(difftime(Sys.time(), t2, units = "secs"))
message(sprintf("Phase 2 完成: %.1f 秒", t2_elapsed))

# ---- Phase 3: Batched Pooled Logistic ----
message("\n>>> Phase 3: Batched Pooled Logistic <<<")
t3 <- Sys.time()
source("03_batched_logistic.R")
t3_elapsed <- as.numeric(difftime(Sys.time(), t3, units = "secs"))
message(sprintf("Phase 3 完成: %.1f 秒", t3_elapsed))

# ---- Phase 4: Quality Check ----
message("\n>>> Phase 4: Quality Check <<<")
t4 <- Sys.time()
source("04_quality_check.R")
t4_elapsed <- as.numeric(difftime(Sys.time(), t4, units = "secs"))
message(sprintf("Phase 4 完成: %.1f 秒", t4_elapsed))

# ---- 汇总 ----
total_time <- t1_elapsed + t2_elapsed + t3_elapsed + t4_elapsed
message("")
message("================================================================")
message("  管线执行完成")
message(sprintf("  总运行时间: %.1f 秒 (%.2f 分钟)", total_time, total_time / 60))
message(sprintf("  完成时间: %s", Sys.time()))
message("================================================================")
message("")
message("输出目录:")
message(sprintf("  Phase 1 (权重): %s", file.path("..", "..", "..", "output", "seqtarget_batch_v2", "weighted")))
message(sprintf("  Phase 2 (展开): %s", file.path("..", "..", "..", "output", "seqtarget_batch_v2", "expanded")))
message(sprintf("  Phase 3-4 (结果): %s", file.path("..", "..", "..", "output", "seqtarget_batch_v2", "results")))
message("")
