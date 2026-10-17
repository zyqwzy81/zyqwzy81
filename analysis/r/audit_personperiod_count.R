# 核验：当前锁定队列（32,330 TP）的 person-period 行数各候选口径
# 无 bootstrap、无模型，仅构建日级表并计数
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")
suppressPackageStartupMessages(library(data.table))

cat("build event offsets...\n"); off <- build_event_offsets()
trials <- load_trials(1L)
daily  <- build_multi_daily(trials, off, H = TAU)

n_stays   <- uniqueN(trials$stay_id)
n_origins <- nrow(trials)
grid_rows <- nrow(daily)                       # 停留连续日网格行数
elig_rows <- sum(daily$eligible == 1L)         # 起点日行数（应=n_origins）

# SEQuential 克隆口径：每个起点展开 day0..28（29 行），网格保证 o+28 不越界
clone29 <- n_origins * 29L
clone28 <- n_origins * 28L

# 死亡截断口径：死亡停留中，对每个起点 o，窗口在 death_day 截断（含死亡日）
dd <- trials[, .(death_any = any(death_28day_num == 1L)), by = stay_id]
# 停留死亡日（与 build_multi_daily 同算法）
dday <- trials[death_28day_num == 1L,
               .(cand = pmin(TAU, pmax(0L, origin_day + as.integer(round(surv_time))))),
               by = stay_id][, .(dday = { tb <- sort(table(cand), decreasing = TRUE);
                                         as.integer(names(tb)[1]) }), by = stay_id]
od <- merge(trials[, .(stay_id, origin_day)], dday, by = "stay_id", all.x = TRUE)
od[is.na(dday), rows29 := 29L]
od[!is.na(dday), rows29 := pmin(29L, dday - origin_day + 1L)]
clone_deathcut <- sum(od$rows29)

cat("\n===== person-period count audit (current locked cohort) =====\n")
cat(sprintf("stays=%d  origins(trial-periods)=%d\n", n_stays, n_origins))
cat(sprintf("eligible rows on daily grid=%d (must equal origins)\n", elig_rows))
cat(sprintf("stay continuous-grid rows nrow(daily)=%d  (mean/stay=%.2f)\n",
            grid_rows, grid_rows / n_stays))
cat(sprintf("origins x 29 = %d (mean %.2f)\n", clone29, clone29 / n_origins))
cat(sprintf("origins x 28 = %d (mean %.2f)\n", clone28, clone28 / n_origins))
cat(sprintf("death-truncated clones = %d (mean %.2f)\n",
            clone_deathcut, clone_deathcut / n_origins))
cat("OLD manuscript claim: 825,532 rows on 32,328 origins, mean 25.5 days\n")
