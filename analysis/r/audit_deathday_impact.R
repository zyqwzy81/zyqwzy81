# 核验：修正 stay_death_day（绝对日历日期共识 vs 旧 round 众数）的影响面
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_common.R")
suppressPackageStartupMessages(library(data.table))

trials <- load_trials(1L)
mode1 <- function(x) { tb <- sort(table(x), decreasing = TRUE); as.integer(names(tb)[1]) }

# 旧规则：origin_day + round(surv) 众数
old_dd <- trials[death_28day_num == 1L,
    .(old_day = mode1(pmin(TAU, pmax(0L, origin_day + as.integer(round(surv_time)))))),
    by = stay_id]

# 新规则：绝对时刻加 surv → 死亡日历日 → 众数，再相对最早入组日
trials[, cand_date := as.Date(time_zero + surv_time * 86400)]
new_dd <- trials[death_28day_num == 1L,
    .(new_day = as.integer(mode1(as.integer(cand_date)) -
                           min(as.integer(as.Date(time_zero))))),
    by = stay_id]

cmp <- merge(old_dd, new_dd, by = "stay_id")
cmp[, delta := new_day - old_day]
cat("dead stays:", nrow(cmp), "\n")
cat("death_day changed stays:", nrow(cmp[delta != 0]), "\n")
print(cmp[, .N, by = delta][order(delta)])

# 对入模保留/剔除的影响
dt <- merge(trials, old_dd, by = "stay_id", all.x = TRUE)
kept_old <- dt[is.na(old_day) | origin_day <= old_day, .N]
dt2 <- merge(trials, new_dd, by = "stay_id", all.x = TRUE)
kept_new <- dt2[is.na(new_day) | origin_day <= new_day, .N]
cat("origins kept old/new:", kept_old, "/", kept_new, "\n")
cat("TP deaths kept old/new:",
    dt[is.na(old_day) | origin_day <= old_day, sum(death_28day_num)], "/",
    dt2[is.na(new_day) | origin_day <= new_day, sum(death_28day_num)], "\n")

# 两个目标停留新规则下的死亡日与保留
cat("\ntarget stays:\n")
print(new_dd[stay_id %in% c(36046816, 39360004)])
print(trials[stay_id %in% c(36046816, 39360004),
      .(stay_id, trial_num, time_zero, origin_day, death_28day_num, surv_time)][order(stay_id,trial_num)])

# 变化最大的停留示例（delta 绝对值最大的前10）
cmp[, absd := abs(delta)]
print(cmp[absd > 0][order(-absd)][1:10])
