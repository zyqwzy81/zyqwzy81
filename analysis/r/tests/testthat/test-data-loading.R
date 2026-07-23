# test-data-loading.R — 数据加载测试
# 验证核心Parquet文件可加载且结构符合预期

library(testthat)
library(arrow)

test_that("sequential_trials.parquet exists and is loadable", {
  parquet_path <- file.path("TTE/v26/output/sequential_trials.parquet")
  expect_true(file.exists(parquet_path))
  df <- read_parquet(parquet_path)
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 30000)  # 预期 >32,000
})

test_that("core columns are present with correct types", {
  parquet_path <- file.path("TTE/v26/output/sequential_trials.parquet")
  df <- read_parquet(parquet_path)

  required_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
                     "itt_group", "baseline_phosphate", "death_28day_num",
                     "surv_time", "age", "gender")
  for (col in required_cols) {
    expect_true(col %in% names(df), info = paste("Missing column:", col))
  }
})

test_that("itt_group has expected values", {
  parquet_path <- file.path("TTE/v26/output/sequential_trials.parquet")
  df <- read_parquet(parquet_path)
  expect_setequal(unique(df$itt_group), c("treatment", "conservative"))
})

test_that("surv_time is non-negative and death_28day_num is binary", {
  parquet_path <- file.path("TTE/v26/output/sequential_trials.parquet")
  df <- read_parquet(parquet_path)
  expect_true(all(df$surv_time >= 0, na.rm = TRUE))
  expect_setequal(unique(df$death_28day_num[!is.na(df$death_28day_num)]), c(0, 1))
})

test_that("N and events match expected values", {
  parquet_path <- file.path("TTE/v26/output/sequential_trials.parquet")
  df <- read_parquet(parquet_path)
  expect_equal(nrow(df), 32325, tolerance = 50)
  n_events <- sum(df$death_28day_num, na.rm = TRUE)
  expect_equal(n_events, 5610, tolerance = 50)
})
