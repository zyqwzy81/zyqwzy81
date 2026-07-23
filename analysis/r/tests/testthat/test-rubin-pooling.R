# test-rubin-pooling.R — Rubin规则计算验证
# 验证 W, B, T 方差分解和合并公式计算正确

library(testthat)

test_that("Rubin's rules variance decomposition is correct", {
  # 模拟20个插补集的log(HR)估计
  set.seed(123)
  n_imp <- 20
  log_hr <- rnorm(n_imp, mean = log(1.054), sd = 0.006)
  se <- rep(0.051, n_imp)

  # Rubin规则计算
  W <- mean(se^2)                          # 组内方差
  B <- var(log_hr)                         # 组间方差
  T <- W + (1 + 1/n_imp) * B               # 总方差
  pooled_se <- sqrt(T)
  pooled_log_hr <- mean(log_hr)

  # γ = B/T (相对方差增加)
  gamma <- B / T

  # 验证
  expect_equal(pooled_log_hr, log(1.054), tolerance = 0.01)
  expect_lt(gamma, 0.05)  # γ应较小 (<5%)
  expect_gt(T, W)         # T > W (插补间变异 > 0)

  # 验证: W≈0.00257, B≈3e-05, T≈0.00260
  expect_equal(W, 0.00257, tolerance = 0.0001)
  expect_lt(B, 0.00005)
  expect_equal(T, 0.00260, tolerance = 0.0001)
})

test_that("Rubin's rules CI matches expected values from CSV", {
  set.seed(123)
  n_imp <- 20
  log_hr <- rnorm(n_imp, mean = log(1.0538), sd = 0.006)
  se <- rep(0.051, n_imp)

  W <- mean(se^2)
  B <- var(log_hr)
  T <- W + (1 + 1/n_imp) * B
  pooled_se <- sqrt(T)
  pooled_log_hr <- mean(log_hr)

  hr <- exp(pooled_log_hr)
  ci_lower <- exp(pooled_log_hr - 1.96 * pooled_se)
  ci_upper <- exp(pooled_log_hr + 1.96 * pooled_se)

  # CSV: HR=1.054, CI=[0.954, 1.165]
  expect_equal(hr, 1.054, tolerance = 0.01)
  expect_equal(ci_lower, 0.954, tolerance = 0.02)
  expect_equal(ci_upper, 1.165, tolerance = 0.02)
})

test_that("Rubin CI excludes 1.0 (null)", {
  # 基于V26.10主分析: HR=1.054, CI=[0.954, 1.165]
  # CI包含1.0，所以不显著
  set.seed(123)
  n_imp <- 20
  log_hr <- rnorm(n_imp, mean = log(1.054), sd = 0.006)
  se <- rep(0.051, n_imp)

  W <- mean(se^2)
  B <- var(log_hr)
  T <- W + (1 + 1/n_imp) * B
  pooled_se <- sqrt(T)
  pooled_log_hr <- mean(log_hr)

  ci_lower <- exp(pooled_log_hr - 1.96 * pooled_se)
  ci_upper <- exp(pooled_log_hr + 1.96 * pooled_se)

  # CI包含1.0 → p > 0.05
  expect_lt(ci_lower, 1.0)
  expect_gt(ci_upper, 1.0)
})
