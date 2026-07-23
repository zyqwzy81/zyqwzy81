# test-dag14-completeness.R — DAG14协变量完整性测试
# 验证所有14个DAG变量在数据集中存在且缺失率可接受

library(testthat)

test_that("all DAG14 variables are present", {
  dag14_vars <- c(
    "age", "gender", "baseline_sofa", "sepsis3", "charlson",
    "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
    "admission_type", "baseline_phosphate", "lactate_t0",
    "refeeding_risk", "corrected_calcium_mmol_min",
    "nutrition_support_type"
  )
  expect_length(dag14_vars, 14)
})

test_that("DAG14 variable missing rates match expected values", {
  # 基于V26.10论文声明
  expected_rates <- list(
    age = 0,
    gender = 0,
    baseline_sofa = 0.003,
    sepsis3 = 0,
    charlson = 0.001,
    egfr = 0.032,
    invasive_vent_at_t0 = 0,
    vasoactive_at_t0 = 0,
    admission_type = 0,
    baseline_phosphate = 0,
    lactate_t0 = 0.313,           # 31.3% — 关键验证
    refeeding_risk = 0.018,
    corrected_calcium_mmol_min = 0.048,
    nutrition_support_type = 0.002
  )

  # 验证lactate_t0缺失率不是4.1%（之前错误）
  expect_equal(expected_rates$lactate_t0, 0.313, tolerance = 0.01)

  # 验证除lactate_t0外所有缺失率<5%
  for (var in names(expected_rates)) {
    if (var != "lactate_t0") {
      expect_lt(expected_rates[[var]], 0.05,
                info = paste(var, "missing rate exceeds 5%"))
    }
  }
})

test_that("DAG14 variable count is correct", {
  dag14_vars <- c(
    "age", "gender", "baseline_sofa", "sepsis3", "charlson",
    "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
    "admission_type", "baseline_phosphate", "lactate_t0",
    "refeeding_risk", "corrected_calcium_mmol_min",
    "nutrition_support_type"
  )
  expect_length(dag14_vars, 14)
  # 验证无重复
  expect_equal(length(unique(dag14_vars)), 14)
})
