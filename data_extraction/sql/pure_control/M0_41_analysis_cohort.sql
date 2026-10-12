-- =====================================================================
-- M0_41: Flat analysis cohort — exposure + DAG14 covariates (stay-level)
-- Matches dim tables on FIRST enrollment time_zero (= t0)
-- =====================================================================
SET work_mem='1GB';

DROP TABLE IF EXISTS public.pure_control_cohort;
CREATE TABLE public.pure_control_cohort AS
SELECT e.*,
       d1.age, d1.gender, d1.bmi,
       d3.sepsis3, d3.invasive_vent_at_t0 AS vent_dim3,
       d4.charlson_comorbidity_index AS charlson,
       d5.egfr, d5.lactate_t0, d5.sofa_t0 AS baseline_sofa,
       d5.vasoactive_at_t0, d5.sodium_t0, d5.potassium_t0,
       d6.corrected_calcium_mmol_min,
       sc.admission_type
FROM public.pure_control_exposure e
LEFT JOIN public.baseline_dim1_demographics d1
       ON e.stay_id=d1.stay_id AND d1.time_zero=e.t0
LEFT JOIN public.baseline_dim3_disease_severity d3
       ON e.stay_id=d3.stay_id AND d3.time_zero=e.t0
LEFT JOIN public.baseline_dim4_comorbidities d4
       ON e.stay_id=d4.stay_id AND d4.time_zero=e.t0
LEFT JOIN public.baseline_dim5_physio_organ d5
       ON e.stay_id=d5.stay_id AND d5.time_zero=e.t0
LEFT JOIN public.baseline_dim6_hypo_electrolyte d6
       ON e.stay_id=d6.stay_id AND d6.time_zero=e.t0
LEFT JOIN (SELECT DISTINCT ON (stay_id) stay_id, admission_type
           FROM public.study_cohort ORDER BY stay_id, time_zero) sc
       ON e.stay_id=sc.stay_id;

CREATE INDEX idx_pcc_stay ON public.pure_control_cohort(stay_id);

\echo '===== cohort rows (must 22,685) ====='
SELECT count(*) FROM public.pure_control_cohort;

\echo '===== covariate availability in C1 (A_IV vs D_none) ====='
SELECT exposure_w48, count(*) n,
  round(100.0*count(age)/count(*),1) age_pct,
  round(100.0*count(baseline_sofa)/count(*),1) sofa_pct,
  round(100.0*count(egfr)/count(*),1) egfr_pct,
  round(100.0*count(lactate_t0)/count(*),1) lactate_pct,
  round(100.0*count(corrected_calcium_mmol_min)/count(*),1) ca_pct,
  round(100.0*count(charlson)/count(*),1) charlson_pct,
  round(100.0*count(sepsis3)/count(*),1) sepsis_pct
FROM public.pure_control_cohort
WHERE exposure_w48 IN ('A_IV','D_none','B_PO_only')
GROUP BY 1 ORDER BY 1;

\echo '===== variable types quick check ====='
SELECT exposure_w48,
  round(avg(age)::numeric,1) mean_age,
  string_agg(DISTINCT gender, '/') genders,
  round(avg(baseline_phosphate)::numeric,3) mean_p,
  round(avg(baseline_sofa)::numeric,1) mean_sofa,
  round(100.0*avg(sepsis3::int),1) sepsis_pct,
  round(100.0*avg(vent_dim3::int),1) vent_pct,
  round(100.0*avg(vasoactive_at_t0::int),1) vaso_pct
FROM public.pure_control_cohort
WHERE exposure_w48 IN ('A_IV','D_none','B_PO_only')
GROUP BY 1 ORDER BY 1;
