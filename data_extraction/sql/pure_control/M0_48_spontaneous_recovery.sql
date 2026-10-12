-- D_none spontaneous phosphate recovery: mechanism-oriented descriptive analysis
\pset pager off
-- 0. checks: time_zero alignment & window labels
SELECT 'win labels' AS info, win, win_center, count(*) FROM public.pure_control_phos_windows GROUP BY win, win_center ORDER BY win_center;
SELECT count(*) AS nutrition_days,
       count(*) FILTER (WHERE n.time_zero = c.t0) AS aligned_t0
FROM public.nutrition_daily n
JOIN (SELECT stay_id, t0 FROM public.pure_control_cohort WHERE exposure_w48='D_none') c
  ON n.stay_id=c.stay_id;

DROP TABLE IF EXISTS tmp_dnone_mech;
CREATE TEMP TABLE tmp_dnone_mech AS
WITH d AS (
  SELECT stay_id, t0, icu_intime, baseline_phosphate, refeeding_risk, rrt_at_t0,
         nutrition_support_type, baseline_sofa, vasoactive_at_t0, lactate_t0, egfr
  FROM public.pure_control_cohort WHERE exposure_w48='D_none'
),
pw AS (
  SELECT stay_id,
         max(phos) FILTER (WHERE win='w48') AS phos48,
         max(phos) FILTER (WHERE win='w24') AS phos24,
         max(actual_hph) FILTER (WHERE win='w48') AS h48
  FROM public.pure_control_phos_windows WHERE win IN ('w24','w48') GROUP BY stay_id
),
nd0 AS (
  SELECT DISTINCT n.stay_id, n.nutrition_day, n.en_kcal, n.tpn_kcal, n.protein_g
  FROM public.nutrition_daily n
  JOIN (SELECT stay_id, t0 FROM public.pure_control_cohort WHERE exposure_w48='D_none') c
    ON n.stay_id=c.stay_id
   AND n.nutrition_day >= c.t0 AND n.nutrition_day < c.t0 + interval '48 hours'
),
nd AS (
  SELECT stay_id,
         coalesce(sum(en_kcal),0)  AS en_kcal_48,
         coalesce(sum(tpn_kcal),0) AS tpn_kcal_48,
         coalesce(sum(protein_g),0) AS protein_g_48
  FROM nd0 GROUP BY stay_id
),
cr AS (
  SELECT c.stay_id,
         avg(le.valuenum) FILTER (WHERE le.charttime >= c.t0 - interval '24 hours' AND le.charttime < c.t0) AS cr0,
         avg(le.valuenum) FILTER (WHERE le.charttime >= c.t0 + interval '24 hours' AND le.charttime < c.t0 + interval '48 hours') AS cr48
  FROM (SELECT stay_id, hadm_id, t0 FROM public.pure_control_cohort WHERE exposure_w48='D_none') c
  JOIN mimiciv_hosp.labevents le
    ON le.hadm_id=c.hadm_id AND le.charttime >= c.t0 - interval '24 hours' AND le.charttime < c.t0 + interval '48 hours'
  WHERE le.itemid=50912 AND le.valuenum IS NOT NULL AND le.valuenum>0 AND le.valuenum<30
  GROUP BY c.stay_id
),
hb AS (
  SELECT c.stay_id,
         avg(le.valuenum) FILTER (WHERE le.charttime >= c.t0 - interval '24 hours' AND le.charttime < c.t0) AS hb0,
         avg(le.valuenum) FILTER (WHERE le.charttime >= c.t0 + interval '24 hours' AND le.charttime < c.t0 + interval '48 hours') AS hb48
  FROM (SELECT stay_id, hadm_id, t0 FROM public.pure_control_cohort WHERE exposure_w48='D_none') c
  JOIN mimiciv_hosp.labevents le
    ON le.hadm_id=c.hadm_id AND le.charttime >= c.t0 - interval '24 hours' AND le.charttime < c.t0 + interval '48 hours'
  WHERE le.itemid=51222 AND le.valuenum IS NOT NULL AND le.valuenum>0 AND le.valuenum<25
  GROUP BY c.stay_id
),
pp AS (
  SELECT stay_id, avg(phos) AS pre24
  FROM public.pure_control_prephos GROUP BY stay_id
)
SELECT d.stay_id, d.baseline_phosphate, d.t0, d.icu_intime,
       extract(epoch FROM (d.t0-d.icu_intime))/3600 AS h_icu_to_t0,
       d.refeeding_risk, d.rrt_at_t0, d.baseline_sofa, d.vasoactive_at_t0, d.egfr,
       pw.phos24, pw.phos48, (pw.phos48-d.baseline_phosphate) AS d48,
       CASE WHEN pw.phos48>=0.80 THEN 1 ELSE 0 END AS norm48,
       nd.en_kcal_48, nd.tpn_kcal_48, nd.protein_g_48,
       CASE WHEN nd.en_kcal_48>0 AND nd.tpn_kcal_48=0 THEN 'EN'
            WHEN nd.tpn_kcal_48>0 AND nd.en_kcal_48=0 THEN 'TPN'
            WHEN nd.en_kcal_48>0 AND nd.tpn_kcal_48>0 THEN 'EN+TPN'
            ELSE 'none' END AS feed48,
       cr.cr0, cr.cr48, (cr.cr48-cr.cr0) AS dcr,
       hb.hb0, hb.hb48, (hb.hb48-hb.hb0) AS dhb,
       pp.pre24
FROM d
LEFT JOIN pw ON d.stay_id=pw.stay_id
LEFT JOIN nd ON d.stay_id=nd.stay_id
LEFT JOIN cr ON d.stay_id=cr.stay_id
LEFT JOIN hb ON d.stay_id=hb.stay_id
LEFT JOIN pp ON d.stay_id=pp.stay_id;

-- population denominators
SELECT count(*) n_dnone,
       count(phos48) n_p48,
       round(100.0*count(phos48)/count(*),1) pct_meas48,
       round(avg(baseline_phosphate)::numeric,3) base_mean,
       round(avg(phos48)::numeric,3) p48_mean,
       round(avg(d48)::numeric,3) d48_mean,
       round(100.0*avg(norm48)::numeric,1) norm_pct
FROM tmp_dnone_mech;

-- A. setpoint rebound (pre-episode vs 48h among those with both)
SELECT count(*) n,
       round(avg(pre24)::numeric,3) pre24,
       round(avg(baseline_phosphate)::numeric,3) t0,
       round(avg(phos48)::numeric,3) p48,
       round(avg(phos48-pre24)::numeric,3) gap48
FROM tmp_dnone_mech WHERE pre24 IS NOT NULL AND phos48 IS NOT NULL;

-- B1. feeding exposure prevalence (0-48h)
SELECT feed48, count(*) n,
       round(100.0*count(*)/sum(count(*)) OVER (),1) pct,
       round(avg(en_kcal_48+tpn_kcal_48)::numeric,0) kcal_mean,
       round(avg(protein_g_48)::numeric,1) protein_mean
FROM tmp_dnone_mech GROUP BY feed48 ORDER BY n DESC;

-- B2. recovery by feeding (measured p48)
SELECT feed48, count(*) n,
       round(avg(d48)::numeric,3) d48_mean,
       stddev(d48)::numeric(6,3) d48_sd,
       round(100.0*avg(norm48)::numeric,1) norm_pct,
       round(avg(phos48)::numeric,3) p48_mean
FROM tmp_dnone_mech WHERE phos48 IS NOT NULL GROUP BY feed48 ORDER BY d48_mean DESC;

-- B3. recovery by kcal adequacy tertiles among fed
SELECT width_bucket(en_kcal_48+tpn_kcal_48, ARRAY[0,200,800,100000]) AS kcal_bucket,
       count(*) n, round(avg(d48)::numeric,3) d48, round(100.0*avg(norm48)::numeric,1) norm_pct
FROM tmp_dnone_mech WHERE phos48 IS NOT NULL AND (en_kcal_48+tpn_kcal_48)>0
GROUP BY 1 ORDER BY 1;

-- C. renal: dCr strata vs phosphate recovery (no RRT at t0)
SELECT CASE WHEN dcr <= -0.1 THEN '1_improved' WHEN dcr >= 0.1 THEN '3_worsened' ELSE '2_stable' END cr_stratum,
       count(*) n, round(avg(dcr)::numeric,3) dcr_mean,
       round(avg(d48)::numeric,3) d48, round(100.0*avg(norm48)::numeric,1) norm_pct
FROM tmp_dnone_mech WHERE phos48 IS NOT NULL AND dcr IS NOT NULL AND rrt_at_t0=0
GROUP BY 1 ORDER BY 1;
-- correlation dCr vs d48
SELECT round(corr(dcr,d48)::numeric,3) corr_dcr_d48, count(*) n FROM tmp_dnone_mech WHERE dcr IS NOT NULL AND d48 IS NOT NULL;
-- AKI recovery: t0 eGFR strata
SELECT CASE WHEN egfr>=60 THEN '1_egfr>=60' WHEN egfr>=30 THEN '2_30-59' WHEN egfr IS NULL THEN '4_unknown' ELSE '3_<30' END egfr_stratum,
       count(*) n, round(avg(d48)::numeric,3) d48, round(100.0*avg(norm48)::numeric,1) norm_pct
FROM tmp_dnone_mech WHERE phos48 IS NOT NULL GROUP BY 1 ORDER BY 1;

-- D1. hemodilution proxy: Hb change strata vs phosphate recovery
SELECT CASE WHEN dhb >= 0.5 THEN '1_rise>=0.5' WHEN dhb <= -0.5 THEN '3_drop<=-0.5' ELSE '2_flat' END hb_stratum,
       count(*) n, round(avg(dhb)::numeric,3) dhb_mean,
       round(avg(d48)::numeric,3) d48, round(100.0*avg(norm48)::numeric,1) norm_pct
FROM tmp_dnone_mech WHERE phos48 IS NOT NULL AND dhb IS NOT NULL GROUP BY 1 ORDER BY 1;
SELECT round(corr(dhb,d48)::numeric,3) corr_dhb_d48, count(*) n FROM tmp_dnone_mech WHERE dhb IS NOT NULL AND d48 IS NOT NULL;

-- D2. timing: early-ICU vs later-ICU eligible (acute transient dip hypothesis)
SELECT CASE WHEN h_icu_to_t0 <= 24 THEN '1_<=24h' WHEN h_icu_to_t0 <= 72 THEN '2_24-72h' ELSE '3_>72h' END icu_timing,
       count(*) n, round(avg(h_icu_to_t0)::numeric,1) h_mean,
       round(avg(baseline_phosphate)::numeric,3) base,
       round(avg(d48)::numeric,3) d48, round(100.0*avg(norm48)::numeric,1) norm_pct
FROM tmp_dnone_mech WHERE phos48 IS NOT NULL GROUP BY 1 ORDER BY 1;

-- D3. baseline severity strata & vasopressor resolution proxy via sofa
SELECT CASE WHEN baseline_sofa<=5 THEN '1_sofa<=5' WHEN baseline_sofa<=10 THEN '2_6-10' ELSE '3_>10' END sofa_stratum,
       count(*) n, round(avg(d48)::numeric,3) d48, round(100.0*avg(norm48)::numeric,1) norm_pct
FROM tmp_dnone_mech WHERE phos48 IS NOT NULL AND baseline_sofa IS NOT NULL GROUP BY 1 ORDER BY 1;

-- baseline phosphate strata (regression-to-mean / severity)
SELECT CASE WHEN baseline_phosphate<0.5 THEN '1_<0.5' WHEN baseline_phosphate<0.65 THEN '2_0.5-0.65' ELSE '3_0.65-0.8' END base_stratum,
       count(*) n, round(avg(baseline_phosphate)::numeric,3) base,
       round(avg(d48)::numeric,3) d48, round(100.0*avg(norm48)::numeric,1) norm_pct,
       round(avg(phos48)::numeric,3) p48
FROM tmp_dnone_mech WHERE phos48 IS NOT NULL GROUP BY 1 ORDER BY 1;

-- crude multivariable: OLS of d48 on factors (scaled coefficients)
SELECT 'OLS' AS model,
       round(regr_slope(d48,(en_kcal_48+tpn_kcal_48)/100.0)::numeric,4) AS beta_kcal_per100,
       round(regr_slope(d48,dcr)::numeric,4) AS beta_dcr,
       round(regr_slope(d48,dhb)::numeric,4) AS beta_dhb,
       round(regr_slope(d48,h_icu_to_t0/24.0)::numeric,4) AS beta_per_day_icu,
       round(regr_slope(d48,baseline_sofa)::numeric,4) AS beta_sofa,
       round(regr_slope(d48,baseline_phosphate)::numeric,4) AS beta_base,
       count(*) n
FROM tmp_dnone_mech
WHERE d48 IS NOT NULL AND dcr IS NOT NULL AND dhb IS NOT NULL AND baseline_sofa IS NOT NULL;

-- EN phosphate supply estimate: among EN-fed, kcal distribution & implied P intake (0.75 mg/kcal standard polymeric)
SELECT count(*) n_fed,
       round(avg(en_kcal_48+tpn_kcal_48)::numeric,0) kcal_mean,
       round(percentile_cont(0.5) WITHIN GROUP (ORDER BY en_kcal_48+tpn_kcal_48)::numeric,0) kcal_median,
       round((sum(en_kcal_48+tpn_kcal_48)/48/20/1000.0)::numeric,2) AS implied_total_kcals_k
FROM tmp_dnone_mech WHERE (en_kcal_48+tpn_kcal_48)>0;
