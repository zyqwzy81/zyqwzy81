-- =====================================================================
-- M0_43: Safety electrolyte events within (T0, T0+48h]
-- K itemid 50971 mmol/L; Mg 50960 mg/dL (threshold 1.7 mg/dL = 0.7 mmol/L);
-- total Ca 50893 mg/dL x0.25 -> mmol/L
-- Phosphate >1.45 uses pure_control_phos_long
-- =====================================================================
SET work_mem='1GB';

DROP TABLE IF EXISTS public.pure_control_safety;
CREATE TABLE public.pure_control_safety AS
WITH elec AS (
  SELECT a.stay_id, l.itemid, l.v
  FROM public.pure_control_anchor a
  JOIN public.labs_key_electrolytes l ON l.hadm_id=a.hadm_id
  WHERE l.charttime > a.t0 AND l.charttime <= a.t0 + INTERVAL '48 hours'
    AND l.itemid IN (50971,50960,50893)
),
agg AS (
  SELECT stay_id,
    min(v) FILTER (WHERE itemid=50971)              AS k_min,
    min(v) FILTER (WHERE itemid=50960)              AS mg_min,
    min(v*0.25) FILTER (WHERE itemid=50893)         AS ca_min_mmol,
    count(*) FILTER (WHERE itemid=50971)            AS k_n,
    count(*) FILTER (WHERE itemid=50960)            AS mg_n,
    count(*) FILTER (WHERE itemid=50893)            AS ca_n
  FROM elec GROUP BY stay_id
),
phyp AS (
  SELECT stay_id, max(phos) AS phos_max48, count(*) phos_n48
  FROM public.pure_control_phos_long WHERE hph<=48 GROUP BY stay_id
)
SELECT a.stay_id,
  p.phos_max48, p.phos_n48,
  e.k_min, e.mg_min, e.ca_min_mmol, e.k_n, e.mg_n, e.ca_n,
  CASE WHEN p.phos_max48 > 1.45 THEN 1 ELSE 0 END        AS hyperphos,
  CASE WHEN e.k_min < 3.0  THEN 1 ELSE 0 END             AS hypoK,
  -- Mg 50960 is mg/dL (p50~2.0); 0.7 mmol/L = 1.70 mg/dL
  CASE WHEN e.mg_min < 1.7 THEN 1 ELSE 0 END             AS hypoMg,
  -- total (uncorrected) Ca mg/dL x0.25; <2.0 mmol/L = <8.0 mg/dL; subject to
  -- hypoalbuminaemia confounding in ICU, reported as sensitivity only
  CASE WHEN e.ca_min_mmol < 2.0 THEN 1 ELSE 0 END        AS hypoCa,
  CASE WHEN e.ca_min_mmol < 1.88 THEN 1 ELSE 0 END       AS hypoCa_severe
FROM public.pure_control_anchor a
LEFT JOIN phyp p ON a.stay_id=p.stay_id
LEFT JOIN agg  e ON a.stay_id=e.stay_id;

\echo '===== safety event rates by exposure (denominator = measured) ====='
SELECT e.exposure_w48,
  count(*) n,
  round(100.0*sum(s.hyperphos) FILTER (WHERE s.phos_n48>0)/nullif(sum(CASE WHEN s.phos_n48>0 THEN 1 END),0),1) pct_hyperphos,
  round(100.0*sum(s.hypoK)   FILTER (WHERE s.k_n>0)/nullif(sum(CASE WHEN s.k_n>0 THEN 1 END),0),1)  pct_hypoK,
  round(100.0*sum(s.hypoMg)  FILTER (WHERE s.mg_n>0)/nullif(sum(CASE WHEN s.mg_n>0 THEN 1 END),0),1) pct_hypoMg,
  round(100.0*sum(s.hypoCa)  FILTER (WHERE s.ca_n>0)/nullif(sum(CASE WHEN s.ca_n>0 THEN 1 END),0),1) pct_hypoCa
FROM public.pure_control_exposure e
JOIN public.pure_control_safety s ON e.stay_id=s.stay_id
WHERE e.exposure_w48 IN ('A_IV','B_PO_only','D_none')
GROUP BY 1 ORDER BY 1;

\echo '===== 48h normalization (>=0.80) among those measured at w48 ====='
SELECT e.exposure_w48, count(*) n_measured,
  round(100.0*sum(CASE WHEN w.phos>=0.80 THEN 1 ELSE 0 END)/count(*),1) pct_norm48
FROM public.pure_control_phos_windows w
JOIN public.pure_control_exposure e ON w.stay_id=e.stay_id
WHERE w.win='w48' AND e.exposure_w48 IN ('A_IV','B_PO_only','D_none')
GROUP BY 1 ORDER BY 1;
