-- =====================================================================
-- M0_42: Post-T0 phosphate measurements (long) + nearest-value windows
-- Unit mmol/L (itemid 50970 mg/dL x 0.323); strictly AFTER t0; 0.2-3.0
-- =====================================================================
SET work_mem='1GB';

-- 1. Long table: every valid post-T0 phosphate measurement within 78h
DROP TABLE IF EXISTS public.pure_control_phos_long;
CREATE TABLE public.pure_control_phos_long AS
SELECT a.stay_id, e.exposure_w48, a.t0,
       l.charttime,
       EXTRACT(EPOCH FROM (l.charttime - a.t0))/3600.0 AS hph,
       round(l.v*0.323, 3) AS phos
FROM public.pure_control_anchor a
JOIN public.pure_control_exposure e ON a.stay_id=e.stay_id
JOIN public.labs_key_electrolytes l ON l.hadm_id=a.hadm_id AND l.itemid=50970
WHERE l.charttime > a.t0
  AND l.charttime <= a.t0 + INTERVAL '78 hours'
  AND l.v*0.323 BETWEEN 0.2 AND 3.0;
CREATE INDEX idx_pcpl_stay ON public.pure_control_phos_long(stay_id);
CREATE INDEX idx_pcpl_hph ON public.pure_control_phos_long(hph);

-- 2. Window table: one nearest-to-center value per stay per window
DROP TABLE IF EXISTS public.pure_control_phos_windows;
CREATE TABLE public.pure_control_phos_windows AS
WITH wins AS (
  SELECT 'w6'  win,  3 lo,  9 hi,  6 center UNION ALL
  SELECT 'w12',      9,    15,    12 UNION ALL
  SELECT 'w24',     18,    30,    24 UNION ALL
  SELECT 'w48',     42,    54,    48 UNION ALL
  SELECT 'w72',     66,    78,    72
),
tagged AS (
  SELECT p.*, w.win, w.center,
         row_number() OVER (PARTITION BY p.stay_id, w.win
                            ORDER BY abs(p.hph - w.center)) rn
  FROM public.pure_control_phos_long p
  JOIN wins w ON p.hph >= w.lo AND p.hph < w.hi
)
SELECT stay_id, exposure_w48, win, center AS win_center,
       hph AS actual_hph, phos
FROM tagged WHERE rn=1;
CREATE INDEX idx_pcpw_stay ON public.pure_control_phos_windows(stay_id);

\echo '===== long table ====='
SELECT count(*) n_rows, count(DISTINCT stay_id) n_stays FROM public.pure_control_phos_long;

\echo '===== window coverage by exposure group ====='
SELECT e.exposure_w48,
  count(*) AS group_n,
  count(DISTINCT w6.stay_id)  AS n6,
  count(DISTINCT w12.stay_id) AS n12,
  count(DISTINCT w24.stay_id) AS n24,
  count(DISTINCT w48.stay_id) AS n48,
  count(DISTINCT w72.stay_id) AS n72
FROM public.pure_control_exposure e
LEFT JOIN (SELECT stay_id FROM public.pure_control_phos_windows WHERE win='w6') w6 ON e.stay_id=w6.stay_id
LEFT JOIN (SELECT stay_id FROM public.pure_control_phos_windows WHERE win='w12') w12 ON e.stay_id=w12.stay_id
LEFT JOIN (SELECT stay_id FROM public.pure_control_phos_windows WHERE win='w24') w24 ON e.stay_id=w24.stay_id
LEFT JOIN (SELECT stay_id FROM public.pure_control_phos_windows WHERE win='w48') w48 ON e.stay_id=w48.stay_id
LEFT JOIN (SELECT stay_id FROM public.pure_control_phos_windows WHERE win='w72') w72 ON e.stay_id=w72.stay_id
WHERE e.exposure_w48 IN ('A_IV','B_PO_only','D_none')
GROUP BY 1 ORDER BY 1;
