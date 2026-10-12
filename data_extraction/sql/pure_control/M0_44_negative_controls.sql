-- =====================================================================
-- M0_44: Negative controls
--  (1) Sodium 50983 post-T0 same 5 windows (plausible range 110-170)
--  (2) Pre-T0 phosphate slope: windows [-24,-12) c=-18, [-12,0) c=-6
--      strictly BEFORE t0 (charttime < t0)
-- =====================================================================
SET work_mem='1GB';

DROP TABLE IF EXISTS public.pure_control_sodium_windows;
CREATE TABLE public.pure_control_sodium_windows AS
WITH sl AS (
  SELECT a.stay_id, l.valuenum AS sodium,
         EXTRACT(EPOCH FROM (l.charttime-a.t0))/3600.0 hph
  FROM public.pure_control_anchor a
  JOIN mimiciv_hosp.labevents l ON l.hadm_id=a.hadm_id AND l.itemid=50983
  WHERE l.charttime > a.t0 AND l.charttime <= a.t0+INTERVAL '78 hours'
    AND l.valuenum BETWEEN 110 AND 170
), wins AS (
  SELECT 'w6' win,3 lo,9 hi,6 c UNION ALL SELECT 'w12',9,15,12 UNION ALL
  SELECT 'w24',18,30,24 UNION ALL SELECT 'w48',42,54,48 UNION ALL SELECT 'w72',66,78,72
), t AS (
  SELECT s.*, w.win, row_number() OVER (PARTITION BY s.stay_id,w.win ORDER BY abs(s.hph-w.c)) rn
  FROM sl s JOIN wins w ON s.hph>=w.lo AND s.hph<w.hi
)
SELECT stay_id,win,hph AS actual_hph,sodium FROM t WHERE rn=1;
CREATE INDEX idx_pcna_stay ON public.pure_control_sodium_windows(stay_id);

DROP TABLE IF EXISTS public.pure_control_prephos;
CREATE TABLE public.pure_control_prephos AS
WITH pl AS (
  SELECT a.stay_id, l.v*0.323 phos,
         EXTRACT(EPOCH FROM (l.charttime-a.t0))/3600.0 hph
  FROM public.pure_control_anchor a
  JOIN public.labs_key_electrolytes l ON l.hadm_id=a.hadm_id AND l.itemid=50970
  WHERE l.charttime >= a.t0-INTERVAL '24 hours' AND l.charttime < a.t0
    AND l.v*0.323 BETWEEN 0.2 AND 3.0
), wins AS (
  SELECT 'pre24' win,-24 lo,-12 hi,-18 c UNION ALL SELECT 'pre6',-12,0,-6
), t AS (
  SELECT p.*,w.win, row_number() OVER (PARTITION BY p.stay_id,w.win ORDER BY abs(p.hph-w.c)) rn
  FROM pl p JOIN wins w ON p.hph>=w.lo AND p.hph<w.hi
)
SELECT stay_id,win,hph AS actual_hph,phos FROM t WHERE rn=1;
CREATE INDEX idx_pcpre_stay ON public.pure_control_prephos(stay_id);

\echo '===== sodium window coverage ====='
SELECT e.exposure_w48, w.win, count(*) n, round(avg(w.sodium)::numeric,2) mean
FROM public.pure_control_exposure e
JOIN public.pure_control_sodium_windows w ON e.stay_id=w.stay_id
WHERE e.exposure_w48 IN ('A_IV','D_none') GROUP BY 1,2 ORDER BY 1,2;

\echo '===== pre-T0 phosphate (slope / parallel trend) ====='
SELECT e.exposure_w48, w.win, count(*) n,
  round(avg(w.phos)::numeric,3) mean,
  round(avg(w.actual_hph)::numeric,1) mean_h
FROM public.pure_control_exposure e
JOIN public.pure_control_prephos w ON e.stay_id=w.stay_id
WHERE e.exposure_w48 IN ('A_IV','D_none') GROUP BY 1,2 ORDER BY 1,2;
