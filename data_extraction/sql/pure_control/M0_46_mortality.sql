-- =====================================================================
-- M0_46: 28-day all-cause mortality from T0 (gold standard)
-- death date = COALESCE(admissions.deathtime::date, patients.dod)
-- follow to 28d regardless of discharge (out-of-hospital death via dod)
-- =====================================================================
SET work_mem='1GB';

DROP TABLE IF EXISTS public.pure_control_mortality;
CREATE TABLE public.pure_control_mortality AS
WITH icu AS (
  SELECT i.stay_id, i.subject_id, i.hadm_id, i.outtime
  FROM mimiciv_icu.icustays i
  WHERE i.stay_id IN (SELECT stay_id FROM public.pure_control_anchor)
),
adm AS (
  SELECT hadm_id, deathtime FROM mimiciv_hosp.admissions
),
pat AS (
  SELECT subject_id, dod FROM mimiciv_hosp.patients
),
j AS (
  SELECT a.stay_id, a.t0, icu.subject_id, icu.outtime,
    adm.deathtime, pat.dod,
    COALESCE(adm.deathtime::date, pat.dod) AS death_date
  FROM public.pure_control_anchor a
  JOIN icu ON a.stay_id=icu.stay_id
  LEFT JOIN adm ON icu.hadm_id=adm.hadm_id
  LEFT JOIN pat ON icu.subject_id=pat.subject_id
)
SELECT stay_id, subject_id, t0, dod, deathtime,
  CASE WHEN death_date IS NOT NULL
       THEN (death_date - t0::date) ELSE NULL END AS death_day,
  CASE WHEN death_date IS NOT NULL
        AND (death_date - t0::date) BETWEEN 0 AND 28
       THEN 1 ELSE 0 END AS death28,
  -- continuous follow-up time in days (min 1h so t>0)
  CASE WHEN death_date IS NOT NULL AND (death_date - t0::date) BETWEEN 0 AND 28
       THEN GREATEST(
              EXTRACT(EPOCH FROM (COALESCE(deathtime::timestamp,(death_date+INTERVAL '12 hours'))-t0))/86400.0,
              1/24.0)
       ELSE 28.0 END AS surv_days
FROM j;
CREATE INDEX idx_pcm_stay ON public.pure_control_mortality(stay_id);

\echo '===== C1/C3 28d mortality ====='
SELECT e.exposure_w48, count(*) n,
  sum(m.death28) deaths,
  round(100.0*avg(m.death28),2) pct
FROM public.pure_control_exposure e
JOIN public.pure_control_mortality m ON e.stay_id=m.stay_id
WHERE e.exposure_w48 IN ('A_IV','B_PO_only','D_none')
GROUP BY 1 ORDER BY 1;

\echo '===== cross-check: ALL 22685 anchor stays ====='
SELECT count(*) n, sum(death28) deaths, round(100.0*avg(death28),2) pct
FROM public.pure_control_mortality;
