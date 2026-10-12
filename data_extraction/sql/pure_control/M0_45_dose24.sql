-- M0_45: first-24h IV phosphate dose (mmol) per stay, cumulative at 48h
SET work_mem='512MB';
DROP TABLE IF EXISTS public.pure_control_dose;
CREATE TABLE public.pure_control_dose AS
SELECT a.stay_id,
  COALESCE(sum(i.mmol) FILTER (WHERE i.ev_time>a.t0 AND i.ev_time<=a.t0+INTERVAL '24 hours'),0) AS dose24,
  COALESCE(sum(i.mmol) FILTER (WHERE i.ev_time>a.t0 AND i.ev_time<=a.t0+INTERVAL '48 hours'),0) AS dose48,
  count(*) FILTER (WHERE i.ev_time>a.t0 AND i.ev_time<=a.t0+INTERVAL '48 hours') AS n_doses
FROM public.pure_control_anchor a
LEFT JOIN public.pure_control_iv_ie i ON a.stay_id=i.stay_id
GROUP BY a.stay_id;

\echo '===== dose24 distribution among A_IV (W48) ====='
SELECT e.exposure_w48, count(*) n,
  round(percentile_cont(0.25) WITHIN GROUP (ORDER BY d.dose24)::numeric,1) q1,
  round(percentile_cont(0.5) WITHIN GROUP (ORDER BY d.dose24)::numeric,1) median,
  round(percentile_cont(0.75) WITHIN GROUP (ORDER BY d.dose24)::numeric,1) q3,
  round(max(d.dose24)::numeric,1) max
FROM public.pure_control_dose d
JOIN public.pure_control_exposure e ON d.stay_id=e.stay_id
WHERE e.exposure_w48='A_IV' GROUP BY 1;
