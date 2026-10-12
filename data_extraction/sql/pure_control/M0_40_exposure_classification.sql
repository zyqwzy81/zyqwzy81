-- =====================================================================
-- M0_40: Pure-control exposure classification (post-hoc mechanistic)
-- Spec: .trae/specs/iv-vs-pure-control-phosphate/spec.md (§3-§4)
-- Unit: one row per unique ICU stay (n=22,685), T0 = first enrollment
-- Windows: W12=(T0,+12h], W48=(T0,+48h] (primary), W72=(T0,+72h]
-- =====================================================================
SET work_mem = '1GB';

-- ---------------------------------------------------------------------
-- 0. Anchor: one row per stay, T0 = first enrollment time_zero
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS public.pure_control_anchor;
CREATE TABLE public.pure_control_anchor AS
SELECT stay_id,
       min(hadm_id)    AS hadm_id,
       min(time_zero)  AS t0,
       min(icu_intime) AS icu_intime,
       max(icu_outtime) AS icu_outtime,
       min(baseline_phosphate) AS baseline_phosphate
FROM public.study_cohort
GROUP BY stay_id;
CREATE INDEX idx_pca_stay ON public.pure_control_anchor(stay_id);

-- ---------------------------------------------------------------------
-- 1. IV phosphate ACTUAL INFUSIONS — inputevents (stay-level, mmol)
--    itemid 225834 "K Phos" + 225835 "Na Phos" (both actual administrations)
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS public.pure_control_iv_ie;
CREATE TABLE public.pure_control_iv_ie AS
SELECT ie.stay_id, ie.starttime AS ev_time, ie.amount::numeric AS mmol,
       CASE WHEN ie.itemid=225834 THEN 'K Phos' ELSE 'Na Phos' END AS drug
FROM mimiciv_icu.inputevents ie
WHERE ie.itemid IN (225834, 225835)
  AND ie.amount BETWEEN 0.006 AND 500
  AND ie.starttime IS NOT NULL;
CREATE INDEX idx_pcie_stay ON public.pure_control_iv_ie(stay_id);
CREATE INDEX idx_pcie_time ON public.pure_control_iv_ie(ev_time);

-- ---------------------------------------------------------------------
-- 2. IV phosphate events — prescriptions (hadm-level, no reliable dose)
--    Whitelist: Potassium/Sodium/Sodium Glycerophosphate, IV routes only
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS public.pure_control_iv_rx;
CREATE TABLE public.pure_control_iv_rx AS
SELECT a.stay_id, p.starttime AS ev_time
FROM public.pure_control_anchor a
JOIN mimiciv_hosp.prescriptions p ON p.hadm_id = a.hadm_id
WHERE p.drug IN ('Potassium Phosphate','Sodium Phosphate','Sodium Glycerophosphate')
  AND p.route IN ('IV','IV BOLUS')
  AND p.starttime IS NOT NULL;
CREATE INDEX idx_pcirx_stay ON public.pure_control_iv_rx(stay_id);

-- ---------------------------------------------------------------------
-- 3. Oral/enteral phosphate events — prescriptions explicit whitelist
--    Neutra-Phos / Phospha 250 Neutral / Phospha / Phos / phosp / Phosphorus
--    Routes: PO, PO/NG, NG, ORAL, J TUBE
--    EXCLUDES Fleet Phospho-Soda (laxative), Caphosol (mouthwash),
--             PhosLo (calcium acetate binder) — those go to audit table.
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS public.pure_control_po_rx;
CREATE TABLE public.pure_control_po_rx AS
SELECT a.stay_id, p.starttime AS ev_time, p.drug, p.route
FROM public.pure_control_anchor a
JOIN mimiciv_hosp.prescriptions p ON p.hadm_id = a.hadm_id
WHERE p.drug IN ('Neutra-Phos','Phospha 250 Neutral','Phospha','Phos','phosp','Phosphorus')
  AND p.route IN ('PO','PO/NG','NG','ORAL','J TUBE')
  AND p.starttime IS NOT NULL;
CREATE INDEX idx_pcpo_stay ON public.pure_control_po_rx(stay_id);

-- ---------------------------------------------------------------------
-- 4. Audit table: ALL phosphate-like drugs NOT in supplement whitelist
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS public.pure_control_excluded_drugs;
CREATE TABLE public.pure_control_excluded_drugs AS
SELECT p.drug, p.route, count(*) AS n_orders, count(DISTINCT p.hadm_id) AS n_hadm
FROM mimiciv_hosp.prescriptions p
WHERE p.drug ILIKE '%phos%'
  AND p.drug NOT IN ('Potassium Phosphate','Sodium Phosphate','Sodium Glycerophosphate',
                     'Neutra-Phos','Phospha 250 Neutral','Phospha','Phos','phosp','Phosphorus')
GROUP BY p.drug, p.route
ORDER BY n_orders DESC;

-- ---------------------------------------------------------------------
-- 5. Per-stay exposure summaries across windows + ever-in-admission
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS public.pure_control_exposure;
CREATE TABLE public.pure_control_exposure AS
WITH iv_agg AS (  -- ACTUAL IV INFUSIONS only (inputevents K Phos + Na Phos)
  SELECT stay_id,
    min(ev_time) FILTER (WHERE ev_time >  t0 AND ev_time <= t0 + interval '12 hours') AS iv_first_w12,
    min(ev_time) FILTER (WHERE ev_time >  t0 AND ev_time <= t0 + interval '48 hours') AS iv_first_w48,
    min(ev_time) FILTER (WHERE ev_time >  t0 AND ev_time <= t0 + interval '72 hours') AS iv_first_w72,
    min(ev_time) AS iv_first_any,
    sum(mmol)   FILTER (WHERE ev_time > t0 AND ev_time <= t0 + interval '24 hours') AS iv_mmol_24h,
    sum(mmol)   FILTER (WHERE ev_time > t0 AND ev_time <= t0 + interval '48 hours') AS iv_mmol_48h
  FROM public.pure_control_iv_ie JOIN public.pure_control_anchor USING (stay_id)
  GROUP BY stay_id
),
rx_agg AS (  -- IV prescriptions (orders) WITHOUT matching infusion — ambiguity flag only
  SELECT stay_id,
    min(ev_time) FILTER (WHERE ev_time >  t0 AND ev_time <= t0 + interval '12 hours') AS rx_first_w12,
    min(ev_time) FILTER (WHERE ev_time >  t0 AND ev_time <= t0 + interval '48 hours') AS rx_first_w48,
    min(ev_time) FILTER (WHERE ev_time >  t0 AND ev_time <= t0 + interval '72 hours') AS rx_first_w72,
    min(ev_time) AS rx_first_any
  FROM public.pure_control_iv_rx JOIN public.pure_control_anchor USING (stay_id)
  GROUP BY stay_id
),
po_agg AS (
  SELECT po.stay_id,
    min(po.ev_time) FILTER (WHERE po.ev_time > a.t0 AND po.ev_time <= a.t0 + interval '12 hours') AS po_first_w12,
    min(po.ev_time) FILTER (WHERE po.ev_time > a.t0 AND po.ev_time <= a.t0 + interval '48 hours') AS po_first_w48,
    min(po.ev_time) FILTER (WHERE po.ev_time > a.t0 AND po.ev_time <= a.t0 + interval '72 hours') AS po_first_w72,
    min(po.ev_time) AS po_first_any,
    string_agg(DISTINCT po.drug, '; ') FILTER (WHERE po.ev_time > a.t0 AND po.ev_time <= a.t0 + interval '48 hours') AS po_drugs_w48
  FROM public.pure_control_po_rx po
  JOIN public.pure_control_anchor a ON po.stay_id = a.stay_id
  GROUP BY po.stay_id
),
d7u AS (  -- dim tables are trial-period level; take one row per stay
  SELECT DISTINCT ON (stay_id) stay_id, nutrition_support_type, refeeding_risk
  FROM public.baseline_dim7_nutrition ORDER BY stay_id
),
d5u AS (
  SELECT DISTINCT ON (stay_id) stay_id, rrt_at_t0
  FROM public.baseline_dim5_physio_organ ORDER BY stay_id
),
base AS (
  SELECT a.*,
         d7.nutrition_support_type, d7.refeeding_risk,
         d5.rrt_at_t0
  FROM public.pure_control_anchor a
  LEFT JOIN d7u d7 ON a.stay_id = d7.stay_id
  LEFT JOIN d5u d5 ON a.stay_id = d5.stay_id
)
SELECT b.stay_id, b.hadm_id, b.t0, b.icu_intime, b.icu_outtime, b.baseline_phosphate,
       b.nutrition_support_type, b.refeeding_risk, b.rrt_at_t0,
       i.iv_first_w12, i.iv_first_w48, i.iv_first_w72, i.iv_first_any,
       COALESCE(i.iv_mmol_24h,0) AS iv_mmol_24h,
       COALESCE(i.iv_mmol_48h,0) AS iv_mmol_48h,
       r.rx_first_w12, r.rx_first_w48, r.rx_first_w72, r.rx_first_any,
       p.po_first_w12, p.po_first_w48, p.po_first_w72, p.po_first_any, p.po_drugs_w48,
       (i.iv_first_w12 IS NOT NULL) AS iv_w12,
       (i.iv_first_w48 IS NOT NULL) AS iv_w48,
       (i.iv_first_w72 IS NOT NULL) AS iv_w72,
       (p.po_first_w12 IS NOT NULL) AS po_w12,
       (p.po_first_w48 IS NOT NULL) AS po_w48,
       (p.po_first_w72 IS NOT NULL) AS po_w72,
       (r.rx_first_w48 IS NOT NULL AND i.iv_first_w48 IS NULL) AS iv_rx_only_w48,
       (i.iv_first_any IS NOT NULL OR p.po_first_any IS NOT NULL OR r.rx_first_any IS NOT NULL) AS any_supplement_ever,
       CASE
         WHEN i.iv_first_w48 IS NOT NULL THEN 'A_IV'
         WHEN p.po_first_w48 IS NOT NULL THEN 'B_PO_only'
         WHEN r.rx_first_w48 IS NOT NULL THEN 'C_rx_only_uninfused'
         WHEN i.iv_first_any IS NOT NULL OR p.po_first_any IS NOT NULL OR r.rx_first_any IS NOT NULL THEN 'C_outside_window'
         ELSE 'D_none'
       END AS exposure_w48,
       CASE
         WHEN i.iv_first_w12 IS NOT NULL THEN 'A_IV'
         WHEN p.po_first_w12 IS NOT NULL THEN 'B_PO_only'
         WHEN r.rx_first_w12 IS NOT NULL THEN 'C_rx_only_uninfused'
         WHEN i.iv_first_any IS NOT NULL OR p.po_first_any IS NOT NULL OR r.rx_first_any IS NOT NULL THEN 'C_outside_window'
         ELSE 'D_none'
       END AS exposure_w12,
       CASE
         WHEN i.iv_first_w72 IS NOT NULL THEN 'A_IV'
         WHEN p.po_first_w72 IS NOT NULL THEN 'B_PO_only'
         WHEN r.rx_first_w72 IS NOT NULL THEN 'C_rx_only_uninfused'
         WHEN i.iv_first_any IS NOT NULL OR p.po_first_any IS NOT NULL OR r.rx_first_any IS NOT NULL THEN 'C_outside_window'
         ELSE 'D_none'
       END AS exposure_w72
FROM base b
LEFT JOIN iv_agg i ON b.stay_id = i.stay_id
LEFT JOIN rx_agg r ON b.stay_id = r.stay_id
LEFT JOIN po_agg p ON b.stay_id = p.stay_id;

CREATE INDEX idx_pce_stay ON public.pure_control_exposure(stay_id);

-- ---------------------------------------------------------------------
-- 6. QC output
-- ---------------------------------------------------------------------
\echo '===== Row count (must = 22,685) ====='
SELECT count(*) AS n_stays FROM public.pure_control_exposure;

\echo '===== W48 four-class distribution ====='
SELECT exposure_w48, count(*) AS n, round(100.0*count(*)/sum(count(*)) OVER (),1) AS pct
FROM public.pure_control_exposure GROUP BY 1 ORDER BY 1;

\echo '===== W12 / W72 distributions ====='
SELECT win, cls, n FROM (
  SELECT 'W12' win, exposure_w12 cls, count(*) n FROM public.pure_control_exposure GROUP BY exposure_w12
  UNION ALL
  SELECT 'W48', exposure_w48, count(*) FROM public.pure_control_exposure GROUP BY exposure_w48
  UNION ALL
  SELECT 'W72', exposure_w72, count(*) FROM public.pure_control_exposure GROUP BY exposure_w72
) q ORDER BY win, cls;

\echo '===== A_IV vs study_cohort itt_group_48h cross-check ====='
WITH sc AS (
  SELECT stay_id, max(itt_group_48h) AS grp48
  FROM public.study_cohort GROUP BY stay_id
)
SELECT e.exposure_w48, sc.grp48, count(*) n
FROM public.pure_control_exposure e JOIN sc USING (stay_id)
GROUP BY 1,2 ORDER BY 1,2;

\echo '===== Audit table row count + top exclusions ====='
SELECT count(*) AS n_distinct_excluded_drug_routes FROM public.pure_control_excluded_drugs;
SELECT drug, route, n_orders FROM public.pure_control_excluded_drugs ORDER BY n_orders DESC LIMIT 12;
