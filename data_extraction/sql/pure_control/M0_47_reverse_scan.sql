-- =====================================================================
-- M0_47: D_none zero-contamination REVERSE SCAN audit (review m5 artifact)
-- Post-hoc pure-control analysis, read-only. Every check in section A
-- must return 0 rows; section B lists remaining '%phos%' prescription
-- names in D_none for manual audit (black-list / non-supplement only).
-- Run log: output/pure_control/qc/M0_47_reverse_scan.txt
-- =====================================================================
\echo '=== A1. D_none stays with actual IV phosphate infusion (225834/225835, any time, 0.006-500 mmol) — expect 0 ==='
SELECT count(*) AS iv_rows, count(DISTINCT ie.stay_id) AS iv_stays
FROM mimiciv_icu.inputevents ie
JOIN public.pure_control_cohort c ON c.stay_id = ie.stay_id
WHERE c.exposure_w48 = 'D_none'
  AND ie.itemid IN (225834,225835)
  AND ie.amount BETWEEN 0.006 AND 500;

\echo '=== A2. D_none stays with white-list IV phosphate PRESCRIPTION (any time) — expect 0 ==='
SELECT count(*) AS rx_rows, count(DISTINCT p.hadm_id) AS rx_hadm
FROM mimiciv_hosp.prescriptions p
JOIN public.pure_control_cohort c ON c.hadm_id = p.hadm_id
WHERE c.exposure_w48 = 'D_none'
  AND p.drug IN ('Potassium Phosphate','Sodium Phosphate','Sodium Glycerophosphate')
  AND p.route IN ('IV','IV BOLUS');

\echo '=== A3. D_none stays with white-list ORAL/ENTERAL phosphate PRESCRIPTION (any time) — expect 0 ==='
SELECT count(*) AS po_rows, count(DISTINCT p.hadm_id) AS po_hadm
FROM mimiciv_hosp.prescriptions p
JOIN public.pure_control_cohort c ON c.hadm_id = p.hadm_id
WHERE c.exposure_w48 = 'D_none'
  AND p.drug IN ('Neutra-Phos','Phospha 250 Neutral','Phospha','Phos','phosp','Phosphorus')
  AND p.route IN ('PO','PO/NG','NG','ORAL','J TUBE');

\echo '=== A4. D_none rows under alternative phosphate-like input itemids (220448,225925,225612,225677) — expect 0 ==='
SELECT ie.itemid, count(*) AS n_rows, count(DISTINCT ie.stay_id) AS n_stays
FROM mimiciv_icu.inputevents ie
JOIN public.pure_control_cohort c ON c.stay_id = ie.stay_id
WHERE c.exposure_w48 = 'D_none'
  AND ie.itemid IN (220448,225925,225612,225677)
GROUP BY ie.itemid ORDER BY ie.itemid;

\echo '=== A5. Classification-internal flag any_supplement_ever among D_none — expect 0 ==='
SELECT count(*) AS flagged_ever
FROM public.pure_control_cohort
WHERE exposure_w48 = 'D_none' AND any_supplement_ever IS TRUE;

\echo '=== B. All %phos% prescriptions present in D_none hospitalizations (manual audit; expect ONLY black-list/non-repletion products) ==='
SELECT p.drug, p.route, count(*) AS n_orders, count(DISTINCT p.hadm_id) AS n_hadm
FROM mimiciv_hosp.prescriptions p
JOIN public.pure_control_cohort c ON c.hadm_id = p.hadm_id
WHERE c.exposure_w48 = 'D_none' AND p.drug ILIKE '%phos%'
GROUP BY p.drug, p.route
ORDER BY n_orders DESC;
