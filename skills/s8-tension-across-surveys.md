---
name: s8-tension-across-surveys
description: Published results (for comparison only) - quantify the S8 tension of DES Y6, DES Y3, KiDS-1000, KiDS-Legacy and HSC Y3 against Planck 2018 from the surveys' published chains; contours and a whisker plot with sigma tensions
---

# S8 tension across surveys (published results, for comparison only)

Every number here is the surveys' PUBLISHED inference. Use it to put
YOUR own analysis of catalogues / summary statistics in context, e.g. after `measure-kids1000-xi-from-catalogue`
or a fit to `fetch_wl_datavector` output. Never treat it as a measurement or
as a prior in your own fit.

Goal: a defensible, sourced statement like "DES Y3 3x2pt sits 2.5 sigma below
Planck while KiDS-Legacy is consistent at 0.8 sigma". One `output_dir`.

## 1. Re-analyse the public chains
DES Y6 first: `fetch_des_y6_chain(output_dir, "3x2pt_lcdm")` (370 MB, ~15 s parse) gives
S8 0.789 +/- 0.011, Om 0.325 +/- 0.023, sigma8 0.759 +/- 0.033 (paper 0.789 +/- 0.012, Om 0.333 +0.023/-0.028,
arXiv:2601.14559) and writes `des_y6_3x2pt_lcdm_chain.csv` + `_marginals.json` for steps 2-3. Also
`1x2pt_lcdm` (shear only, S8 0.783 +/- 0.017) and `3x2pt_cmb_lcdm` (S8 0.811 +/- 0.008).

For each (dataset, chain) in (des_y3, 3x2pt_lcdm), (kids_1000, cosebis),
(kids_legacy, cosebis):
`fetch_wl_chain(dataset, chain, output_dir)` then `analyze_wl_chain(dataset, chain, output_dir)`.
KiDS downloads its release tarball (17 MB / 131 MB) once. `analyze_wl_chain`
is heavy: with `set_dispatch` it runs on the facility (see hpc-dispatch-handoff).

Expected (mean +/- std | mode, HPD68):
- DES Y3 3x2pt: S8 0.776 +/- 0.018 (paper 0.776 +/- 0.017, mean convention).
- KiDS-1000 COSEBIs: mean 0.751 +/- 0.022, mode 0.760, max-post sample 0.760
  (paper 0.759 +0.024/-0.021 is MAP + PJ-HPD - quote the mode/MAP to compare).
- KiDS-Legacy COSEBIs: 0.814 +/- 0.019 (paper mean 0.813 +/- 0.018; mode 0.815 +0.016/-0.021).
- HSC Y3: chains are credentialed (tool raises with manual-download steps);
  use the published 0.769 +0.031/-0.034 (real space, Li+23) and 0.776 +0.032/-0.033
  (Fourier, Dalal+23), which compare_s8_with_published includes by default.

## 2. Contours
`plot_s8_om_contours([<4 chain CSVs>], output_dir, labels=["DES Y6 3x2pt","DES Y3 3x2pt","KiDS-1000","KiDS-Legacy"])`.
Check: banana degeneracy in sigma8-Om collapses to a horizontal band in S8-Om;
KiDS-1000 prefers low Om (~0.26-0.29), DES ~0.34.

## 3. Tension
`compare_s8_with_published(output_dir, marginals_json=[<4 JSONs>])` (keeps the
published entries). Expected sigma vs Planck 0.832 +/- 0.013:
DES Y6 3x2pt ~2.4 (1D Gaussian; paper quotes 1.8 sigma from CMB with its full metric), DES Y3 3x2pt ~2.5-2.6, KiDS-1000 ~2.7 (published) / ~3.1 (chain mean), KiDS-Legacy ~0.8,
HSC Y3 real ~1.9, Fourier ~1.6.

## Report
Table of S8 per survey with convention (mean/std vs MAP/HPD), sigma tension,
and figure paths. Caveats: 1D Gaussian difference-in-means (papers use
full-posterior metrics, differences of a few 0.1 sigma); the Planck value is
TT,TE,EE+lowE+lensing; DES Y3 and Y6 overlap heavily (not independent of each other), DES Y3 and KiDS-1000 share sky (not independent);
the KiDS-Legacy shift is mostly redshift calibration, not a statistical fluke.
