---
name: s8-tension-across-surveys
description: Quantify the S8 tension of DES Y3, KiDS-1000, KiDS-Legacy and HSC Y3 against Planck 2018 - re-analyse public chains, overlay Om-S8 contours, whisker plot with sigma tensions
---

# S8 tension across surveys

Goal: a defensible, sourced statement like "DES Y3 3x2pt sits 2.5 sigma below
Planck while KiDS-Legacy is consistent at 0.8 sigma". One `output_dir`.

## 1. Re-analyse the public chains
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
  (Fourier, Dalal+23), which compare_s8_constraints includes by default.

## 2. Contours
`plot_s8_om_contours([<3 chain CSVs>], output_dir, labels=["DES Y3 3x2pt","KiDS-1000","KiDS-Legacy"])`.
Check: banana degeneracy in sigma8-Om collapses to a horizontal band in S8-Om;
KiDS-1000 prefers low Om (~0.26-0.29), DES ~0.34.

## 3. Tension
`compare_s8_constraints(output_dir, marginals_json=[<3 JSONs>])` (keeps the
published entries). Expected sigma vs Planck 0.832 +/- 0.013:
DES Y3 3x2pt ~2.5-2.6, KiDS-1000 ~2.7 (published) / ~3.1 (chain mean), KiDS-Legacy ~0.8,
HSC Y3 real ~1.9, Fourier ~1.6.

## Report
Table of S8 per survey with convention (mean/std vs MAP/HPD), sigma tension,
and figure paths. Caveats: 1D Gaussian difference-in-means (papers use
full-posterior metrics, differences of a few 0.1 sigma); the Planck value is
TT,TE,EE+lowE+lensing; DES Y3 and KiDS-1000 share sky (not independent);
the KiDS-Legacy shift is mostly redshift calibration, not a statistical fluke.
