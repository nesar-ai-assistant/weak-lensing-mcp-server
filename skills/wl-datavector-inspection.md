---
name: wl-datavector-inspection
description: Inspect a weak-lensing data vector before using it - statistics available, n(z) per bin, covariance health (PD, condition number, correlations) and detection S/N, compared across DES Y3, KiDS and HSC Y3
---

# Data-vector inspection (n(z), covariance, S/N)

One `output_dir`. Datasets: des_y3, kids_1000, kids_legacy, hsc_y3 (des_y6 not public).

## 1. What is there
`describe_wl_dataset(name)` -> `statistics` (DES: xip, xim, gammat, wtheta;
KiDS: xip, xim, cosebis, bandpowers (+ cosebis_b, bandpowers_b B-modes for Legacy);
HSC: xip, xim (real space, Li+23), cl_ee (Fourier, Dalal+23)).

## 2. n(z)
`plot_redshift_distributions(["des_y3","kids_1000","kids_legacy","hsc_y3"], output_dir)`.
Mean z per bin (`metadata.mean_z`), roughly:
DES Y3 0.33/0.52/0.74/0.93; KiDS-1000 0.26/0.40/0.56/0.79/0.98; KiDS-Legacy 0.34/0.48/0.59/0.79/0.94/1.22;
HSC Y3 0.45/0.77/1.08/1.33. `sample="lens"` only for des_y3 (6 MagLim bins, 0.30-0.97).

## 3. Covariance health
`fetch_wl_datavector(ds, stat, output_dir)` -> `<ds>_<stat>_cov.npy`, then
`plot_covariance(cov_npy, output_dir)`. Require `positive_definite` true.
Measured: DES Y3 xi+ condition number 2.4e4, max off-diagonal correlation 0.75,
mean |r| 0.10; KiDS-1000 COSEBIs condition number 2.1e4, max |r| 0.92 between
neighbouring modes but mean |r| only 0.03 (most mode pairs nearly independent).

## 4. S/N
`compute_signal_to_noise([...], ["xip","xim"], output_dir)` - joint with cross-covariance.
Values measured here (uncut, full published vectors):
DES Y3 xi+- 44.7; KiDS-1000 xi+- 25.5; KiDS-Legacy xi+- 34.8; HSC Y3 xi+- 34.6 (CosmoSIS-library file, wider than the Li+23 cuts);
DES Y3 3x2pt (xip+xim+gammat+wtheta) 264. `snr_diag_only` ignores correlations and is
~2x too high - quote `snr`.
B-mode null test (KiDS-Legacy): `compute_signal_to_noise(["kids_legacy"], ["cosebis_b"], ...)`;
chi2 = snr^2 vs n_points. Measured on the FULL uncut vectors: COSEBIs B_n chi2 = 459 / 420
(p = 0.09, consistent with zero); band-power B chi2 = 233 / 168 (p = 7e-4) - the paper's
null tests use the analysis scale/mode range, so report this as "uncut, not the paper test".

## Report
Statistics, binning, mean z per bin, covariance PD + condition number, S/N with and
without correlations, plus caveats: DES vectors are uncut (S/N overstated relative to
the analysis), the HSC file is wider than the Li+23 cuts, KiDS covariances are analytic.
