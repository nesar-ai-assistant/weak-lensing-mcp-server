---
name: reproduce-des-y3-shear-figure
description: Reproduce the DES Y3 cosmic-shear xi+/xi- data figure (Amon+22 / Secco+22 Fig. 1 style), plus the KiDS-1000 / KiDS-Legacy / HSC Y3 equivalents
---

# Reproduce the DES Y3 xi+- figure

1. `fetch_wl_datavector("des_y3", "xip", output_dir)` and `("des_y3", "xim", output_dir)`.
   200 + 200 points; theta 2.5-250 arcmin bin edges (centres 2.8-225), 20 log bins; bins (i,j) with i<=j, i,j in 1..4.
2. `plot_shear_correlation(xip_csv, output_dir, xim_csv=xim_csv, label="DES Y3 cosmic shear")`.
   Upper triangle xi+ (blue), lower xi- (red), y = theta*xi*1e4 (the paper convention, in arcmin).
3. Sanity checks against the paper figure:
   - theta*xi+*1e4 rises with source redshift: max over theta is 0.6 (1,1), 1.0 (2,2), 2.4 (3,3), 3.3 (4,4).
   - theta*xi- is small at small theta and peaks at large theta (~180 arcmin for (4,4)).
   - errors grow at small theta in the low-z bins (shape noise) and at large theta (cosmic variance).
   - The published file has NO scale cuts; the paper greys out the small-scale points
     (roughly theta < 10-30 arcmin for xi-, depending on bin). Say so when presenting.
4. Optional: same plot for `kids_1000` (5 bins, 9 angles 0.5-300 arcmin), `kids_legacy`
   (6 bins, 9 angles) and `hsc_y3` (4 bins; the CosmoSIS-library real-space file holds
   xi+ at 3.5-120 arcmin (13 angles) and xi- at 15-216 arcmin (10 angles); Li+23 used
   xi+ 7.1-56.6 and xi- 31.2-248 arcmin, so apply those cuts yourself when comparing).
5. `plot_covariance(xip_cov_npy, output_dir)` to show the strong correlation between
   neighbouring large-scale angular bins.

Report: PNG path, number of points per statistic, and the scale-cut caveat.
