---
name: wl-tour
description: First-contact tour - one call per tool family (data vectors, n(z), covariance, S/N, chains, S8 comparison) on DES Y3, all local, ~1-3 min incl. ~60 MB of downloads
---

# Weak-lensing server tour

Use one `output_dir` (e.g. `output/tour`). Everything runs locally; first
use downloads the DES Y3 2pt FITS (28 MB) and one chain (~10 MB) into
`~/.cache/cosmo-mcp/weak-lensing-mcp-server`.

1. `list_wl_datasets()` - five keys: des_y3, des_y6 (unreleased), kids_1000,
   kids_legacy, hsc_y3. `describe_wl_dataset("des_y3")` lists statistics and chain keys.
2. `fetch_wl_datavector("des_y3", "xip", output_dir)` and the same with `"xim"`.
   Expect 200 points each, 10 bin pairs (4 source bins x 20 angles, centres 2.8-225 arcmin).
   Files: `des_y3_xip.csv` (bin1, bin2, angbin, x, value, error) + `des_y3_xip_cov.npy`.
3. `plot_shear_correlation(xip_csv, output_dir, xim_csv=xim_csv, label="DES Y3")` -
   xi+ in the upper triangle, xi- in the lower, theta*xi*1e4 on the y axis.
4. `plot_redshift_distributions(["des_y3"], output_dir)` - mean z of the
   4 source bins 0.33, 0.52, 0.74, 0.93 (`metadata.mean_z`).
5. `plot_covariance(xip_cov_npy, output_dir)` - `positive_definite` true.
6. `compute_signal_to_noise(["des_y3"], ["xip","xim"], output_dir)` -
   S/N about 44.7 (uncut; 3x2pt with gammat+wtheta about 264).
7. `fetch_wl_chain("des_y3", "3x2pt_lcdm", output_dir)` then
   `analyze_wl_chain("des_y3", "3x2pt_lcdm", output_dir)` -
   S8 = 0.776 +/- 0.018, Om = 0.339 +/- 0.031, sigma8 = 0.733 +/- 0.046
   (paper: S8 = 0.776 +/- 0.017).
8. `compare_s8_constraints(output_dir, marginals_json=[<step 7 json>])` -
   whisker plot vs Planck 0.832 +/- 0.013; DES Y3 3x2pt tension about 2.5-2.6 sigma.

Report: the S/N, the S8 value, the tension, and the PNG paths. Next steps:
`s8-tension-across-surveys`, `reproduce-des-y3-shear-figure`, `wl-datavector-inspection`.
