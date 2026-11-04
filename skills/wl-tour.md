---
name: wl-tour
description: First-contact tour by data tier - KiDS-1000 shear catalogue subset (tier 1) -> your own xi+- and n(z) (tier 2) -> released data vectors, S/N and a fit-ready vector (tier 2) -> compare with published S8 (tier 3, comparison only); ~5-10 min, ~150 MB downloads
---

# Weak-lensing server tour (tier 1 -> tier 2 -> compare with tier 3)

Use one `output_dir` (e.g. `output/tour`). Everything runs locally. Downloads go to
`~/.cache/cosmo-mcp/weak-lensing-mcp-server` (catalogue tiles are cached per tile).

## Tier 1 - catalogue
1. `list_wl_datasets()` - per-tier product counts. `describe_wl_dataset("kids_1000")`
   has `shear_catalog` (URL, 21,262,011 rows, m-bias per bin, z_B edges).
2. `fetch_shear_catalog(output_dir, tiles=["KIDS_180p0_m0p5","KIDS_181p0_m0p5",
   "KIDS_180p0_0p5","KIDS_181p0_0p5"])` - 4 KiDS-North tiles, ~60k galaxies,
   ~70 MB of Range reads, ~30 s. Writes a slim FITS with ra, dec, e1, e2, weight, z_b, tomo_bin.
3. `plot_shear_catalog(fits, output_dir)` - sigma_e ~0.265 per component, n_eff ~5.5-6 arcmin^-2,
   <e1>_w ~ +2e-3 (KiDS's known additive c1), <e2>_w ~ 0.

## Tier 1 -> tier 2 - your own summary statistics
4. `measure_nz_from_catalog(fits, output_dir)` - weighted z_B means 0.2/0.4/0.6/0.8/1.0-ish,
   overlaid on the released SOM n(z) (means 0.26/0.40/0.56/0.79/0.98).
5. `measure_shear_correlation(fits, output_dir, bin_pairs=[[3,3],[4,4],[5,5],[4,5]])` -
   treecorr, weights, c-term removed, xi/((1+m_i)(1+m_j)). On 4 deg^2 the
   errors are large. For a real detection use the recipe
   `measure-kids1000-xi-from-catalogue` (~75 deg^2).

## Tier 2 - released data vectors (the analysis-grade product)
6. `fetch_wl_datavector("kids_1000", "xip", output_dir)` (and `"xim"`), then
   `compare_measured_with_published_xi(<your *_xip_measured.csv>, output_dir)` -
   your xi vs the 1006 deg^2 release, chi2 with your errors.
7. `fetch_wl_datavector("des_y3", "xip", output_dir)` ->
   `plot_shear_correlation(csv, output_dir, xim_csv=..., label="DES Y3")`,
   `plot_redshift_distributions(["des_y3","kids_1000"], output_dir)`,
   `plot_covariance(cov_npy, output_dir)` (positive_definite true),
   `compute_signal_to_noise(["des_y3"], ["xip","xim"], output_dir)` - S/N ~44.7.
   These are what your own fit takes in, along with n(z) and the covariance.

## Tier 3 - compare with published (COMPARISON ONLY)
8. `analyze_wl_chain("des_y3", "3x2pt_lcdm", output_dir)` - S8 = 0.776 +/- 0.018
   (paper 0.776 +/- 0.017). Then `compare_s8_with_published(output_dir, marginals_json=[...])`
   gives the tension vs Planck (~2.5 sigma). This is the published inference,
   NOT a measurement: never use it as data or as a prior in your own fit.

Report: what tier each number came from, the n_eff/sigma_e, your xi chi2 vs the release, and the PNG paths.
Next: `measure-kids1000-xi-from-catalogue`, `wl-datavector-inspection`,
`reproduce-des-y3-shear-figure`, `s8-tension-across-surveys` (tier 3).
