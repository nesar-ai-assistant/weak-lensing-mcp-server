---
name: measure-kids1000-xi-from-catalogue
description: Tier 1 -> tier 2 - fetch a ~75 deg^2 KiDS-1000 gold shear-catalogue subset, QA it, measure xi+/xi- in 5 tomographic bins with treecorr (weights, m-bias, c-term, jackknife), histogram n(z), and validate against the released KiDS-1000 xi+- with the area-scaled covariance; ~10 min, ~800 MB of Range reads
---

# Measure KiDS-1000 cosmic shear from the catalogue (tier 1 -> tier 2)

Goal: your own xi+- data vector from the per-galaxy shear catalogue,
validated against the released KiDS-1000 xi+-. The release covers 1006 deg^2
(777.4 deg^2 effective, Giblin+21). A subset reproduces its SHAPE and
AMPLITUDE within errors that are about sqrt(777/A_sub) larger. It does not
reproduce it point by point.

1. `fetch_shear_catalog(output_dir, ra_range=[170, 190], dec_range=[-2, 2], max_tiles=200)`
   - KiDS-North block of 79 tiles, 1,714,314 galaxies, ~75 deg^2 occupied (~64 deg^2 unmasked),
     per-bin counts 150k/302k/494k/364k/403k. Cost ~790 MB of Range reads,
     ~6 min the first time, then seconds (tiles cached in `kids_1000/shearcat_tiles/`).
   - Smaller test: `tiles=[...4 names...]` or `max_rows=50000`. Bin cut: `tomo_bins=[3,4,5]`.
2. `plot_shear_catalog(fits, output_dir)` - expect sigma_e = 0.265 per component,
   n_eff = 5.1 arcmin^-2 (occupied-area estimate; Giblin+21 quote 6.2 over the unmasked area),
   c1 = +4e-4, c2 = -1e-4 overall. Per bin, |c| <~ 1e-3.
3. `measure_nz_from_catalog(fits, output_dir)` - weighted z_B means
   0.22/0.41/0.61/0.80/1.02 vs the released SOM n(z) 0.26/0.40/0.56/0.79/0.98. The
   z_B histogram is confined to the bin edges and lacks the calibrated outlier tails,
   so use the RELEASE n(z) for theory.
4. `measure_shear_correlation(fits, output_dir, npatch=40)` - all 15 pairs, 9 log bins
   over 0.5-300 arcmin (the release binning). Takes 243 s locally on a 4-core laptop (heavy; with
   `set_dispatch` it runs as `catalogs.kernels.shear_correlation` on the facility).
   - Per bin it removes the weighted mean e1/e2 and divides by (1+m_i)(1+m_j) with
     m = -0.009, -0.011, -0.015, +0.002, +0.007 (Asgari+21 Table 1).
   - Sanity: xi+ ~ 1e-5 at 10 arcmin for bins 4-5 (release: 1.8e-5 at 12 arcmin for 4-4/4-5,
     2.5e-5 for 5-5), so theta*xi+ ~ 1-3e-4 arcmin. Jackknife errors on 64 deg^2 are ~1e-5 at
     13 arcmin, giving S/N of order 1 per point. Measured 4-5 pair: 1.2e-5 +/- 0.8e-5 at 13.5 arcmin.
   - All pairs, theta > 2': chi2(xi+ vs 0) = 124.5/105 with diagonal jackknife errors (correlated
     bins, so this is only indicative).
5. `fetch_wl_datavector("kids_1000", "xip", output_dir)` (the release), then
   `compare_measured_with_published_xi(<*_xip_measured.csv>, output_dir, subset_area_deg2=63.9)`:
   - With the release covariance (noise + cosmic variance) x 777.4/63.9: A = 0.86 +/- 0.18,
     chi2(measured - release) = 82.2/105 (xi+). xi-: A = 0.86 +/- 0.19, chi2 = 114.1/105.
     The subset therefore agrees with the release and detects the signal at ~5 sigma.
   - With diagonal jackknife errors only: A = 0.67 +/- 0.10 (too optimistic, since it ignores bin
     correlations and cosmic variance).

Pitfalls
- e1/e2 in the slim FITS are RAW. Always use measure_shear_correlation (or remove c
  and correct m yourself) before you interpret anything.
- c-term removal also removes any true constant shear on the patch. On ~1 deg^2 that
  matters; on >= 10 deg^2 it is negligible.
- Small scales (theta < 2') in a subset are noisy and the release itself cuts them for
  cosmology, so start amplitude fits at 2-5 arcmin.
- Contiguous blocks matter: scattered tiles lose all the large-scale pairs.
- KiDS-North RA 128-238 / Dec -4..+4; KiDS-South RA 330-55 / Dec -36..-26.

Final step (tier 3, comparison only): to put your measurement in cosmological context,
compare with the published KiDS-1000 S8 via `compare_s8_with_published` / `analyze_wl_chain("kids_1000", "xipm")`.
Do NOT use those chains as a prior when fitting your xi.
