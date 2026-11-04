---
name: hpc-dispatch-handoff
description: Run this server's heavy kernels (catalogs.kernels.shear_correlation - treecorr xi+- from a shear catalogue; catalogs.kernels.kids1000_subset - KiDS-1000 catalogue subset; chains.kernels.chain_marginals - chain summary) on ALCF Polaris or NERSC Perlmutter through the hep-genesis IRI engine - client-side pack export (hosted) or server-side set_dispatch (own machine)
---

# HPC dispatch handoff

This server holds NO facility credentials, project, workdir or environment
settings. All facility access goes through IRI (job submission) and Globus /
the IRI filesystem API (staging) via the hep-genesis engine; nothing here
writes batch scripts or opens SSH sessions. Tool names `set_dispatch`,
`get_dispatch`, `auth_status`, `export_dispatch_pack` are a cross-server
convention matched by the hep-genesis harness.

## Client-side (hosted servers - recommended)
1. `export_dispatch_pack` -> the client saves `files` locally (do NOT echo
   their contents) and reads `kernels`: each has `function`, `args_shape`,
   `pip_deps`, `duration_hint_s`, `returns`.
2. On the client's facility server call
   `run_pack_kernel(pack=<saved>/tools, function=<kernel function>,
   args={...}, pip_deps=<kernel pip_deps>, duration>=duration_hint_s,
   artifact_dir=<run dir>)`.
3. Long or parallel jobs: `wait=False`, then `find_dispatched_jobs` /
   `collect_dispatched_job(record=...)`.

## This server's kernels
- `catalogs.kernels.shear_correlation` (used by `measure_shear_correlation`, from shear catalogue to xi+-):
  args `{"catalog_file": <slim FITS from fetch_shear_catalog, staged into the job cwd>,
  "bin_pairs": [[1,1],[1,2],...], "m_bias": {"1": -0.009, ...}, "min_sep": 0.5, "max_sep": 300,
  "nbins": 9, "bin_slop": 0.05, "remove_c": true, "npatch": 0, "out_csv": "xipm_measured.csv"}`;
  pip_deps numpy, astropy, treecorr>=5.0.1; duration_hint_s 1800. Returns per-bin c1/c2/n_eff/sigma_e
  and rows (bin1, bin2, theta, xip, xim, errors, npairs). Scale: 1.7M galaxies / 75 deg^2,
  15 pairs, 40 jackknife patches took 243 s locally (4-core laptop). Full 1006 deg^2 belongs on a facility.
- `catalogs.kernels.kids1000_subset` (the reduction behind `fetch_shear_catalog`): args
  `{"tile_ranges": {tile: [[row_start, n_rows], ...]}, "out_file": "kids1000_subset.fits",
  "ra_range", "dec_range", "tomo_bins", "max_rows"}`; pip_deps numpy, astropy, certifi;
  duration_hint_s 3600. Get `tile_ranges` from `data/kids1000_gold_tile_runs.npz`. ~17 MB of
  remote Range reads per tile, which makes it useful for pulling hundreds of tiles straight onto facility storage.
- `chains.kernels.chain_marginals` (used by `analyze_wl_chain`): args
  `{"url": <DES chain .txt or KiDS release tarball URL>, "member": <path inside
  the tarball or null>, "params": null}`; pip_deps numpy, certifi; duration_hint_s 900.
  Returns S_8/omega_m/sigma_8/h0 marginals, n_eff, corr_S8_Om. Get url/member from
  `describe_wl_dataset(...)["chains"]` (DES: chains base URL + file name;
  KiDS: `archive` + member). HSC Y3 chains are credentialed and cannot be dispatched.

## Server-side (server runs on the user's own machine)
1. `auth_status` -> IRI + transfer tokens must be "available".
2. `set_dispatch("perlmutter" | "polaris")`. Heavy tools now run each call
   as one facility job; light tools stay local. `set_dispatch("local")` to undo.
3. Never blind-retry a failed dispatch: surface the error. A walltime
   rejection may be resubmitted ONCE with a shorter duration.
