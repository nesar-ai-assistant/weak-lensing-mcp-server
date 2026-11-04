---
name: hpc-dispatch-handoff
description: Run this server's heavy kernel (chains.kernels.chain_marginals - download + summarise a DES/KiDS chain) on ALCF Polaris or NERSC Perlmutter through the hep-genesis IRI engine - client-side pack export (hosted) or server-side set_dispatch (own machine)
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
