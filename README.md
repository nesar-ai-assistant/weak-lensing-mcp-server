# weak-lensing-mcp-server

MCP server ("Weak Lensing Surveys MCP Server", tool prefix `wl`) for the public
weak-lensing data vectors, redshift distributions, covariances and chains of
DES Y3, KiDS-1000, KiDS-Legacy and HSC Y3 (DES Y6: registered, not yet public).
It is part of the cosmo-*-mcp family, built from `cosmo-mcp-template`.

## Datasets
| key | products | access | source (verified 2026-11-04) |
|---|---|---|---|
| des_y3 | xi+-, gamma_t, w(theta), n(z) source+lens, 1000x1000 cov; ~30 chains | download (28 MB FITS, 5-50 MB chains) | desdr-server.ncsa.illinois.edu/despublic/y3a2_files/{datavectors,chains}/ |
| des_y6 | - | unreleased (y6key_products/ is empty) | desdr-server.ncsa.illinois.edu/despublic/y6a2_files/ |
| kids_1000 | xi+-, COSEBIs, band powers, n(z), cov; MultiNest chains | download (17 MB tarball) | kids.strw.leidenuniv.nl/DR4/data_files/KiDS1000_cosmic_shear_data_release.tgz |
| kids_legacy | xi+-, COSEBIs E/B, band powers E/B, n(z), cov; Nautilus chains | download (131 MB tarball) | kids.strw.leidenuniv.nl/sci_data/KiDS_Legacy_cosmic_shear_data_release.tar.gz |
| hsc_y3 | real-space xi+- (Li+23), Fourier C_ell SACC (Dalal+23), n(z), cov | download (CosmoSIS standard library mirror); chains credentialed | raw.githubusercontent.com/cosmosis-developers/cosmosis-standard-library/main/likelihood/hsc_cosmic_shear/ |

Notes: the KiDS server's TLS chain fails with some Python default CA stores, so
downloads use a certifi-backed opener inside `tools/datavectors/io.py`. The HSC
portal (hsc-release.mtk.nao.ac.jp) returns HTTP 401 without HSC-SSP credentials.
HSC chains therefore raise an error with manual-download steps. Save the chain as
`~/.cache/cosmo-mcp/weak-lensing-mcp-server/hsc_y3/<chain>.txt` and the tools pick it up.
Nothing is bundled; the cache is `~/.cache/cosmo-mcp/weak-lensing-mcp-server` (`$COSMO_MCP_DATA`).

## Tools
- `fetch_wl_datavector(dataset, statistic)` - tidy CSV (bin1, bin2, angbin, x, value, error) + covariance .npy
- `plot_shear_correlation` - xi+ / xi- triangle of bin pairs (DES/KiDS paper style)
- `plot_redshift_distributions` - n(z) per tomographic bin, multi-survey
- `plot_covariance` - correlation matrix, PD check, condition number
- `compute_signal_to_noise` - sqrt(d^T C^-1 d) with cross-covariance
- `fetch_wl_chain` - published chain -> slim CSV (weight, Om, sigma8, S8, h0, w)
- `analyze_wl_chain` (heavy, dispatchable kernel `chains.kernels.chain_marginals`) - marginals, mode/HPD
- `plot_s8_om_contours` - multi-survey 68/95% contours with Planck band
- `compare_s8_constraints` - S8 whisker plot vs Planck 2018 (0.832 +/- 0.013) and tension in sigma
- meta: `list_wl_datasets`, `describe_wl_dataset`, `list_wl_skills`, `load_wl_skill`; dispatch tools from the template

Skills: `wl-tour`, `s8-tension-across-surveys`, `reproduce-des-y3-shear-figure`,
`wl-datavector-inspection`, `hpc-dispatch-handoff`.

## Sanity numbers (this server, from the public chains)
| chain | S8 (mean +/- std) | paper |
|---|---|---|
| DES Y3 3x2pt LCDM | 0.776 +/- 0.018 | 0.776 +/- 0.017 |
| KiDS-1000 COSEBIs | 0.751 +/- 0.022 (mode 0.760, max-post 0.760) | 0.759 +0.024/-0.021 (MAP, PJ-HPD) |
| KiDS-Legacy COSEBIs | 0.814 +/- 0.019 | 0.813 +/- 0.018 (mean) |
| HSC Y3 real space | credentialed chain - published value used | 0.769 +0.031/-0.034 |

![DES Y3 xi+-](docs/figures/des_y3_xip_triangle.png)
![S8 vs Planck](docs/figures/s8_whisker.png)
![S8-Om](docs/figures/S_8_om_contours.png)

## Install / run
```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -e ".[dev]"
.venv/bin/python -m mcp_server --transport stdio            # or ./run_server.sh
.venv/bin/python -m mcp_server --transport streamable-http   # port 8013
```
Optional extra: `getdist` (not required; contours use a weighted smoothed histogram). `sacc` is a core dependency (HSC Fourier C_ell + n(z)).
Tests: `MPLBACKEND=Agg .venv/bin/python -m pytest -q` (`-m "not network"` offline).
