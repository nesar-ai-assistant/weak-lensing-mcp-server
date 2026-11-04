"""Tier-1 shear catalogues: subset download, QA plots, and tier-1 -> tier-2 reductions (xi+-, n(z)).

KiDS-1000 gold WL catalogue (Giblin+21) is the public catalogue-level product: one 17.7 GB FITS
table read as SUBSETS with HTTP Range requests (bundled tile index: data/kids1000_gold_tile_runs.npz).
DES Y3 (312 GB HDF5, no subset service) and HSC Y3 (credentialed) catalogues are documented but not fetched.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Literal

import numpy as np
from pydantic import Field, validate_call

from ..common import BUNDLED_DATA, ArtifactResult, data_cache_dir, mpl, ok, param_slug, resolve_outdir, save_fig, tier
from ..registry import DATASETS
from . import kernels

__all__ = ["fetch_shear_catalog", "plot_shear_catalog", "measure_shear_correlation", "measure_nz_from_catalog",
           "compare_measured_with_published_xi"]

DISPATCH_KERNELS = {
    "shear_correlation": {
        "function": "catalogs.kernels.shear_correlation",
        "args_shape": "{'catalog_file': <slim FITS name>, 'subset': <kids1000_subset kwargs to rebuild it on the "
                      "node if not staged: tile_ranges, ra_range, dec_range, tomo_bins, max_rows>, "
                      "'bin_pairs': [[1,1],[1,2],...], 'm_bias': {'1': -0.009, ...}, 'min_sep': 0.5, "
                      "'max_sep': 300.0, 'nbins': 9, 'bin_slop': 0.05, 'remove_c': True, 'npatch': 0, "
                      "'out_csv': 'xipm_measured.csv'}",
        "pip_deps": ["numpy", "astropy", "treecorr>=5.0.1"],
        "duration_hint_s": 1800,
        "returns": "per-bin c1/c2/m/n/n_eff/sigma_e and rows (bin1, bin2, theta, xip, xim, errors, npairs); "
                   "writes out_csv in the job cwd",
    },
    "kids1000_subset": {
        "function": "catalogs.kernels.kids1000_subset",
        "args_shape": "{'tile_ranges': {tile: [[row_start, n_rows], ...]}, 'out_file': 'kids1000_subset.fits', "
                      "'ra_range': [lo, hi] | None, 'dec_range': [lo, hi] | None, 'tomo_bins': [1..5] | None, "
                      "'max_rows': int | None}",
        "pip_deps": ["numpy", "astropy", "certifi"],
        "duration_hint_s": 3600,
        "returns": "file, n_rows, tiles, bytes_downloaded, area_deg2, tomo_counts; the slim FITS stays in the job cwd",
    },
}

OutDir = Annotated[str, Field(min_length=1, description="Directory where output files are written.")]
CatalogDataset = Literal["kids_1000", "des_y3", "hsc_y3", "kids_legacy"]
SC = DATASETS["kids_1000"]["shear_catalog"]
KIDS1000_AREA_EFF = 777.4  # deg^2, effective KiDS-1000 shear area (Giblin et al. 2021)


def _not_public(dataset: str) -> None:
    if dataset != "kids_1000":
        note = DATASETS[dataset].get("tier1_note", "no public shear catalogue")
        raise RuntimeError(f"{dataset}: shear catalogue not available through this server. {note}")


def _tile_index() -> dict:
    z = np.load(BUNDLED_DATA / "kids1000_gold_tile_runs.npz")
    return {"tiles": [str(t) for t in z["tiles"]], "runs": z["runs"].astype(np.int64)}


def _tile_ranges(index: dict, tiles: list[str]) -> dict:
    tid = {t: i for i, t in enumerate(index["tiles"])}
    r = index["runs"]
    return {t: [[int(s), int(n)] for s, n in r[r[:, 0] == tid[t]][:, 1:]] for t in tiles}


@tier(1)
@validate_call
def fetch_shear_catalog(
    output_dir: OutDir,
    dataset: Annotated[CatalogDataset, Field(description="Survey. Only kids_1000 has a publicly downloadable "
                                                         "catalogue; the others raise with the reason.")] = "kids_1000",
    tiles: Annotated[list[str] | None, Field(description="Explicit KiDS tile names, e.g. ['KIDS_180p0_m0p5'] "
                                                         "(1x1 deg each; names encode the centre RA/Dec).")] = None,
    patch: Annotated[Literal["N", "S"] | None, Field(description="KiDS-North (equatorial, RA 128-238) or "
                                                                 "KiDS-South (Dec ~ -31, RA 330-55).")] = None,
    ra_range: Annotated[list[float] | None, Field(min_length=2, max_length=2, description="[ra_min, ra_max] deg; "
                                                  "ra_min > ra_max wraps through 0.")] = None,
    dec_range: Annotated[list[float] | None, Field(min_length=2, max_length=2,
                                                   description="[dec_min, dec_max] deg.")] = None,
    tomo_bins: Annotated[list[int] | None, Field(description="Keep only these KiDS-1000 tomographic bins (1-5: "
                                                             "z_B edges 0.1,0.3,0.5,0.7,0.9,1.2).")] = None,
    max_tiles: Annotated[int, Field(ge=1, le=1006, description="Cap on tiles read (~18 MB of the remote file "
                                                               "and ~10-40 s per tile).")] = 4,
    max_rows: Annotated[int | None, Field(ge=100, description="Stop once this many galaxies pass the cuts.")] = None,
) -> ArtifactResult:
    """Download a SUBSET of the KiDS-1000 gold weak-lensing shear catalogue (tier 1) as a slim FITS.

    Use this first for any catalogue-level / measure-it-yourself question.
    The 17.7 GB, 21.3M-galaxy file is never downloaded whole: whole tiles
    are read with HTTP Range requests using the bundled tile->row index,
    then cut to the RA/Dec box, tomographic bins and row limit. Columns:
    ra, dec, e1, e2 (raw lensfit ellipticities, NOT m-corrected), weight,
    z_b (BPZ photo-z) + z_b_min/max, mag_r, snr, psf_e1/2, tomo_bin, tile.
    Tiles are cached per tile in the data cache, so repeats are free.
    Pick a contiguous block (ra_range/dec_range) for correlation work.
    """
    _not_public(dataset)
    index = _tile_index()
    chosen = kernels.select_tiles(index["tiles"], tiles, ra_range, dec_range, patch, max_tiles)
    if not chosen:
        raise ValueError("Subset selects no KiDS-1000 tiles: KiDS-N covers RA 128-238, Dec -4..+4; KiDS-S "
                         "RA 330-55, Dec -36..-26.")
    out = resolve_outdir(output_dir)
    params = {"t": ",".join(chosen), "ra": ra_range, "dec": dec_range, "tb": tomo_bins, "mr": max_rows}
    fits_path = out / f"kids1000_shearcat_{len(chosen)}tiles_{param_slug(params)}.fits"
    res = kernels.kids1000_subset(_tile_ranges(index, chosen), str(fits_path), ra_range=ra_range,
                                  dec_range=dec_range, tomo_bins=tomo_bins, max_rows=max_rows,
                                  cache_dir=str(data_cache_dir() / "kids_1000" / "shearcat_tiles"), url=SC["url"])
    msg = (f"KiDS-1000 gold shear catalogue subset: {res['n_rows']:,} galaxies from {res['n_tiles']} tile(s), "
           f"~{res['area_deg2']:.2f} deg^2 occupied; per-bin counts {res['tomo_counts']}.")
    return ok(msg, [fits_path], dataset=dataset, n_rows=res["n_rows"], tiles=res["tiles"],
              area_deg2=res["area_deg2"], tomo_counts=res["tomo_counts"], columns=res["columns"],
              bytes_downloaded=res["bytes_downloaded"], source_url=SC["url"],
              note="e1/e2 are raw (no m-correction, c-term not removed); measure_shear_correlation applies both.",
              citation=SC["citation"])


fetch_shear_catalog.weight = "heavy"


def _load(path: str) -> dict:
    return kernels.load_catalog(path)


@tier(1)
@validate_call
def plot_shear_catalog(
    catalog_fits: Annotated[str, Field(description="Slim FITS from fetch_shear_catalog.")],
    output_dir: OutDir,
) -> ArtifactResult:
    """QA plots of a tier-1 shear catalogue: e1/e2 and weight histograms, z_B, and the sky footprint.

    Use right after fetch_shear_catalog to check what was downloaded.
    Reports weighted mean e1/e2 (the additive c-terms, ~1e-4..1e-3 for
    KiDS), the per-component weighted ellipticity dispersion sigma_e
    (~0.27 for KiDS-1000, Giblin+21), the effective number density
    n_eff = (sum w)^2 / sum w^2 / area (Heymans+12; ~6-7 arcmin^-2 for
    the KiDS-1000 gold sample) and per-bin counts.
    """
    plt = mpl()
    c = _load(catalog_fits)
    w, e1, e2 = c["weight"], c["e1"], c["e2"]
    area = kernels.footprint_area(c["ra"], c["dec"])
    neff = float(w.sum() ** 2 / (w ** 2).sum() / (area * 3600.0)) if area > 0 else float("nan")
    sig_e = float(np.sqrt(np.sum(w ** 2 * (e1 ** 2 + e2 ** 2)) / (2 * np.sum(w ** 2))))
    c1, c2 = float(np.average(e1, weights=w)), float(np.average(e2, weights=w))
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
    h = np.linspace(-1, 1, 81)
    ax[0, 0].hist(e1, h, weights=w, histtype="step", lw=1.4, label=f"e1 (<e1>_w = {c1:.1e})")
    ax[0, 0].hist(e2, h, weights=w, histtype="step", lw=1.4, label=f"e2 (<e2>_w = {c2:.1e})")
    ax[0, 0].set_xlabel("ellipticity")
    ax[0, 0].set_ylabel("weighted counts")
    ax[0, 0].legend(fontsize=8)
    ax[0, 1].hist(w, 60, color="C2")
    ax[0, 1].set_xlabel("lensfit weight")
    ax[0, 1].set_yscale("log")
    tb = c.get("tomo_bin", np.zeros(len(w)))
    zb_edges = np.arange(0.0, 2.05, 0.05)
    for b in sorted(set(np.unique(tb).astype(int)) - {0}):
        s = tb == b
        ax[1, 0].hist(c["z_b"][s], zb_edges, weights=w[s], histtype="step", lw=1.3, label=f"bin {b}")
    ax[1, 0].hist(c["z_b"], zb_edges, weights=w, histtype="step", color="k", lw=0.8, label="all")
    ax[1, 0].set_xlabel("z_B (BPZ)")
    ax[1, 0].legend(fontsize=7)
    n = len(w)
    sub = np.random.default_rng(0).choice(n, size=min(n, 60000), replace=False)
    ra = (c["ra"] - np.median(c["ra"]) + 180) % 360 - 180 + np.median(c["ra"])
    ax[1, 1].scatter(ra[sub], c["dec"][sub], s=0.2, c=w[sub], cmap="viridis", rasterized=True)
    ax[1, 1].set_xlabel("RA [deg]")
    ax[1, 1].set_ylabel("Dec [deg]")
    ax[1, 1].invert_xaxis()
    ax[1, 1].set_title(f"{n:,} galaxies, ~{area:.2f} deg$^2$", fontsize=9)
    fig.suptitle(f"{Path(catalog_fits).stem}: sigma_e = {sig_e:.3f}, n_eff = {neff:.2f} arcmin$^{{-2}}$")
    fig.tight_layout()
    out = resolve_outdir(output_dir)
    png = save_fig(fig, out, f"{Path(catalog_fits).stem}_qa")
    counts = {str(int(b)): int((tb == b).sum()) for b in np.unique(tb)}
    return ok(f"{n:,} galaxies; sigma_e={sig_e:.3f}, n_eff={neff:.2f}/arcmin^2, c1={c1:.1e}, c2={c2:.1e}.", [png],
              n=n, area_deg2=round(area, 3), sigma_e=round(sig_e, 4), n_eff_arcmin2=round(neff, 3),
              c1=c1, c2=c2, ra_range=[float(ra.min()), float(ra.max())],
              dec_range=[float(c["dec"].min()), float(c["dec"].max())], tomo_counts=counts)


def _subset_spec(cat: dict) -> dict:
    """kids1000_subset kwargs that rebuild exactly this catalogue on a facility node (no input staging):
    same tiles in the same order, the catalogue's own RA/Dec bounds, its bins, and its row count."""
    tiles = list(dict.fromkeys(str(t) for t in cat["tile"]))
    ra = np.asarray(cat["ra"], float)
    dec = np.asarray(cat["dec"], float)
    wraps = ra.max() - ra.min() > 180
    ra_range = [float(ra[ra > 180].min()), float(ra[ra < 180].max())] if wraps else [float(ra.min()), float(ra.max())]
    return {"tile_ranges": _tile_ranges(_tile_index(), tiles), "ra_range": ra_range,
            "dec_range": [float(dec.min()), float(dec.max())],
            "tomo_bins": sorted(int(b) for b in set(np.unique(cat["tomo_bin"]).tolist())), "max_rows": len(ra)}


def _bin_pairs(spec: list[list[int]] | None, bins: list[int]) -> list[list[int]]:
    if spec:
        return [sorted(map(int, p)) for p in spec]
    return [[i, j] for i in bins for j in bins if i <= j]


@tier(1)
@validate_call
def measure_shear_correlation(
    catalog_fits: Annotated[str, Field(description="Slim FITS from fetch_shear_catalog (tier-1 input).")],
    output_dir: OutDir,
    bin_pairs: Annotated[list[list[int]] | None, Field(description="Tomographic bin pairs, e.g. [[3,3],[4,5]] "
                                                                   "(default: all auto + cross pairs present).")] = None,
    min_sep: Annotated[float, Field(gt=0, description="Min separation [arcmin] (KiDS-1000 release: 0.5).")] = 0.5,
    max_sep: Annotated[float, Field(gt=0, description="Max separation [arcmin] (KiDS-1000 release: 300).")] = 300.0,
    nbins: Annotated[int, Field(ge=2, le=60, description="Log-spaced theta bins (KiDS-1000 release: 9).")] = 9,
    apply_m_bias: Annotated[bool, Field(description="Divide by (1+m_i)(1+m_j) with the KiDS-1000 per-bin m "
                                                    "(Asgari+21 Table 1).")] = True,
    remove_c: Annotated[bool, Field(description="Subtract the weighted mean e1, e2 per bin (additive c-term).")] = True,
    npatch: Annotated[int, Field(ge=0, le=200, description="Jackknife patches for errors (0 = shape-noise "
                                                           "only; needs >~ 4 deg^2 for >= 10 patches).")] = 0,
    bin_slop: Annotated[float, Field(ge=0, le=1, description="treecorr bin_slop (accuracy vs speed).")] = 0.05,
) -> ArtifactResult:
    """Measure xi+(theta) and xi-(theta) per tomographic bin pair from a shear catalogue (tier 1 -> tier 2).

    Use to turn fetch_shear_catalog output into your own cosmic-shear
    data vector. Weighted treecorr GG correlation with per-bin additive
    c-term removal and multiplicative-bias correction xi / ((1+m_i)(1+m_j)).
    Output CSV has the fetch_wl_datavector layout (bin1, bin2, angbin, x,
    value, error) for xip and xim, so plot_shear_correlation and
    compare_measured_with_published_xi read it directly. Errors are
    shape noise (or jackknife with npatch) and do not include cosmic
    variance; a few deg^2 subset gives a noisy, not unbiased-in-area,
    estimate vs the 1006 deg^2 release. Sanity: theta*xi+ ~ 1e-4 x
    (1-5) arcmin at 10 arcmin for z~0.5-1 pairs (Asgari+21 Fig. 1).
    Heavy: with set_dispatch it runs catalogs.kernels.shear_correlation on the facility.
    """
    import pandas as pd

    from mcp_server.dispatch import remote_site, run_kernel  # lazy: server-only

    cat = _load(catalog_fits)
    present = sorted(int(b) for b in set(np.unique(cat["tomo_bin"]).tolist()) - {0})
    pairs = _bin_pairs(bin_pairs, present)
    m = {str(b): (SC["m_bias"][b - 1] if apply_m_bias else 0.0) for b in range(1, 6)}
    args = dict(bin_pairs=pairs, m_bias=m, min_sep=min_sep, max_sep=max_sep, nbins=nbins, bin_slop=bin_slop,
                remove_c=remove_c, npatch=npatch)
    site = remote_site()
    if site:
        res = run_kernel("catalogs.kernels.shear_correlation",
                         {"catalog_file": Path(catalog_fits).name, **args, "out_csv": "xipm_measured.csv",
                          "subset": _subset_spec(cat)},
                         pip_deps=DISPATCH_KERNELS["shear_correlation"]["pip_deps"], duration=1800)
        payload, computed_on = res["result"], res.get("host", site)
        info, rows = {"bins": payload["bins"]}, payload["rows"]
    else:
        out_all = kernels.shear_correlation_arrays(cat, num_threads=None, **args)
        info, rows, computed_on = out_all[0], out_all[1:], "local"
    df = pd.DataFrame(rows)
    out = resolve_outdir(output_dir)
    stem = Path(catalog_fits).stem
    files = []
    for s in ("xip", "xim"):
        tidy = pd.DataFrame({"bin1": df.bin1, "bin2": df.bin2, "angbin": df.angbin, "x": df.theta_arcmin,
                             "value": df[s], "error": df[f"{s}_err"], "npairs": df.npairs})
        p = out / f"{stem}_{s}_measured.csv"
        tidy.to_csv(p, index=False)
        files.append(p)
    js = out / f"{stem}_xipm_measured.json"
    js.write_text(json.dumps({"bins": info["bins"], "settings": args, "computed_on": computed_on}, indent=2))
    files.append(js)
    near10 = {}
    for (i, j), g in df.groupby(["bin1", "bin2"]):
        k = int(np.argmin(np.abs(np.log(g.theta_arcmin.values / 10.0))))
        near10[f"{i}{j}"] = {"theta": round(float(g.theta_arcmin.values[k]), 2), "xip": float(g.xip.values[k]),
                             "xip_err": float(g.xip_err.values[k])}
    return ok(f"Measured xi+/xi- for {len(pairs)} bin pair(s) x {nbins} theta bins ({computed_on}); "
              f"xi+ near 10 arcmin: " + ", ".join(f"{k}: {v['xip']:.2e}+/-{v['xip_err']:.1e}" for k, v in near10.items()),
              files, bin_pairs=pairs, xip_near_10arcmin=near10, per_bin=info["bins"], computed_on=computed_on,
              m_bias_applied=m if apply_m_bias else None, c_removed=remove_c,
              error_kind="jackknife" if npatch >= 2 else "shape noise only (no cosmic variance)")


measure_shear_correlation.weight = "heavy"


@tier(1)
@validate_call
def measure_nz_from_catalog(
    catalog_fits: Annotated[str, Field(description="Slim FITS from fetch_shear_catalog.")],
    output_dir: OutDir,
    z_column: Annotated[Literal["z_b"], Field(description="Photo-z column to histogram (BPZ z_B).")] = "z_b",
    dz: Annotated[float, Field(gt=0.005, le=0.2, description="Histogram bin width.")] = 0.05,
    compare_release: Annotated[bool, Field(description="Overlay the released KiDS-1000 SOM n(z) (tier 2).")] = True,
) -> ArtifactResult:
    """Weighted photo-z histogram n(z) per tomographic bin from a shear catalogue (tier 1 -> tier 2).

    Use for 'what redshifts does my subset probe'. This is the naive
    lensfit-weighted z_B histogram; by construction each bin is confined
    to its z_B edges and lacks the outlier tails that the released,
    SOM-calibrated n(z) (Hildebrandt+21, Wright+20) has - so its mean z
    sits below the calibrated one for the high bins. Treat the release
    n(z) as the analysis-grade product; this tool shows the selection.
    """
    import pandas as pd

    from ..datavectors import io as dvio

    c = _load(catalog_fits)
    z, w, tb = np.asarray(c[z_column], float), np.asarray(c["weight"], float), np.asarray(c["tomo_bin"]).astype(int)
    bins = sorted(int(b) for b in set(np.unique(tb).tolist()) - {0})
    edges = np.arange(0.0, 2.0 + dz / 2, dz)
    res = kernels.nz_histogram(z, w, tb, bins, edges)
    out = resolve_outdir(output_dir)
    cols = {"z": res["z"]} | {f"bin{b}": res["bins"][str(b)]["nz"] for b in bins}
    csv = out / f"{Path(catalog_fits).stem}_nz_{z_column}.csv"
    pd.DataFrame(cols).to_csv(csv, index=False)
    rel_mean = {}
    plt = mpl()
    fig, ax = plt.subplots(figsize=(6.4, 4))
    for k, b in enumerate(bins):
        ax.step(res["z"], res["bins"][str(b)]["nz"], where="mid", color=f"C{k}", lw=1.4, label=f"bin {b} (z_B hist)")
    if compare_release:
        try:
            nz = dvio.load_nz("kids_1000", "source")
            for k, b in enumerate(bins):
                n = np.asarray(nz["bins"][b - 1], float)
                zz = np.asarray(nz["z"], float)
                norm = np.trapezoid(n, zz)
                ax.plot(zz, n / norm, color=f"C{k}", ls="--", lw=1)
                rel_mean[str(b)] = round(float(np.trapezoid(zz * n, zz) / norm), 3)
        except Exception as exc:  # noqa: BLE001 - release n(z) optional (needs the 17 MB tarball)
            rel_mean = {"error": str(exc)[:200]}
    ax.set_xlabel("z")
    ax.set_ylabel("n(z) (normalised)")
    ax.set_xlim(0, 2)
    ax.set_title("solid: catalogue z_B histogram; dashed: released SOM n(z)" if compare_release else "")
    ax.legend(fontsize=7)
    png = save_fig(fig, out, f"{Path(catalog_fits).stem}_nz")
    means = {str(b): round(res["bins"][str(b)]["mean_z"], 3) for b in bins}
    return ok("Mean z_B per bin: " + ", ".join(f"{b}: {v}" for b, v in means.items())
              + (f"; released SOM n(z) means: {rel_mean}" if rel_mean else ""), [csv, png],
              mean_z=means, release_mean_z=rel_mean or None, counts={str(b): res["bins"][str(b)]["n"] for b in bins})


@tier(2)
@validate_call
def compare_measured_with_published_xi(
    measured_xip_csv: Annotated[str, Field(description="*_xip_measured.csv (or xim) from measure_shear_correlation.")],
    output_dir: OutDir,
    dataset: Annotated[Literal["kids_1000"], Field(description="Release whose xi to compare with.")] = "kids_1000",
    statistic: Annotated[Literal["xip", "xim"], Field(description="Which xi.")] = "xip",
    subset_area_deg2: Annotated[float | None, Field(gt=0, description="Effective (unmasked) area of your subset; "
                                                    "if given, also fit the amplitude with the RELEASE covariance "
                                                    "(noise + cosmic variance) scaled by 777.4 deg^2 / area.")] = None,
) -> ArtifactResult:
    """Overlay your own catalogue-measured xi+/- on the released KiDS-1000 xi+/- (both tier 2) per bin pair.

    Use after measure_shear_correlation to validate the measurement. The
    release covers ~1006 deg^2 (777 deg^2 effective), a subset only a few,
    so expect agreement within the subset's (much larger) errors plus
    cosmic variance, not point-by-point. Reports chi2 of (measured -
    released) using the measured errors only, per pair and in total, the
    ratio of weighted means over 5-60 arcmin, and one amplitude
    A = <measured / released shape> over all pairs at theta > 2' (A = 1
    means the subset reproduces the release; A = 0 means no signal).
    """
    import pandas as pd

    from ..datavectors import io as dvio

    me = pd.read_csv(measured_xip_csv)
    blk, cov, _ = dvio.load_statistic(dataset, statistic)
    rel = pd.DataFrame({"bin1": blk["bin1"], "bin2": blk["bin2"], "x": blk["x"], "value": blk["value"],
                        "error": np.sqrt(np.clip(np.diag(cov), 0, None))})
    pairs = sorted({(int(a), int(b)) for a, b in zip(me.bin1, me.bin2, strict=True)})
    plt = mpl()
    ncol = min(len(pairs), 5)
    nrow = -(-len(pairs) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.6 * nrow), squeeze=False, sharex=True)
    for ax in axes.flat[len(pairs):]:
        ax.set_visible(False)
    rows = []
    num = den = 0.0
    for ax, (i, j) in zip(axes.flat, pairs, strict=False):
        m = me[(me.bin1 == i) & (me.bin2 == j)].sort_values("x")
        r = rel[(rel.bin1 == min(i, j)) & (rel.bin2 == max(i, j))].sort_values("x")
        rv = np.interp(np.log(m.x), np.log(r.x), r.value)
        ok_ = (m.error > 0) & np.isfinite(m.value)
        chi2 = float(np.sum(((m.value - rv)[ok_] / m.error[ok_]) ** 2))
        fit = ok_ & (m.x > 2)  # amplitude of the released shape, small scales (theta < 2') excluded
        num += float(np.sum((m.value * rv / m.error ** 2)[fit]))
        den += float(np.sum((rv ** 2 / m.error ** 2)[fit]))
        sel = (m.x > 5) & (m.x < 60) & ok_
        wts = 1 / m.error[sel] ** 2
        ratio = float(np.sum(wts * m.value[sel]) / np.sum(wts * rv[sel])) if sel.any() else float("nan")
        rows.append({"bin1": i, "bin2": j, "chi2": round(chi2, 2), "ndof": int(ok_.sum()),
                     "ratio_5_60arcmin": round(ratio, 3)})
        ax.errorbar(r.x, 1e4 * r.x * r.value, 1e4 * r.x * r.error, fmt="s", ms=3, color="0.4",
                    label="released (1006 deg$^2$)")
        ax.errorbar(m.x * 1.03, 1e4 * m.x * m.value, 1e4 * m.x * m.error, fmt="o", ms=3, color="C3",
                    label="measured (subset)")
        ax.axhline(0, color="0.7", lw=0.5)
        ax.set_xscale("log")
        ax.set_title(f"bins {i},{j}: chi2={chi2:.1f}/{int(ok_.sum())}", fontsize=8)
        ax.set_xlabel(r"$\theta$ [arcmin]")
    for r_ in range(nrow):
        axes[r_][0].set_ylabel(rf"$\theta\,\xi_{{{'+' if statistic == 'xip' else '-'}}}\times10^4$ [arcmin]")
    axes[0][0].legend(fontsize=6)
    fig.tight_layout()
    out = resolve_outdir(output_dir)
    png = save_fig(fig, out, f"{Path(measured_xip_csv).stem}_vs_{dataset}")
    csv = out / f"{Path(measured_xip_csv).stem}_vs_{dataset}.csv"
    pd.DataFrame(rows).to_csv(csv, index=False)
    tot = sum(r["chi2"] for r in rows)
    nd = sum(r["ndof"] for r in rows)
    amp, amp_err = (num / den, den ** -0.5) if den > 0 else (float("nan"), float("nan"))
    cv = {}
    if subset_area_deg2:
        idx, mv = [], []
        for k in range(len(blk["value"])):
            i, j, x = int(blk["bin1"][k]), int(blk["bin2"][k]), float(blk["x"][k])
            g = me[(me.bin1 == i) & (me.bin2 == j)].sort_values("x")
            if x > 2 and len(g):
                idx.append(k)
                mv.append(float(np.interp(np.log(x), np.log(g.x), g.value)))
        rv_all, mv = np.asarray(blk["value"], float)[idx], np.asarray(mv)
        ci = np.linalg.inv(cov[np.ix_(idx, idx)] * (KIDS1000_AREA_EFF / subset_area_deg2))
        f = float(rv_all @ ci @ rv_all)
        cv = {"amplitude_cv": round(float(rv_all @ ci @ mv) / f, 3), "amplitude_cv_err": round(f ** -0.5, 3),
              "chi2_cv": round(float((mv - rv_all) @ ci @ (mv - rv_all)), 1), "ndof_cv": len(idx)}
    return ok(f"measured vs released {dataset} {statistic}: chi2 = {tot:.1f} / {nd} (measured errors only); "
              f"amplitude A = measured/released = {amp:.2f} +/- {amp_err:.2f} (theta > 2', diagonal errors)"
              + (f"; with area-scaled release covariance A = {cv['amplitude_cv']:.2f} +/- {cv['amplitude_cv_err']:.2f}, "
                 f"chi2 = {cv['chi2_cv']:.1f}/{cv['ndof_cv']}" if cv else "") + ".",
              [png, csv], per_pair=rows, chi2_total=round(tot, 2), ndof=nd, amplitude=round(amp, 3), **cv,
              amplitude_err=round(amp_err, 3),
              note="Released xi is m-corrected over ~1006 deg^2; subset errors exclude cosmic variance.")
