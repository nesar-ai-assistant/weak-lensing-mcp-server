"""MCMC chains: fetch, marginal summaries (heavy + dispatchable), S8-Om contours, S8 tension."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

import numpy as np
from pydantic import Field, validate_call

from ..common import ArtifactResult, data_cache_dir, mpl, ok, resolve_outdir, save_fig
from ..datavectors import io as dvio
from ..registry import PLANCK18_S8, PUBLISHED_S8
from . import kernels

__all__ = ["fetch_wl_chain", "analyze_wl_chain", "plot_s8_om_contours", "compare_s8_constraints"]

DISPATCH_KERNELS = {
    "chain_marginals": {
        "function": "chains.kernels.chain_marginals",
        "args_shape": "{'url': <chain .txt or release tarball URL>, 'member': <path suffix in tarball or None>, "
                      "'params': ['S_8','omega_m','sigma_8'] or None}",
        "pip_deps": ["numpy", "certifi"],
        "duration_hint_s": 900,
        "returns": "n_samples, n_eff, per-parameter mean/std/median/68%/95%/mode/HPD68, corr_S8_Om, "
                   "max_post_sample; the downloaded chain stays in the job cwd",
    },
}

Dataset = Literal["des_y3", "des_y6", "kids_1000", "kids_legacy", "hsc_y3"]
OutDir = Annotated[str, Field(min_length=1, description="Directory where output files are written.")]
LABELS = {"des_y3": "DES Y3", "des_y6": "DES Y6", "kids_1000": "KiDS-1000",
          "kids_legacy": "KiDS-Legacy", "hsc_y3": "HSC Y3"}


def chain_source(dataset: str, chain: str) -> tuple[str, str | None]:
    """(url, tarball member or None) for a chain key."""
    meta = dvio.check_dataset(dataset)
    chains = meta.get("chains", {})
    if chain not in chains:
        raise ValueError(f"{dataset} chains: {list(chains)}; got {chain!r}")
    rel = chains[chain]
    if rel == "credentialed":
        dest = data_cache_dir() / dataset
        raise RuntimeError(
            f"{dataset} chain {chain!r} is behind HSC-SSP credentials (hsc-release.mtk.nao.ac.jp returns 401). "
            f"Log in at https://hsc-release.mtk.nao.ac.jp/doc/index.php/wly3/, download the chain "
            f"(dalal23/hsc_y3_fourier_cosmic_shear.txt or li23/outputs/...), save it under {dest}/ as "
            f"{chain}.txt, then call fetch_wl_chain again. compare_s8_constraints uses the published values "
            f"meanwhile.")
    if "archive" in meta:
        return meta["archive"], rel
    from ..registry import DES_Y3_CHAINS
    return DES_Y3_CHAINS + rel, None


def chain_local_path(dataset: str, chain: str) -> Path:
    manual = data_cache_dir() / dataset / f"{chain}.txt"
    if manual.is_file():
        return manual
    url, member = chain_source(dataset, chain)
    if member:
        return dvio.archive_member(dataset, member)
    return dvio.download(url, subdir=dataset, timeout=900)


@validate_call
def fetch_wl_chain(
    dataset: Annotated[Dataset, Field(description="Survey release key.")],
    chain: Annotated[str, Field(description="Chain key, see describe_wl_dataset(...)['chains']: des_y3 "
                                            "3x2pt_lcdm, 3x2pt_wcdm, 1x2pt_lcdm, 2x2pt_lcdm, planck_lcdm, ...; "
                                            "kids_* cosebis, bandpowers, xipm.")],
    output_dir: OutDir,
    max_rows: Annotated[int, Field(ge=1000, le=2_000_000, description="Cap on rows written to the CSV "
                                                                      "(highest-weight rows kept).")] = 200_000,
) -> ArtifactResult:
    """Download a published MCMC/nested-sampling chain and write a slim CSV.

    Use before analyze_wl_chain / plot_s8_om_contours. The CSV holds
    weight, omega_m, sigma_8, S_8 (derived as sigma_8 sqrt(Om/0.3) where the
    chain lacks it or has NaN) and h0 / w when present. DES chains are
    5-50 MB each; KiDS chains come from the release tarball. HSC Y3 chains
    need HSC-SSP credentials (clear error + manual-download path).
    """
    import pandas as pd

    path = chain_local_path(dataset, chain)
    cols, data = kernels.parse_cosmosis_text(path)
    sc = kernels.standard_columns(cols, data)
    df = pd.DataFrame(sc)
    if len(df) > max_rows:
        df = df.nlargest(max_rows, "weight")
    out = resolve_outdir(output_dir)
    csv = out / f"{dataset}_{chain}_chain.csv"
    df.to_csv(csv, index=False)
    return ok(f"{LABELS[dataset]} {chain}: {len(df)} samples.", [csv], dataset=dataset, chain=chain,
              source_file=str(path), columns=list(df.columns), n_samples=len(df),
              n_raw_columns=len(cols))


def _fmt(p: dict) -> str:
    return f"{p['mean']:.3f} +/- {p['std']:.3f}"


@validate_call
def analyze_wl_chain(
    dataset: Annotated[Dataset, Field(description="Survey release key.")],
    chain: Annotated[str, Field(description="Chain key (as fetch_wl_chain).")],
    output_dir: OutDir,
) -> ArtifactResult:
    """Weighted marginal constraints on S8, Omega_m, sigma8 (and h, w if sampled) from a chain.

    Use to quote a survey's parameter constraints. Reports mean +/- std,
    median with 16/84% quantiles, 95% intervals, and the marginal mode with
    68% HPD (KiDS papers quote MAP/mode + PJ-HPD; DES quotes mean +/- std).
    Heavy: downloads tens of MB; with set_dispatch it runs the
    chains.kernels.chain_marginals kernel on the facility instead.
    Sanity: des_y3 3x2pt_lcdm S8 = 0.776 +/- 0.017.
    """
    import json

    from mcp_server.dispatch import remote_site, run_kernel  # lazy: server-only

    site = remote_site()
    if site:
        url, member = chain_source(dataset, chain)
        res = run_kernel("chains.kernels.chain_marginals", {"url": url, "member": member, "params": None},
                         pip_deps=DISPATCH_KERNELS["chain_marginals"]["pip_deps"], duration=900)
        payload, computed_on = res["result"], res.get("host", site)
    else:
        path = chain_local_path(dataset, chain)
        cols, data = kernels.parse_cosmosis_text(path)
        payload, computed_on = kernels.summarize_columns(kernels.standard_columns(cols, data)), "local"
    out = resolve_outdir(output_dir)
    js = out / f"{dataset}_{chain}_marginals.json"
    js.write_text(json.dumps(payload, indent=2))
    P = payload["params"]
    msg = f"{LABELS[dataset]} {chain}: S8 = {_fmt(P['S_8'])}, Om = {_fmt(P['omega_m'])}, sigma8 = {_fmt(P['sigma_8'])}"
    return ok(msg, [js], dataset=dataset, chain=chain, computed_on=computed_on, **payload)


analyze_wl_chain.weight = "heavy"


def _kde2(x, y, w, n=120):
    """Weighted 2D histogram smoothed with a Gaussian kernel; returns grid and density."""
    xr = np.quantile(x, [0.001, 0.999])
    yr = np.quantile(y, [0.001, 0.999])
    padx, pady = 0.15 * np.ptp(xr), 0.15 * np.ptp(yr)
    H, xe, ye = np.histogram2d(x, y, bins=n, range=[[xr[0] - padx, xr[1] + padx], [yr[0] - pady, yr[1] + pady]],
                               weights=w)
    from scipy.ndimage import gaussian_filter
    H = gaussian_filter(H, 2.5)
    return 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H.T


def _levels(H, masses=(0.95, 0.68)):
    f = np.sort(H.ravel())[::-1]
    c = np.cumsum(f) / f.sum()
    return sorted(float(f[np.searchsorted(c, m)]) for m in masses)


@validate_call
def plot_s8_om_contours(
    chain_csvs: Annotated[list[str], Field(min_length=1, description="Chain CSVs from fetch_wl_chain.")],
    output_dir: OutDir,
    labels: Annotated[list[str] | None, Field(description="Legend labels (default: file stems).")] = None,
    y_param: Annotated[Literal["S_8", "sigma_8"], Field(description="y-axis parameter.")] = "S_8",
    show_planck: Annotated[bool, Field(description="Draw the Planck 2018 S8 band.")] = True,
) -> ArtifactResult:
    """Overlay 68/95% contours in the Omega_m - S8 (or sigma8) plane for several surveys.

    Use to reproduce the classic multi-survey S8-Om figure (DES Y3 Fig. 7 /
    Asgari+21 Fig. 5 style). Contours come from a weighted, smoothed 2D
    histogram (no getdist needed). The Planck S8 band is 0.832 +/- 0.013.
    """
    import pandas as pd

    plt = mpl()
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    handles = []
    for k, path in enumerate(chain_csvs):
        df = pd.read_csv(path)
        lab = labels[k] if labels and k < len(labels) else Path(path).stem.replace("_chain", "")
        x, y, H = _kde2(df.omega_m.values, df[y_param].values, df.weight.values)
        col = f"C{k}"
        lv = _levels(H)
        ax.contourf(x, y, H, levels=[lv[0], lv[1], H.max() * 1.01], colors=[col, col], alpha=0.25)
        ax.contour(x, y, H, levels=lv, colors=col, linewidths=1.2)
        handles.append(plt.Line2D([], [], color=col, lw=2, label=lab))
    if show_planck and y_param == "S_8":
        m, s = PLANCK18_S8["mean"], PLANCK18_S8["sigma"]
        ax.axhspan(m - s, m + s, color="0.5", alpha=0.25)
        handles.append(plt.Line2D([], [], color="0.5", lw=6, alpha=0.4, label="Planck 2018 S8"))
    ax.set_xlabel(r"$\Omega_{\rm m}$")
    ax.set_ylabel(r"$S_8$" if y_param == "S_8" else r"$\sigma_8$")
    ax.legend(handles=handles, fontsize=8, loc="upper right")
    out = resolve_outdir(output_dir)
    png = save_fig(fig, out, f"{y_param}_om_contours")
    return ok(f"Contours for {len(chain_csvs)} chains.", [png], y_param=y_param)


@validate_call
def compare_s8_constraints(
    output_dir: OutDir,
    marginals_json: Annotated[list[str] | None, Field(description="Optional analyze_wl_chain JSONs to include "
                                                                  "(our own chain re-analysis).")] = None,
    include_published: Annotated[bool, Field(description="Also include the published headline S8 values "
                                                         "(needed for HSC Y3, whose chains are credentialed).")] = True,
    reference_mean: Annotated[float, Field(description="Reference S8 (default Planck 2018).")] = 0.832,
    reference_sigma: Annotated[float, Field(gt=0, description="Reference S8 1-sigma.")] = 0.013,
) -> ArtifactResult:
    """Whisker plot of S8 across surveys vs Planck (0.832 +/- 0.013) and the tension in sigma.

    Use for 'is there an S8 tension' questions. Tension per entry is the
    Gaussian difference-in-means (S8_ref - S8) / sqrt(sigma^2 + sigma_ref^2),
    taking the survey error on the side facing the reference for
    asymmetric errors. This is the simple 1D estimate the papers quote as
    a first number (DES Y3 3x2pt ~2.3 sigma, KiDS-1000 ~3 sigma, KiDS-Legacy
    < 1 sigma); full-posterior metrics differ at the few-0.1 sigma level.
    """
    import json

    import pandas as pd

    rows = []
    if include_published:
        for name, v in PUBLISHED_S8.items():
            rows.append({"label": name, "S8": v["S8"], "plus": v["plus"], "minus": v["minus"],
                         "source": "published: " + v["ref"]})
    for path in marginals_json or []:
        d = json.loads(Path(path).read_text())
        p = d["params"]["S_8"]
        stem = Path(path).stem.replace("_marginals", "")
        rows.append({"label": stem + " (this server)", "S8": p["mean"], "plus": p["std"], "minus": p["std"],
                     "source": "chain re-analysis (mean +/- std)"})
    if not rows:
        raise ValueError("Nothing to compare: pass marginals_json and/or include_published=True")
    for r in rows:
        side = r["plus"] if reference_mean > r["S8"] else r["minus"]
        r["tension_sigma"] = round(float(reference_mean - r["S8"]) / float(np.hypot(side, reference_sigma)), 2)
    df = pd.DataFrame(rows)
    out = resolve_outdir(output_dir)
    csv = out / "s8_comparison.csv"
    df.to_csv(csv, index=False)

    plt = mpl()
    fig, ax = plt.subplots(figsize=(6.2, 0.45 * len(df) + 1.4))
    yy = np.arange(len(df))[::-1]
    ax.errorbar(df.S8, yy, xerr=[df.minus, df.plus], fmt="o", color="k", capsize=3)
    ax.axvspan(reference_mean - reference_sigma, reference_mean + reference_sigma, color="C1", alpha=0.3,
               label=f"Planck 2018 S8 = {reference_mean} +/- {reference_sigma}")
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{lab}  ({t:+.1f}$\\sigma$)" for lab, t in zip(df.label, df.tension_sigma, strict=True)],
                       fontsize=8)
    ax.set_xlabel(r"$S_8 = \sigma_8\sqrt{\Omega_{\rm m}/0.3}$")
    ax.legend(fontsize=7, loc="lower left")
    png = save_fig(fig, out, "s8_whisker")
    summary = {r["label"]: r["tension_sigma"] for r in rows}
    return ok("S8 tension vs reference (sigma): " + ", ".join(f"{k} {v}" for k, v in summary.items()),
              [csv, png], tension_sigma=summary, reference={"mean": reference_mean, "sigma": reference_sigma,
                                                             "citation": PLANCK18_S8["citation"]})

