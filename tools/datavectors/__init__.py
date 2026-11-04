"""Weak-lensing 2pt data vectors: fetch, n(z), covariance, S/N, xi+- triangle plots."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

import numpy as np
from pydantic import Field, validate_call

from ..common import ArtifactResult, mpl, ok, resolve_outdir, save_fig
from ..registry import DATASETS
from . import io

__all__ = ["fetch_wl_datavector", "plot_shear_correlation", "plot_redshift_distributions",
           "plot_covariance", "compute_signal_to_noise"]

DISPATCH_KERNELS: dict = {}

Dataset = Literal["des_y3", "des_y6", "kids_1000", "kids_legacy", "hsc_y3"]
OutDir = Annotated[str, Field(min_length=1, description="Directory where output files are written.")]

LABELS = {"des_y3": "DES Y3", "des_y6": "DES Y6", "kids_1000": "KiDS-1000",
          "kids_legacy": "KiDS-Legacy", "hsc_y3": "HSC Y3"}


def _read_csv(path: str) -> dict:
    import pandas as pd
    return pd.read_csv(path)


@validate_call
def fetch_wl_datavector(
    dataset: Annotated[Dataset, Field(description="Survey release key (see list_wl_datasets).")],
    statistic: Annotated[str, Field(description="Statistic: xip, xim (all); gammat, wtheta (des_y3); cosebis, "
                                                "bandpowers (kids_*; plus cosebis_b, bandpowers_b for "
                                                "kids_legacy); cl_ee (hsc_y3).")],
    output_dir: OutDir,
) -> ArtifactResult:
    """Download one published 2pt statistic and write a tidy CSV.

    Use this first for any data-vector question. Columns: bin1, bin2
    (tomographic bin pair, 1-based; for gammat bin1=lens, bin2=source),
    angbin, x (theta in arcmin, COSEBI mode n, or ell - see x_name),
    value, error (sqrt of the covariance diagonal). Also writes the
    statistic's covariance sub-matrix as .npy for plot_covariance /
    compute_signal_to_noise. Data vectors are the published ones (scale
    cuts NOT applied: DES/KiDS files are the full measured vectors; the HSC real-space
    file is the CosmoSIS-library one, wider than the Li+23 analysis cuts).
    KiDS downloads a 50-130 MB release tarball on first use.
    """
    import pandas as pd

    blk, cov, ext = io.load_statistic(dataset, statistic)
    out = resolve_outdir(output_dir)
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    df = pd.DataFrame({"bin1": blk["bin1"], "bin2": blk["bin2"], "angbin": blk["angbin"],
                       "x": blk["x"], "value": blk["value"], "error": err})
    stem = f"{dataset}_{statistic}"
    csv = out / f"{stem}.csv"
    df.to_csv(csv, index=False)
    covp = out / f"{stem}_cov.npy"
    np.save(covp, cov)
    pairs = sorted(set(zip(df.bin1.tolist(), df.bin2.tolist(), strict=True)))
    return ok(f"{LABELS[dataset]} {statistic}: {len(df)} points, {len(pairs)} bin pairs.", [csv, covp],
              dataset=dataset, statistic=statistic, extension=ext, x_name=io.X_LABEL.get(ext, "x"),
              n_points=len(df), n_bin_pairs=len(pairs), columns=list(df.columns),
              citation=DATASETS[dataset]["citation"])


@validate_call
def plot_shear_correlation(
    xip_csv: Annotated[str, Field(description="CSV from fetch_wl_datavector(statistic='xip').")],
    output_dir: OutDir,
    xim_csv: Annotated[str | None, Field(description="Optional xim CSV; drawn in the lower triangle like the "
                                                     "DES/KiDS papers.")] = None,
    scale_theta: Annotated[bool, Field(description="Plot theta * xi * 1e4 (paper convention) instead of xi.")] = True,
    label: Annotated[str, Field(description="Title / legend label.")] = "",
) -> ArtifactResult:
    """Triangle-of-bin-pairs plot of xi_+ (upper) and xi_- (lower), as in Amon+22 Fig. 1 / Asgari+21.

    Use after fetch_wl_datavector. Each panel is one tomographic bin pair
    (i,j); points carry the diagonal-covariance error bars. Works for any
    statistic CSV (COSEBIs, band powers, C_ell) - then only the upper
    triangle is drawn and x is log-scaled when it spans > 1 decade.
    """
    plt = mpl()
    p = _read_csv(xip_csv)
    m = _read_csv(xim_csv) if xim_csv else None
    nb = int(max(p.bin1.max(), p.bin2.max()))
    fig, axes = plt.subplots(nb, nb, figsize=(2.2 * nb, 2.0 * nb), sharex=True)
    axes = np.atleast_2d(axes)
    for ax in axes.ravel():
        ax.set_visible(False)

    def draw(ax, sub, color, name):
        x, y, e = sub.x.values, sub.value.values, sub.error.values
        f = x * 1e4 if scale_theta else np.ones_like(x)
        ax.set_visible(True)
        ax.errorbar(x, f * y, f * e, fmt="o", ms=2.5, color=color, lw=0.8, label=name)
        ax.axhline(0, color="0.6", lw=0.5)
        if np.nanmax(x) / max(np.nanmin(x), 1e-9) > 10:
            ax.set_xscale("log")
        ax.tick_params(labelsize=6)

    for (i, j), sub in p.groupby(["bin1", "bin2"]):
        a, b = sorted((int(i), int(j)))
        ax = axes[a - 1, b - 1]
        draw(ax, sub, "C0", "xi+")
        ax.text(0.05, 0.85, f"{a},{b}", transform=ax.transAxes, fontsize=7)
    if m is not None:
        for (i, j), sub in m.groupby(["bin1", "bin2"]):
            a, b = sorted((int(i), int(j)))
            ax = axes[b - 1, a - 1]
            draw(ax, sub, "C3", "xi-")
            ax.text(0.05, 0.85, f"{a},{b}", transform=ax.transAxes, fontsize=7, color="C3")
    ylab = r"$\theta\,\xi\times10^4$" if scale_theta else r"value"
    fig.supxlabel(r"$\theta$ [arcmin]" if "theta" in xip_csv or "xi" in xip_csv else "x")
    fig.supylabel(ylab + (r"  ($\xi_+$ upper, $\xi_-$ lower)" if m is not None else ""))
    fig.suptitle(label or Path(xip_csv).stem)
    fig.tight_layout()
    out = resolve_outdir(output_dir)
    png = save_fig(fig, out, f"{Path(xip_csv).stem}_triangle")
    return ok(f"Triangle plot of {nb} tomographic bins.", [png], n_bins=nb, lower_triangle=m is not None)


@validate_call
def plot_redshift_distributions(
    datasets: Annotated[list[Dataset], Field(min_length=1, description="One or more dataset keys.")],
    output_dir: OutDir,
    sample: Annotated[Literal["source", "lens"], Field(description="'lens' exists only for des_y3.")] = "source",
) -> ArtifactResult:
    """Plot n(z) per tomographic bin (one panel per dataset) and write them to CSV.

    Use to inspect or compare redshift distributions across surveys.
    Curves are normalised to unit integral; the metadata reports each
    bin's mean redshift (sanity: DES Y3 sources ~0.34, 0.52, 0.74, 0.96;
    KiDS-1000 ~0.26..0.93; HSC Y3 ~0.44..1.3).
    """
    import pandas as pd

    plt = mpl()
    out = resolve_outdir(output_dir)
    fig, axes = plt.subplots(1, len(datasets), figsize=(4.2 * len(datasets), 3.2), squeeze=False)
    files, means = [], {}
    for ax, ds in zip(axes[0], datasets, strict=True):
        nz = io.load_nz(ds, sample)
        z = nz["z"]
        cols, zbar = {"z": z}, []
        for k, n in enumerate(nz["bins"], 1):
            norm = np.trapezoid(n, z)
            nn = n / norm if norm > 0 else n
            cols[f"bin{k}"] = nn
            zbar.append(round(float(np.trapezoid(z * nn, z)), 3))
            ax.plot(z, nn, label=f"bin {k}", lw=1.3)
        means[ds] = zbar
        csv = out / f"{ds}_nz_{sample}.csv"
        pd.DataFrame(cols).to_csv(csv, index=False)
        files.append(csv)
        ax.set_title(f"{LABELS[ds]} {sample} n(z)")
        ax.set_xlabel("z")
        ax.set_xlim(0, min(float(z.max()), 3.0))
        ax.legend(fontsize=7)
    axes[0][0].set_ylabel("n(z) (normalised)")
    fig.tight_layout()
    files.append(save_fig(fig, out, f"nz_{sample}_{'_'.join(datasets)}"))
    return ok(f"n(z) for {len(datasets)} dataset(s).", files, mean_z=means)


@validate_call
def plot_covariance(
    cov_npy: Annotated[str, Field(description="Covariance .npy from fetch_wl_datavector.")],
    output_dir: OutDir,
    title: Annotated[str, Field(description="Plot title.")] = "",
) -> ArtifactResult:
    """Plot the correlation matrix r_ij = C_ij / sqrt(C_ii C_jj) of a data vector.

    Use to inspect off-diagonal structure (strong positive correlation
    between neighbouring angular bins and across bin pairs at large theta
    is expected for cosmic shear).
    """
    plt = mpl()
    cov = np.load(cov_npy)
    d = np.sqrt(np.clip(np.diag(cov), 1e-300, None))
    corr = cov / np.outer(d, d)
    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, origin="lower")
    fig.colorbar(im, ax=ax, label="correlation")
    ax.set_title(title or Path(cov_npy).stem)
    ax.set_xlabel("data-vector index")
    ax.set_ylabel("data-vector index")
    out = resolve_outdir(output_dir)
    png = save_fig(fig, out, f"{Path(cov_npy).stem}_corr")
    off = corr[~np.eye(len(corr), dtype=bool)]
    eig = np.linalg.eigvalsh(cov)
    return ok(f"Correlation matrix {corr.shape[0]}x{corr.shape[0]}.", [png], shape=list(corr.shape),
              max_offdiag=float(off.max()) if off.size else 0.0, mean_abs_offdiag=float(np.abs(off).mean()) if off.size else 0.0,
              positive_definite=bool(eig.min() > 0), condition_number=float(eig.max() / eig.min()) if eig.min() > 0 else None)


@validate_call
def compute_signal_to_noise(
    datasets: Annotated[list[Dataset], Field(min_length=1, description="Dataset keys.")],
    statistics: Annotated[list[str], Field(min_length=1, description="Statistics to combine per dataset, e.g. "
                                                                    "['xip','xim'] (joint, with cross-covariance).")],
    output_dir: OutDir,
) -> ArtifactResult:
    """S/N = sqrt(d^T C^-1 d) of the published data vector, using its full covariance.

    Use for 'how significant is the detection' questions. Statistics
    listed together are combined with their cross-covariance (they must
    live in the same file, e.g. xip+xim, or all four DES 3x2pt blocks).
    No scale cuts are applied (DES files are uncut), so this is the raw
    measurement S/N. Writes a CSV summary.
    """
    import pandas as pd

    rows = []
    for ds in datasets:
        io.check_dataset(ds)
        stats = [s for s in statistics if s in io.STATISTICS[ds]]
        if not stats:
            raise ValueError(f"{ds}: none of {statistics} available (has {list(io.STATISTICS[ds])})")
        groups: dict[str, list[str]] = {}
        for s in stats:
            groups.setdefault(io.STATISTICS[ds][s][0], []).append(s)
        for key, ss in groups.items():
            dv = io.load_product(ds, key)
            idx = np.concatenate([np.arange(dv["cov"].shape[0])[io.block_slice(dv, io.STATISTICS[ds][s][1])]
                                  for s in ss])
            d = np.concatenate([dv["blocks"][io.STATISTICS[ds][s][1]]["value"] for s in ss])
            c = dv["cov"][np.ix_(idx, idx)]
            chi2 = float(d @ np.linalg.solve(c, d))
            rows.append({"dataset": ds, "statistics": "+".join(ss), "n_points": len(d),
                         "snr": round(float(np.sqrt(max(chi2, 0))), 2),
                         "snr_diag_only": round(float(np.sqrt(np.sum(d ** 2 / np.diag(c)))), 2)})
    out = resolve_outdir(output_dir)
    csv = out / "signal_to_noise.csv"
    pd.DataFrame(rows).to_csv(csv, index=False)
    msg = "; ".join(f"{r['dataset']} {r['statistics']}: S/N={r['snr']}" for r in rows)
    return ok(msg, [csv], results=rows)
