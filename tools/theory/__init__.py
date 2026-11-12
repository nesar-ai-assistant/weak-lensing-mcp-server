"""Cosmic-shear theory predictions: xi+/xi- from n(z) + cosmology, optional baryonic suppression."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import numpy as np
from pydantic import Field, validate_call

from ..common import ArtifactResult, data_kind, mpl, ok, resolve_outdir, save_fig
from . import kernels

__all__ = ["predict_shear_correlation"]

DISPATCH_KERNELS: dict = {}

OutDir = Annotated[str, Field(min_length=1, description="Directory where output files are written.")]
SERVED_NZ = ("kids_1000", "kids_legacy", "des_y3", "hsc_y3")


def load_user_nz(path: str) -> dict:
    """Text/CSV n(z): first column z (bin centres), one column per tomographic bin ('#' comments allowed)."""
    import pandas as pd

    p = Path(path)
    if p.suffix.lower() == ".csv":
        df = pd.read_csv(p, comment="#")
        arr = df.to_numpy(float)
    else:
        arr = np.loadtxt(p, comments="#", ndmin=2)
    if arr.shape[1] < 2:
        raise ValueError(f"{path}: need columns z, n_1(z), n_2(z), ...")
    return {"z": arr[:, 0], "bins": [arr[:, i] for i in range(1, arr.shape[1])]}


def _nz(nz_source: str) -> tuple[dict, str]:
    if nz_source in SERVED_NZ:
        from ..datavectors import io
        return io.load_nz(nz_source), nz_source
    if Path(nz_source).expanduser().is_file():
        return load_user_nz(str(Path(nz_source).expanduser())), "user"
    raise ValueError(f"nz_source must be one of {SERVED_NZ} or an existing file (z, n_1, n_2, ...); got {nz_source!r}")


@data_kind("summary_statistic")
@validate_call
def predict_shear_correlation(
    output_dir: OutDir,
    nz_source: Annotated[str, Field(description="n(z): 'kids_1000', 'kids_legacy', 'des_y3', 'hsc_y3' (the "
                                                "release n(z) this server serves) or a path to a text/CSV file "
                                                "with columns z, n_1(z), n_2(z), ...")] = "kids_1000",
    theta_arcmin: Annotated[list[float] | None, Field(description="Angular separations [arcmin]; default 20 "
                                                                  "log-spaced values 0.5-300.")] = None,
    Om: Annotated[float, Field(gt=0.05, lt=0.9, description="Omega_m (total matter, incl. neutrinos).")] = 0.3,
    Ob: Annotated[float, Field(gt=0.01, lt=0.2, description="Omega_b.")] = 0.049,
    h: Annotated[float, Field(gt=0.4, lt=1.0, description="H0 / 100.")] = 0.674,
    ns: Annotated[float, Field(gt=0.8, lt=1.2, description="Scalar spectral index.")] = 0.965,
    sigma8: Annotated[float | None, Field(gt=0.3, lt=1.5, description="sigma8 (give this or S8).")] = None,
    S8: Annotated[float | None, Field(gt=0.3, lt=1.5, description="S8 = sigma8 sqrt(Om/0.3); used if sigma8 "
                                                                  "is not given.")] = None,
    mnu: Annotated[float, Field(ge=0, le=1, description="Sum of neutrino masses [eV].")] = 0.06,
    bin_pairs: Annotated[list[list[int]] | None, Field(description="1-based tomographic pairs, e.g. [[1,1],[1,2]]; "
                                                                   "default all i<=j.")] = None,
    baryon_file: Annotated[str | None, Field(description="Optional P_hydro/P_dmo ratio CSV multiplied into "
                                                         "P(k): hydroemu-mcp-server Pk-ratio output ('# ' header "
                                                         "lines, columns x = k [h/Mpc], y_mean = ratio) or plain "
                                                         "k[h/Mpc], z, ratio (z-dependent) / k, ratio. Without z "
                                                         "the ratio is applied at all redshifts.")] = None,
    halofit: Annotated[str, Field(description="CAMB halofit_version: takahashi, mead2020, mead (HMcode2016), "
                                              "original, bird.")] = "takahashi",
    compare_dataset: Annotated[str | None, Field(description="Overlay the served xi+ data vector of this dataset "
                                                             "(same n(z) recommended) and fit a single amplitude "
                                                             "A = data / theory.")] = None,
    plot: bool = True,
) -> ArtifactResult:
    """Predict cosmic-shear xi+(theta), xi-(theta) per tomographic bin pair (CAMB halofit + Limber + Hankel), optional baryons.

    Pure GG lensing signal: no intrinsic alignments, shear bias, or photo-z
    shifts. P_NL(k,z) from CAMB (halofit, normalised to sigma8);
    C_ell by the Limber approximation (k = (ell+1/2)/chi) in flat
    cosmology; xi+- by dense Bessel J0/J4 quadrature to ell = 1e5. With
    baryon_file the ratio S(k[,z]) multiplies P(k) (S -> 1 below the
    tabulated k, held at its last value above), and the output also gives
    the DMO prediction and the ratio xi_baryon / xi_dmo. Writes
    theory_xipm.csv (bin1, bin2, theta_arcmin, xip, xim[, xip_dmo, xim_dmo,
    ratio_xip, ratio_xim]) and, with plot, a triangle figure. Needs camb
    (pip install camb). ~5-15 s. Validation: KiDS-1000 n(z), Om=0.27,
    S8=0.76 -> served KiDS-1000 xi+ / theory = 0.90 +/- 0.03 (0.5-75
    arcmin; the remaining deficit is mostly intrinsic alignments, absent
    here).
    """
    import pandas as pd

    if sigma8 is None and S8 is None:
        sigma8 = 0.811
    if sigma8 is None:
        sigma8 = S8 / np.sqrt(Om / 0.3)
    th = np.asarray(theta_arcmin if theta_arcmin else np.geomspace(0.5, 300.0, 20), float)
    if np.any(th <= 0):
        raise ValueError("theta_arcmin must be positive")
    nz, nz_label = _nz(nz_source)
    nb = len(nz["bins"])
    if bin_pairs:
        for p in bin_pairs:
            if len(p) != 2 or not (1 <= p[0] <= nb and 1 <= p[1] <= nb):
                raise ValueError(f"bin pair {p} out of range 1..{nb}")
        pairs = [(min(p) - 1, max(p) - 1) for p in bin_pairs]
    else:
        pairs = None
    bary = kernels.read_baryon_ratio(baryon_file) if baryon_file else None
    cosmo = {"Om": Om, "Ob": Ob, "h": h, "ns": ns, "sigma8": float(sigma8), "mnu": mnu}
    r = kernels.predict_xi(nz["z"], nz["bins"], th, cosmo, pairs, bary, halofit=halofit)
    rows = []
    for n, (i, j) in enumerate(r["pairs"]):
        for k, t in enumerate(th):
            row = {"bin1": i + 1, "bin2": j + 1, "theta_arcmin": t, "xip": r["xip"][n, k], "xim": r["xim"][n, k]}
            if bary is not None:
                row.update(xip_dmo=r["xip_dmo"][n, k], xim_dmo=r["xim_dmo"][n, k],
                           ratio_xip=r["xip"][n, k] / r["xip_dmo"][n, k],
                           ratio_xim=r["xim"][n, k] / r["xim_dmo"][n, k])
            rows.append(row)
    df = pd.DataFrame(rows)
    out = resolve_outdir(output_dir)
    tag = f"{nz_label}_Om{Om:.3f}_s8{sigma8:.3f}" + ("_baryons" if bary is not None else "")
    csv = out / f"theory_xipm_{tag}.csv"
    with open(csv, "w") as fh:
        fh.write(f"# predict_shear_correlation: nz={nz_source} Om={Om} Ob={Ob} h={h} ns={ns} sigma8={sigma8:.4f} "
                 f"S8={sigma8 * np.sqrt(Om / 0.3):.4f} mnu={mnu} halofit={halofit} baryons={baryon_file}\n")
        df.to_csv(fh, index=False)
    files = [csv]
    meta = {"cosmology": {**cosmo, "S8": round(float(sigma8 * np.sqrt(Om / 0.3)), 4)}, "nz_source": nz_source,
            "n_bins": nb, "bin_pairs": [[i + 1, j + 1] for i, j in r["pairs"]], "theta_arcmin": th.tolist(),
            "As": r["camb"]["As"], "halofit": halofit,
            "model": "GG only (no IA / shear bias / photo-z nuisance); CAMB halofit + Limber + J0/J4 quadrature"}
    msg_extra = ""
    if bary is not None:
        last = df[(df.bin1 == nb) & (df.bin2 == nb)] if bin_pairs is None else df[(df.bin1 == r["pairs"][-1][0] + 1)
                                                                                  & (df.bin2 == r["pairs"][-1][1] + 1)]
        k5 = int(np.argmin(np.abs(th - 5.0)))
        meta["baryons"] = {"file": baryon_file, "z_dependent": bary["z"] is not None,
                           "min_ratio_xip": round(float(df.ratio_xip.min()), 4),
                           "min_ratio_xim": round(float(df.ratio_xim.min()), 4),
                           "ratio_xip_near_5arcmin_last_pair": round(float(last.ratio_xip.iloc[k5]), 4)}
        msg_extra = (f" Baryons: xi+ suppressed to {meta['baryons']['min_ratio_xip']:.3f}x DMO at most, "
                     f"xi- to {meta['baryons']['min_ratio_xim']:.3f}x.")
    comp = None
    if compare_dataset:
        comp = _compare(r, th, compare_dataset)
        meta["comparison"] = comp
        msg_extra += (f" vs {compare_dataset} xi+: A = data/theory = {comp['A']:.3f} +/- {comp['A_err']:.3f} "
                      f"({comp['theta_range'][0]}-{comp['theta_range'][1]} arcmin, diag. errors).")
    if plot:
        files.append(_plot(r, th, nb, out, f"theory_xipm_{tag}", bary is not None, comp))
    return ok(f"Predicted xi+/xi- for {len(r['pairs'])} bin pair(s) x {len(th)} theta (S8 = "
              f"{meta['cosmology']['S8']}, n(z) {nz_source}).{msg_extra}", files, **meta)


def _compare(r: dict, th: np.ndarray, dataset: str, tmin: float = 0.5, tmax: float = 75.0) -> dict:
    from ..datavectors import io
    blk, cov, _ = io.load_statistic(dataset, "xip")
    err = np.sqrt(np.diag(cov))
    num = den = 0.0
    per = {}
    for n, (i, j) in enumerate(r["pairs"]):
        s = (blk["bin1"] == i + 1) & (blk["bin2"] == j + 1)
        if not s.any():
            continue
        x, d, e = blk["x"][s], blk["value"][s], err[s]
        sel = (x >= tmin) & (x <= tmax) & (x >= th.min()) & (x <= th.max())
        if not sel.any():
            continue
        m = np.interp(np.log(x[sel]), np.log(th), r["xip"][n])
        num += np.sum(d[sel] * m / e[sel] ** 2)
        den += np.sum(m ** 2 / e[sel] ** 2)
        per[f"{i + 1}-{j + 1}"] = round(float(np.sum(d[sel] * m / e[sel] ** 2) / np.sum(m ** 2 / e[sel] ** 2)), 3)
        r.setdefault("_data", {})[n] = (x, d, e)
    if den == 0:
        raise ValueError(f"no overlapping bin pairs / theta between the prediction and {dataset} xi+")
    return {"dataset": dataset, "A": round(float(num / den), 4), "A_err": round(float(den ** -0.5), 4),
            "A_per_pair": per, "theta_range": [tmin, tmax]}


def _plot(r, th, nb, out, stem, has_b, comp):
    plt = mpl()
    pairs = r["pairs"]
    fig, axes = plt.subplots(nb, nb, figsize=(2.3 * nb + 1, 2.0 * nb + 1), sharex=True, squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for n, (i, j) in enumerate(pairs):
        ax = axes[j, i]
        ax.set_visible(True)
        ax.loglog(th, th * r["xip"][n] * 1e4, "C0", label=r"$\xi_+$")
        ax.loglog(th, th * np.abs(r["xim"][n]) * 1e4, "C1", label=r"$|\xi_-|$")
        if has_b:
            ax.loglog(th, th * r["xip_dmo"][n] * 1e4, "C0:", lw=1, label=r"$\xi_+$ DMO")
            ax.loglog(th, th * np.abs(r["xim_dmo"][n]) * 1e4, "C1:", lw=1)
        if comp and n in r.get("_data", {}):
            x, d, e = r["_data"][n]
            ax.errorbar(x, x * d * 1e4, x * e * 1e4, fmt="k.", ms=3, label=comp["dataset"])
        ax.text(0.05, 0.08, f"{i + 1}-{j + 1}", transform=ax.transAxes, fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel(r"$\theta$ [arcmin]")
    for row in axes:
        row[0].set_ylabel(r"$\theta\,\xi\times10^4$")
    axes[pairs[0][1], pairs[0][0]].legend(fontsize=6)
    return save_fig(fig, out, stem)
