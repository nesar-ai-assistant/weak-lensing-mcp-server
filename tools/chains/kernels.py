"""Dispatchable chain-analysis kernels (JSON in -> JSON out).

Imports only numpy (+ stdlib; certifi optional). ``chain_marginals`` takes a
download URL (and, for tarballs, the member path), so it runs unchanged on
an HPC node; ``summarize_columns`` is the pure core shared with the local
branch of analyze_wl_chain.
"""

from __future__ import annotations

import tarfile
import urllib.request
from pathlib import Path

import numpy as np

QUANTILES = (0.025, 0.16, 0.5, 0.84, 0.975)


def parse_cosmosis_text(text_or_path) -> tuple[list[str], np.ndarray]:
    """CosmoSIS text chain: first line '#col1<TAB>col2...', '#' comments, whitespace data."""
    p = Path(text_or_path)
    with open(p) as fh:
        header = fh.readline().lstrip("#").strip()
    cols = header.split("\t") if "\t" in header else header.split()
    data = np.loadtxt(p, comments="#", ndmin=2)
    if data.shape[1] != len(cols):
        raise ValueError(f"{p.name}: {len(cols)} header columns but {data.shape[1]} data columns")
    return [c.lower() for c in cols], data


def _find(cols: list[str], *cands: str) -> int | None:
    for c in cands:
        if c in cols:
            return cols.index(c)
    return None


def standard_columns(cols: list[str], data: np.ndarray) -> dict[str, np.ndarray]:
    """Extract weight, Om, sigma8, S8 (derived when missing/NaN), h, w from a CosmoSIS chain."""
    iw = _find(cols, "weight")
    if iw is not None:
        w = data[:, iw].astype(float)
    else:
        il = _find(cols, "log_weight")
        if il is None:
            w = np.ones(len(data))
        else:
            lw = data[:, il]
            w = np.exp(lw - np.nanmax(lw[np.isfinite(lw)]))
    om = data[:, _find(cols, "cosmological_parameters--omega_m")]
    i8 = _find(cols, "cosmological_parameters--sigma_8")
    if i8 is None:
        raise ValueError("chain has no sigma_8 column")
    s8 = data[:, i8]
    iS = _find(cols, "cosmological_parameters--s_8")
    S8 = data[:, iS] if iS is not None else np.full(len(data), np.nan)
    derived = s8 * np.sqrt(om / 0.3)
    S8 = np.where(np.isfinite(S8), S8, derived)
    out = {"weight": w, "omega_m": om, "sigma_8": s8, "S_8": S8}
    ih = _find(cols, "cosmological_parameters--h0")
    if ih is not None:
        out["h0"] = data[:, ih]
    iwz = _find(cols, "cosmological_parameters--w")
    if iwz is not None:
        out["w"] = data[:, iwz]
    ip = _find(cols, "post")
    if ip is not None:
        out["post"] = data[:, ip]
    good = np.isfinite(w) & (w > 0) & np.isfinite(om) & np.isfinite(s8)
    return {k: v[good] for k, v in out.items()}


def _wquantiles(x, w, qs):
    i = np.argsort(x)
    c = np.cumsum(w[i])
    c = (c - 0.5 * w[i]) / c[-1]
    return np.interp(qs, c, x[i])


def _mode_hpd(x, w, mass=0.68, nbins=200):
    """Marginal mode and highest-posterior-density interval from a smoothed histogram."""
    lo, hi = _wquantiles(x, w, [0.0005, 0.9995])
    h, e = np.histogram(x, bins=nbins, range=(lo, hi), weights=w)
    k = np.exp(-0.5 * (np.arange(-8, 9) / 3.0) ** 2)
    h = np.convolve(h, k / k.sum(), mode="same")
    c = 0.5 * (e[1:] + e[:-1])
    order = np.argsort(h)[::-1]
    cum = np.cumsum(h[order]) / h.sum()
    sel = order[: np.searchsorted(cum, mass) + 1]
    return float(c[np.argmax(h)]), float(c[sel].min()), float(c[sel].max())


def summarize_columns(cols: dict[str, np.ndarray], params: list[str] | None = None) -> dict:
    """Weighted marginal summaries for each parameter (mean, std, quantiles, mode + 68% HPD)."""
    w = np.asarray(cols["weight"], float)
    neff = float(w.sum() ** 2 / np.sum(w ** 2))
    res = {"n_samples": int(len(w)), "n_eff": round(neff, 1), "params": {}}
    for p in params or [k for k in ("S_8", "omega_m", "sigma_8", "h0", "w") if k in cols]:
        x = np.asarray(cols[p], float)
        m = np.isfinite(x)
        x, ww = x[m], w[m]
        mean = float(np.sum(ww * x) / ww.sum())
        std = float(np.sqrt(np.sum(ww * (x - mean) ** 2) / ww.sum()))
        q = _wquantiles(x, ww, QUANTILES)
        mode, hlo, hhi = _mode_hpd(x, ww)
        res["params"][p] = {"mean": round(mean, 4), "std": round(std, 4), "median": round(float(q[2]), 4),
                            "lo68": round(float(q[1]), 4), "hi68": round(float(q[3]), 4),
                            "lo95": round(float(q[0]), 4), "hi95": round(float(q[4]), 4),
                            "mode": round(mode, 4), "hpd68": [round(hlo, 4), round(hhi, 4)]}
    if "post" in cols:
        res["max_post_sample"] = {k: round(float(cols[k][np.argmax(cols["post"])]), 4)
                                  for k in ("S_8", "omega_m", "sigma_8")}
    om, S8 = np.asarray(cols["omega_m"]), np.asarray(cols["S_8"])
    m = np.isfinite(om) & np.isfinite(S8)
    c = np.cov(om[m], S8[m], aweights=w[m])
    res["corr_S8_Om"] = round(float(c[0, 1] / np.sqrt(c[0, 0] * c[1, 1])), 3)
    return res


def _get(url: str, dest: Path) -> Path:
    if dest.is_file() and dest.stat().st_size > 0:
        return dest
    try:
        import ssl

        import certifi
        ctx = ssl.create_default_context(cafile=certifi.where())
    except Exception:
        ctx = None
    req = urllib.request.Request(url, headers={"User-Agent": "cosmo-mcp-server/0.1"})
    with urllib.request.urlopen(req, timeout=900, context=ctx) as r, open(dest, "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    return dest


def chain_marginals(url: str, member: str | None = None, params: list[str] | None = None) -> dict:
    """Download a CosmoSIS chain (or a member of a release tarball) to cwd and summarise it.

    Args (JSON): url - chain .txt URL or tarball URL; member - path suffix inside
    the tarball (None for a plain chain); params - subset of S_8, omega_m,
    sigma_8, h0, w. Returns summarize_columns(...) plus the source.
    """
    src = _get(url, Path.cwd() / url.rstrip("/").split("/")[-1])
    if member:
        with tarfile.open(src) as tf:
            m = next(x for x in tf.getmembers() if x.isfile() and x.name.endswith(member))
            text = tf.extractfile(m).read().decode()
        path = Path.cwd() / Path(member).name
        path.write_text(text)
    else:
        path = src
    cols, data = parse_cosmosis_text(path)
    res = summarize_columns(standard_columns(cols, data), params)
    res["source"] = {"url": url, "member": member}
    return res

