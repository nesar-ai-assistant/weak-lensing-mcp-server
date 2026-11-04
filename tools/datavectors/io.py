"""Readers and downloaders for weak-lensing 2pt data vectors (pure python).

Formats:
- CosmoSIS "2pt" FITS (DES Y3, KiDS-1000, KiDS-Legacy, HSC Y3 real space):
  a COVMAT image with NAME_i/STRT_i header keys giving block order, one
  binary table per statistic (BIN1, BIN2, ANGBIN, VALUE, ANG) and NZ_* tables
  (Z_LOW, Z_MID, Z_HIGH, BIN1..BINn).
- SACC (HSC Y3 Fourier space C_ell), read with the optional ``sacc`` package.

KiDS products ship as one tarball per release; members are extracted on
demand into the cache. kids.strw.leidenuniv.nl serves a certificate chain
that Python's default CA store on some systems cannot verify, so downloads
go through a certifi-backed opener (falls back to the default if certifi is
missing).
"""

from __future__ import annotations

import contextlib
import ssl
import tarfile
import urllib.request
from pathlib import Path

import numpy as np

from ..common import data_cache_dir, fetch_file, get_cached
from ..registry import DATASETS

# statistic -> (file key in DATASETS[ds]["files"], FITS extension or SACC data type)
STATISTICS: dict[str, dict[str, tuple[str, str]]] = {
    "des_y3": {"xip": ("datavector", "xip"), "xim": ("datavector", "xim"),
               "gammat": ("datavector", "gammat"), "wtheta": ("datavector", "wtheta")},
    "kids_1000": {"xip": ("xipm", "xiP"), "xim": ("xipm", "xiM"),
                  "cosebis": ("cosebis", "En"), "bandpowers": ("bandpowers", "PeeE")},
    "kids_legacy": {"xip": ("xipm", "xiP"), "xim": ("xipm", "xiM"),
                    "cosebis": ("cosebis", "En"), "cosebis_b": ("cosebis", "Bn"),
                    "bandpowers": ("bandpowers", "PeeE"), "bandpowers_b": ("bandpowers", "PeeB")},
    "hsc_y3": {"xip": ("real", "xi_plus"), "xim": ("real", "xi_minus"),
               "cl_ee": ("fourier", "cl_ee")},
}

# x-axis meaning for each extension
X_LABEL = {"xip": "theta_arcmin", "xim": "theta_arcmin", "xiP": "theta_arcmin", "xiM": "theta_arcmin",
           "xi_plus": "theta_arcmin", "xi_minus": "theta_arcmin", "gammat": "theta_arcmin",
           "wtheta": "theta_arcmin", "En": "n_mode", "Bn": "n_mode", "PeeE": "ell", "PeeB": "ell",
           "cl_ee": "ell"}


def check_dataset(dataset: str) -> dict:
    if dataset not in DATASETS:
        raise ValueError(f"Unknown dataset {dataset!r}. Available: {', '.join(DATASETS)}")
    meta = DATASETS[dataset]
    if meta["access"] == "unreleased":
        raise RuntimeError(f"{dataset}: {meta['note']}")
    return meta


@contextlib.contextmanager
def _certifi_opener():
    """Temporarily install a urllib opener that verifies against certifi's CA bundle."""
    try:
        import certifi
        ctx = ssl.create_default_context(cafile=certifi.where())
    except Exception:  # pragma: no cover - certifi absent
        yield
        return
    old = urllib.request._opener  # noqa: SLF001 - restore whatever was installed
    urllib.request.install_opener(urllib.request.build_opener(urllib.request.HTTPSHandler(context=ctx)))
    try:
        yield
    finally:
        urllib.request.install_opener(old)


def download(url: str, subdir: str, timeout: float = 600.0) -> Path:
    """fetch_file with the certifi opener and a manual-download hint on failure."""
    with _certifi_opener():
        try:
            return fetch_file(url, subdir=subdir, timeout=timeout)
        except RuntimeError as exc:
            dest = data_cache_dir() / subdir / url.rstrip("/").split("/")[-1]
            raise RuntimeError(f"{exc}\nManual fallback: download {url} in a browser and save it as {dest}") from exc


def archive_member(dataset: str, relpath: str) -> Path:
    """Path of ``relpath`` inside a KiDS release tarball, extracting it on first use."""
    meta = DATASETS[dataset]
    out_root = data_cache_dir() / dataset / "extracted"
    hits = list(out_root.glob(f"*/{relpath}")) if out_root.is_dir() else []
    if hits:
        return hits[0]
    tgz = download(meta["archive"], subdir=dataset)
    with tarfile.open(tgz) as tf:
        members = [m for m in tf.getmembers() if m.isfile() and m.name.endswith("/" + relpath)]
        if not members:
            raise FileNotFoundError(f"{relpath} not found in {tgz.name}")
        tf.extractall(out_root, members=members, filter="data")
    return out_root / members[0].name


def product_path(dataset: str, file_key: str) -> Path:
    """Local path of a data-vector file (download / extract as needed)."""
    meta = check_dataset(dataset)
    rel = meta["files"][file_key]
    if "archive" in meta:
        return archive_member(dataset, rel)
    return download(rel, subdir=dataset)


# ---------------------------------------------------------------- 2pt FITS
def read_2pt_fits(path: Path) -> dict:
    """Parse a CosmoSIS 2pt FITS file -> {"blocks": {ext: table dict}, "cov", "order", "nz": {...}}."""
    from astropy.io import fits

    def _load():
        out = {"blocks": {}, "nz": {}, "order": []}
        with fits.open(path, memmap=False) as h:
            names = [hdu.name for hdu in h]
            cov_hdu = h["COVMAT"]
            out["cov"] = np.asarray(cov_hdu.data, dtype=float)
            hd = cov_hdu.header
            i = 0
            while f"NAME_{i}" in hd:
                out["order"].append((hd[f"NAME_{i}"], int(hd[f"STRT_{i}"])))
                i += 1
            for nm in names:
                hdu = h[nm]
                if nm.upper().startswith("NZ_") and "REALISATION" not in nm.upper():
                    cols = [c for c in hdu.columns.names if c.startswith("BIN")]
                    out["nz"][nm.lower()] = {"z": np.asarray(hdu.data["Z_MID"], float),
                                             "bins": [np.asarray(hdu.data[c], float) for c in cols]}
            for name, _ in out["order"]:
                d = h[name].data
                out["blocks"][name] = {
                    "bin1": np.asarray(d["BIN1"], int), "bin2": np.asarray(d["BIN2"], int),
                    "angbin": np.asarray(d["ANGBIN"], int), "value": np.asarray(d["VALUE"], float),
                    "x": np.asarray(d["ANG"], float)}
        return out

    return get_cached(f"2pt:{path}", _load)


def read_sacc(path: Path) -> dict:
    """Parse an HSC-style shear SACC into the same structure as read_2pt_fits."""
    try:
        import sacc
    except ImportError as exc:  # pragma: no cover
        raise ImportError("The HSC Fourier data vector needs the 'sacc' package: pip install sacc") from exc

    def _load():
        s = sacc.Sacc.load_fits(str(path))
        dt = "cl_ee"
        b1, b2, ab, val, x = [], [], [], [], []
        for t1, t2 in s.get_tracer_combinations(dt):
            ell, cl = s.get_ell_cl(dt, t1, t2)
            i1, i2 = int(t1.split("_")[-1]) + 1, int(t2.split("_")[-1]) + 1
            b1 += [i1] * len(ell)
            b2 += [i2] * len(ell)
            ab += list(range(1, len(ell) + 1))
            val += list(cl)
            x += list(ell)
        idx = np.concatenate([s.indices(dt, (t1, t2)) for t1, t2 in s.get_tracer_combinations(dt)])
        cov = np.asarray(s.covariance.covmat)[np.ix_(idx, idx)]
        names = sorted(s.tracers)
        nz = {"nz_source": {"z": np.asarray(s.tracers[names[0]].z, float),
                            "bins": [np.asarray(s.tracers[n].nz, float) for n in names]}}
        block = {"bin1": np.array(b1), "bin2": np.array(b2), "angbin": np.array(ab),
                 "value": np.array(val), "x": np.array(x)}
        return {"blocks": {dt: block}, "cov": cov, "order": [(dt, 0)], "nz": nz}

    return get_cached(f"sacc:{path}", _load)


def load_product(dataset: str, file_key: str) -> dict:
    path = product_path(dataset, file_key)
    return read_sacc(path) if path.suffix == ".sacc" else read_2pt_fits(path)


def block_slice(dv: dict, ext: str) -> slice:
    """Index range of extension ``ext`` in the full covariance."""
    starts = dict(dv["order"])
    if ext not in starts:
        raise KeyError(f"{ext} not in data vector (has {list(starts)})")
    n = len(dv["blocks"][ext]["value"])
    return slice(starts[ext], starts[ext] + n)


def load_statistic(dataset: str, statistic: str) -> tuple[dict, np.ndarray, str]:
    """(block table, its covariance sub-matrix, extension name) for one statistic."""
    check_dataset(dataset)
    stats = STATISTICS[dataset]
    if statistic not in stats:
        raise ValueError(f"{dataset} has statistics {list(stats)}; got {statistic!r}")
    key, ext = stats[statistic]
    dv = load_product(dataset, key)
    sl = block_slice(dv, ext)
    return dv["blocks"][ext], dv["cov"][sl, sl], ext


def load_nz(dataset: str, sample: str = "source") -> dict:
    """n(z) per tomographic bin: {"z": array, "bins": [array, ...]}."""
    check_dataset(dataset)
    if dataset == "hsc_y3":
        if sample != "source":
            raise ValueError("HSC Y3 cosmic shear has source n(z) only")
        return load_product(dataset, "fourier")["nz"]["nz_source"]
    key = next(iter(DATASETS[dataset]["files"]))
    if dataset.startswith("kids"):
        key = "xipm"
    nz = load_product(dataset, key)["nz"]
    name = f"nz_{sample}"
    if name not in nz:
        raise ValueError(f"{dataset} provides {list(nz)}; no {name}")
    return nz[name]
