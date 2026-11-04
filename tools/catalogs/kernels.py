"""Tier-1 shear-catalogue kernels (JSON in -> JSON out, files to cwd).

Imports only numpy + astropy (+ stdlib, certifi optional); ``shear_correlation``
also needs treecorr (declared pip_dep). Everything here runs unchanged on an
HPC node: the KiDS-1000 subset reader works from the public catalogue URL with
HTTP (multi-)Range requests, so nobody needs the 17.7 GB file.

KiDS-1000 gold WL catalogue layout (verified 2026-11-04 from its FITS header):
21,262,011 rows x 833 bytes, table data starting at byte 169920. Rows are NOT
grouped by tile: the file interleaves tiles in short runs (1 to a few hundred
rows), so a tile is ~100 byte ranges. ``data/kids1000_gold_tile_runs.npz``
(built by scripts/build_kids1000_tile_index.py from Range reads of THELI_NAME)
maps every row run to its tile; callers turn it into explicit row ranges.
"""

from __future__ import annotations

import re
import ssl
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

KIDS1000_URL = "https://kids.strw.leidenuniv.nl/DR4/data_files/KiDS_DR4.1_ugriZYJHKs_SOM_gold_WL_cat.fits"
KIDS1000_LAYOUT = {"header_bytes": 169920, "row_bytes": 833, "n_rows": 21262011}
# Asgari et al. 2021 (A&A 645, A104) Table 1: z_B bin edges (lo < z_B <= hi), multiplicative bias m +/- sigma_m.
KIDS1000_ZB_EDGES = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
KIDS1000_M_BIAS = [-0.009, -0.011, -0.015, 0.002, 0.007]
KIDS1000_M_SIGMA = [0.019, 0.020, 0.017, 0.012, 0.010]
# byte offsets inside one 833-byte row (from the TFORMs of the BINTABLE header) -> slim output columns
KIDS1000_COLUMNS = {  # name: (offset, big-endian dtype) - checked against fits_table_layout() 2026-11-04
    "ra": (89, ">f8"), "dec": (97, ">f8"), "z_b": (565, ">f8"), "z_b_min": (573, ">f8"),
    "z_b_max": (581, ">f8"), "mag_r": (395, ">f4"), "snr": (703, ">f4"),
    "psf_e1": (711, ">f4"), "psf_e2": (715, ">f4"),
    "e1": (821, ">f4"), "e2": (825, ">f4"), "weight": (829, ">f4"),
}
# catalogue column names behind the slim names (used by the offset self-check)
KIDS1000_SOURCE_NAMES = {"ra": "RAJ2000", "dec": "DECJ2000", "z_b": "Z_B", "z_b_min": "Z_B_MIN",
                         "z_b_max": "Z_B_MAX", "mag_r": "MAG_GAAP_r", "snr": "model_SNratio", "psf_e1": "PSF_e1",
                         "psf_e2": "PSF_e2", "e1": "e1", "e2": "e2", "weight": "weight"}
THELI_NAME = (669, 16)
_FMT = {"J": ">i4", "E": ">f4", "D": ">f8", "I": ">i2", "K": ">i8", "B": "u1", "L": "u1"}


# ------------------------------------------------------------------ HTTP range reads
def _ctx():
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:  # pragma: no cover - certifi absent
        return ssl.create_default_context()


def _get(url: str, rng: str, timeout: float = 300.0, retries: int = 5) -> tuple[bytes, str, int]:
    req = urllib.request.Request(url, headers={"User-Agent": "weak-lensing-mcp-server/0.2", "Range": "bytes=" + rng})
    err = None
    for k in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=_ctx()) as f:
                if f.status != 206:
                    raise RuntimeError(f"server ignored the Range header (HTTP {f.status})")
                return f.read(), f.headers.get("Content-Type", ""), f.status
        except Exception as exc:  # noqa: BLE001 - retry transient network errors
            err = exc
            time.sleep(1.5 * (k + 1))
    raise RuntimeError(f"Range read of {url} failed after {retries} tries: {err}")


def http_range(url: str, start: int, stop: int) -> bytes:
    """Bytes [start, stop) of ``url`` via one HTTP Range request."""
    data, _, _ = _get(url, f"{start}-{stop - 1}")
    if len(data) != stop - start:
        raise RuntimeError(f"short read {len(data)} of {stop - start} bytes")
    return data


def http_multirange(url: str, spans: list[tuple[int, int]]) -> dict[int, bytes]:
    """{start: bytes} for several [start, stop) spans in ONE request (multipart/byteranges).

    The KiDS server accepts up to ~200 ranges per request (1000 -> HTTP 400).
    """
    if len(spans) == 1:
        a, b = spans[0]
        return {a: http_range(url, a, b)}
    body, ctype, _ = _get(url, ",".join(f"{a}-{b - 1}" for a, b in spans))
    if "boundary=" not in ctype:
        raise RuntimeError(f"expected multipart/byteranges, got {ctype!r}")
    bnd = ctype.split("boundary=")[1].strip().strip('"').encode()
    out = {}
    for seg in body.split(b"--" + bnd):
        m = re.search(rb"Content-range: bytes (\d+)-(\d+)/", seg, re.I)
        if m:
            a, b = int(m.group(1)), int(m.group(2)) + 1
            payload = seg.split(b"\r\n\r\n", 1)[1]
            out[a] = payload[: b - a]
    if len(out) != len(spans):
        raise RuntimeError(f"multipart reply had {len(out)} parts for {len(spans)} ranges")
    return out


def fits_table_layout(url: str) -> dict:
    """Read the first BINTABLE header of a remote FITS file with Range reads.

    Returns header_bytes (data offset), row_bytes, n_rows and columns
    [(name, numpy dtype, byte offset)] - used to verify KIDS1000_COLUMNS.
    """
    from astropy.io import fits

    blob = http_range(url, 0, 2880 * 120)
    cards = [blob[i:i + 80].decode("ascii", "replace") for i in range(0, len(blob), 80)]
    ends = [i for i, c in enumerate(cards) if c.startswith("END ")]
    p_end = ((ends[0] + 1) * 80 + 2879) // 2880 * 2880
    ext = cards[p_end // 80:]
    e1 = next(i for i, c in enumerate(ext) if c.startswith("END "))
    h = fits.Header.fromstring("".join(ext[:e1 + 1]))
    cols, off = [], 0
    for i in range(1, h["TFIELDS"] + 1):
        form = h[f"TFORM{i}"].strip()
        n, code = int(form[:-1] or 1), form[-1]
        cols.append((h[f"TTYPE{i}"].strip(), f"S{n}" if code == "A" else _FMT[code], off))
        off += n if code == "A" else n * np.dtype(_FMT[code]).itemsize
    return {"header_bytes": p_end + ((e1 + 1) * 80 + 2879) // 2880 * 2880, "row_bytes": int(h["NAXIS1"]),
            "n_rows": int(h["NAXIS2"]), "columns": cols}


def _decode_rows(buf: bytes, row_bytes: int) -> dict:
    a = np.frombuffer(buf, dtype=np.uint8).reshape(-1, row_bytes)
    out = {}
    for name, (off, dt) in KIDS1000_COLUMNS.items():
        w = np.dtype(dt).itemsize
        out[name] = np.ascontiguousarray(a[:, off:off + w]).view(dt).ravel().astype(np.float64)
    return out


def read_row_ranges(row_ranges: list, url: str = KIDS1000_URL, layout: dict | None = None,
                    ranges_per_request: int = 150, workers: int = 8) -> dict:
    """Read rows given as [[row_start, n_rows], ...] and return the slim columns.

    Ranges are batched into multi-range requests and fetched in parallel.
    """
    lay = layout or KIDS1000_LAYOUT
    hb, rb = lay["header_bytes"], lay["row_bytes"]
    spans = [(hb + int(s) * rb, hb + (int(s) + int(n)) * rb) for s, n in row_ranges if int(n) > 0]
    batches = [spans[k:k + ranges_per_request] for k in range(0, len(spans), ranges_per_request)]
    with ThreadPoolExecutor(max(1, workers)) as ex:
        parts = list(ex.map(lambda b: http_multirange(url, b), batches))
    bufs = [p[a] for batch, p in zip(batches, parts, strict=True) for a, _ in batch]
    if not bufs:
        return {k: np.zeros(0) for k in KIDS1000_COLUMNS}
    return _decode_rows(b"".join(bufs), rb)


# ------------------------------------------------------------------ tile index -> row ranges
def tile_centre(name: str) -> tuple[float, float]:
    """KIDS_129p0_m1p5 -> (129.0, -1.5)."""
    m = re.match(r"KIDS_(\d+)p(\d+)_(m?)(\d+)p(\d+)", name)
    if not m:
        raise ValueError(f"not a KiDS tile name: {name!r}")
    return float(f"{m.group(1)}.{m.group(2)}"), float(f"{m.group(4)}.{m.group(5)}") * (-1 if m.group(3) else 1)


def select_tiles(tiles_all: list[str], tiles=None, ra_range=None, dec_range=None, patch=None,
                 max_tiles=None) -> list[str]:
    """Choose tiles by explicit names, RA/Dec box overlap (~1x1 deg tiles) and/or patch 'N'/'S'."""
    if tiles:
        missing = set(tiles) - set(tiles_all)
        if missing:
            raise ValueError(f"unknown KiDS-1000 tiles {sorted(missing)[:5]}; tile names look like KIDS_129p0_m1p5")
        sel = list(dict.fromkeys(tiles))
    else:
        sel = list(tiles_all)
    cen = {t: tile_centre(t) for t in sel}
    if patch:
        north = patch.upper().startswith("N")
        sel = [t for t in sel if (cen[t][1] > -10) == north]
    if ra_range:
        def ok_ra(t):
            ra, dec = cen[t]
            pad = 0.5 / max(np.cos(np.deg2rad(dec)), 0.2)  # tiles are ~1 x 1 deg: strict overlap test
            lo, hi = (float(ra_range[0]) - pad) % 360, (float(ra_range[1]) + pad) % 360
            return lo < ra < hi if lo <= hi else (ra > lo or ra < hi)  # box through RA = 0
        sel = [t for t in sel if ok_ra(t)]
    if dec_range:
        sel = [t for t in sel if dec_range[0] - 0.5 < cen[t][1] < dec_range[1] + 0.5]
    if not tiles:
        sel.sort(key=lambda t: (cen[t][1], cen[t][0]))
    return sel[:max_tiles] if max_tiles else sel


def tile_row_ranges(index: dict, tile: str) -> list[list[int]]:
    """[[row_start, n_rows], ...] of one tile from the run index {'tiles': [...], 'runs': (N,3)}."""
    tid = list(index["tiles"]).index(tile)
    r = index["runs"]
    sub = r[r[:, 0] == tid]
    return [[int(s), int(n)] for s, n in sub[:, 1:]]


# ------------------------------------------------------------------ subsets
def tomo_bin_of(z_b, edges=KIDS1000_ZB_EDGES) -> np.ndarray:
    """Tomographic bin 1..N for lo < z_B <= hi; 0 outside the edges."""
    z = np.asarray(z_b, float)
    out = np.zeros(z.shape, dtype=np.int16)
    for k in range(len(edges) - 1):
        out[(z > edges[k]) & (z <= edges[k + 1])] = k + 1
    return out


def footprint_area(ra, dec, pix_arcmin: float = 4.0) -> float:
    """Occupied area (deg^2): number of pix_arcmin cells (RA cos dec x Dec) holding >= 1 galaxy."""
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    if ra.size == 0:
        return 0.0
    ra = (ra - np.median(ra) + 180.0) % 360.0 - 180.0
    x = ra * np.cos(np.deg2rad(np.median(dec)))
    s = pix_arcmin / 60.0
    cells = np.unique(np.stack([np.floor(x / s), np.floor(dec / s)], axis=1), axis=0)
    return float(len(cells) * s * s)


def kids1000_subset(tile_ranges: dict, out_file: str, ra_range=None, dec_range=None, tomo_bins=None,
                    max_rows: int | None = None, cache_dir: str | None = None, url: str = KIDS1000_URL) -> dict:
    """Download a KiDS-1000 gold-catalogue subset and write a slim FITS table.

    ``tile_ranges``: {tile: [[row_start, n_rows], ...]} (from the bundled run
    index). Whole tiles are read (cached per tile as .npz in ``cache_dir``),
    then cut to the RA/Dec box, tomographic bins (Asgari+21 z_B edges) and
    max_rows; downloading stops once max_rows objects pass the cuts.
    """
    from astropy.table import Table

    cache = Path(cache_dir) if cache_dir else None
    if cache:
        cache.mkdir(parents=True, exist_ok=True)
    chunks, used, n_total, bytes_read = [], [], 0, 0
    for tile, ranges in tile_ranges.items():
        f = cache / f"{tile}.npz" if cache else None
        if f is not None and f.is_file():
            d = dict(np.load(f))
        else:
            d = read_row_ranges(ranges, url=url)
            bytes_read += sum(n for _, n in ranges) * KIDS1000_LAYOUT["row_bytes"]
            if f is not None:
                np.savez(f, **d)
        d["tomo_bin"] = tomo_bin_of(d["z_b"])
        keep = np.ones(len(d["ra"]), bool)
        if ra_range:
            lo, hi = float(ra_range[0]) % 360, float(ra_range[1]) % 360
            keep &= ((d["ra"] >= lo) & (d["ra"] <= hi)) if lo <= hi else ((d["ra"] >= lo) | (d["ra"] <= hi))
        if dec_range:
            keep &= (d["dec"] >= dec_range[0]) & (d["dec"] <= dec_range[1])
        if tomo_bins:
            keep &= np.isin(d["tomo_bin"], list(tomo_bins))
        d = {k: v[keep] for k, v in d.items()}
        d["tile"] = np.full(len(d["ra"]), tile)
        chunks.append(d)
        used.append(tile)
        n_total += len(d["ra"])
        if max_rows and n_total >= max_rows:
            break
    cat = {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]}
    if max_rows:
        cat = {k: v[:max_rows] for k, v in cat.items()}
    tab = Table(cat)
    tab.meta.update({"SURVEY": "KiDS-1000", "CATALOG": "KiDS_DR4.1_ugriZYJHKs_SOM_gold_WL_cat", "NTILES": len(used),
                     "MCORR": "none (raw lensfit e1,e2; apply m per tomo bin)"})
    tab.write(out_file, overwrite=True)
    tb = cat["tomo_bin"]
    return {"file": str(out_file), "n_rows": len(tab), "tiles": used, "n_tiles": len(used),
            "bytes_downloaded": int(bytes_read), "columns": list(tab.colnames),
            "area_deg2": round(footprint_area(cat["ra"], cat["dec"]), 3),
            "tomo_counts": {str(int(b)): int((tb == b).sum()) for b in np.unique(tb)}}


# ------------------------------------------------------------------ tier 1 -> tier 2: xi+/xi-
def load_catalog(path: str) -> dict:
    """Slim catalogue (FITS from fetch_shear_catalog) -> dict of numpy arrays."""
    from astropy.table import Table
    t = Table.read(path)
    out = {c: np.asarray(t[c]) for c in t.colnames}
    if "tile" in out and out["tile"].dtype.kind == "S":
        out["tile"] = out["tile"].astype(str)
    return out


def shear_correlation_arrays(cat: dict, bin_pairs: list, m_bias: dict | None = None, min_sep: float = 0.5,
                             max_sep: float = 300.0, nbins: int = 9, bin_slop: float = 0.05,
                             remove_c: bool = True, num_threads: int | None = None, npatch: int = 0) -> list[dict]:
    """xi+/xi- per tomographic bin pair with weights, m-correction and c-term removal (treecorr).

    Per bin: c_k = weighted mean of e_k (removed when remove_c); the
    multiplicative bias enters as xi -> xi / ((1+m_i)(1+m_j)) - KiDS-1000
    practice for a per-bin scalar m. Errors are treecorr shape-noise
    (varxip) or jackknife over ``npatch`` k-means patches when npatch >= 2.
    """
    import treecorr

    m_bias = m_bias or {}
    tb = np.asarray(cat["tomo_bin"]).astype(int)
    cats, info = {}, {}
    need = sorted({int(b) for p in bin_pairs for b in p})
    patch_centers = None
    for b in need:
        s = tb == b
        if s.sum() < 10:
            raise ValueError(f"tomographic bin {b} has only {int(s.sum())} galaxies in this catalogue")
        w = np.asarray(cat["weight"], float)[s]
        e1, e2 = np.asarray(cat["e1"], float)[s], np.asarray(cat["e2"], float)[s]
        c1, c2 = (float(np.average(e1, weights=w)), float(np.average(e2, weights=w))) if remove_c else (0.0, 0.0)
        kw = dict(ra=np.asarray(cat["ra"], float)[s], dec=np.asarray(cat["dec"], float)[s], g1=e1 - c1, g2=e2 - c2,
                  w=w, ra_units="deg", dec_units="deg")
        if npatch >= 2:
            if patch_centers is None:
                c0 = treecorr.Catalog(**kw, npatch=npatch)
                patch_centers = c0.patch_centers
            kw["patch_centers"] = patch_centers
        cats[b] = treecorr.Catalog(**kw)
        neff = float(w.sum() ** 2 / (w ** 2).sum())
        info[b] = {"n": int(s.sum()), "n_eff_weighted": round(neff, 1), "c1": c1, "c2": c2,
                   "m": float(m_bias.get(b, m_bias.get(str(b), 0.0))),
                   "sigma_e_per_comp": float(np.sqrt(np.sum(w ** 2 * (e1 ** 2 + e2 ** 2)) / (2 * np.sum(w ** 2))))}
    rows = []
    for i, j in bin_pairs:
        i, j = int(i), int(j)
        gg = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, sep_units="arcmin",
                                    bin_slop=bin_slop, var_method="jackknife" if npatch >= 2 else "shot",
                                    num_threads=num_threads)
        if i == j:
            gg.process(cats[i])
        else:
            gg.process(cats[i], cats[j])
        corr = 1.0 / ((1.0 + info[i]["m"]) * (1.0 + info[j]["m"]))
        if npatch >= 2:
            n = len(gg.xip)
            cov = gg.cov
            ep, em = np.sqrt(np.diag(cov)[:n]), np.sqrt(np.diag(cov)[n:])
        else:
            ep, em = np.sqrt(gg.varxip), np.sqrt(gg.varxim)
        for k in range(len(gg.meanr)):
            rows.append({"bin1": i, "bin2": j, "angbin": k, "theta_arcmin": float(gg.meanr[k]),
                         "theta_nominal": float(gg.rnom[k]), "xip": float(gg.xip[k] * corr),
                         "xim": float(gg.xim[k] * corr), "xip_err": float(ep[k] * corr),
                         "xim_err": float(em[k] * corr), "xip_im": float(gg.xip_im[k] * corr),
                         "npairs": float(gg.npairs[k]), "weight": float(gg.weight[k])})
    return [{"bins": {str(k): v for k, v in info.items()}}] + rows


def shear_correlation(catalog_file: str, bin_pairs: list, m_bias: dict | None = None, min_sep: float = 0.5,
                      max_sep: float = 300.0, nbins: int = 9, bin_slop: float = 0.05, remove_c: bool = True,
                      npatch: int = 0, out_csv: str = "xipm_measured.csv", subset: dict | None = None) -> dict:
    """Dispatch entry point: catalogue FITS -> xi+/xi- CSV in cwd + JSON summary.

    If ``catalog_file`` is not present in the job cwd and ``subset`` is given
    (kids1000_subset kwargs: tile_ranges, ra_range, dec_range, tomo_bins,
    max_rows), the identical subset is first rebuilt on the node from the
    public KiDS-1000 file - so nothing has to be staged in.
    """
    import csv

    if not Path(catalog_file).is_file():
        if not subset:
            raise FileNotFoundError(f"{catalog_file} not in the job directory and no subset spec given")
        kids1000_subset(out_file=catalog_file, **subset)
    res = shear_correlation_arrays(load_catalog(catalog_file), bin_pairs, m_bias, min_sep, max_sep, nbins,
                                   bin_slop, remove_c, None, npatch)
    info, rows = res[0], res[1:]
    with open(out_csv, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    return {"csv": str(out_csv), "n_rows": len(rows), **info, "rows": rows}


# ------------------------------------------------------------------ tier 1 -> tier 2: n(z) histogram
def nz_histogram(z: np.ndarray, w: np.ndarray, tomo: np.ndarray, bins: list[int], z_edges: np.ndarray) -> dict:
    """Weighted, unit-normalised histogram of a redshift column per tomographic bin."""
    zc = 0.5 * (z_edges[1:] + z_edges[:-1])
    out = {"z": zc.tolist(), "bins": {}}
    for b in bins:
        s = tomo == b
        h, _ = np.histogram(z[s], bins=z_edges, weights=w[s])
        norm = h.sum() * np.diff(z_edges)[0]
        nz = h / norm if norm > 0 else h
        out["bins"][str(b)] = {"nz": nz.tolist(), "n": int(s.sum()),
                               "mean_z": float(np.average(z[s], weights=w[s])) if s.any() else float("nan")}
    return out
