"""Rebuild data/kids1000_gold_tile_runs.npz: which rows of the 17.7 GB KiDS-1000 gold WL catalogue belong
to which tile, using only HTTP multi-range reads of the 16-byte THELI_NAME column (~3 min at 16 threads).
Rows are interleaved across tiles in short runs, so the index stores every run (tile_id, row_start,
n_rows); together they cover all 21,262,011 rows exactly once.

Usage: python scripts/build_kids1000_tile_index.py [sample_step=64]
"""
import re
import ssl
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import certifi
import numpy as np

ctx = ssl.create_default_context(cafile=certifi.where())
U = "https://kids.strw.leidenuniv.nl/DR4/data_files/KiDS_DR4.1_ugriZYJHKs_SOM_gold_WL_cat.fits"
H, L, N, NO = 169920, 833, 21262011, 669
STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 64


def fetch_chunk(chunk):
    rs = ",".join(f"{H + r * L + NO}-{H + r * L + NO + 15}" for r in chunk)
    req = urllib.request.Request(U, headers={"User-Agent": "wl-mcp-index", "Range": "bytes=" + rs})
    for k in range(6):
        try:
            with urllib.request.urlopen(req, timeout=120, context=ctx) as f:
                body, ctype = f.read(), f.headers.get("Content-Type", "")
            break
        except Exception as e:
            print("retry", k, e, flush=True)
            time.sleep(2 + 3 * k)
    else:
        raise RuntimeError("failed")
    out = {}
    if len(chunk) == 1:
        return {chunk[0]: body[:16].decode().strip("\x00 ")}
    bnd = ctype.split("boundary=")[1].encode()
    for seg in body.split(b"--" + bnd):
        m = re.search(rb"Content-range: bytes (\d+)-(\d+)/", seg, re.I)
        if m:
            out[(int(m.group(1)) - H - NO) // L] = seg.split(b"\r\n\r\n", 1)[1][:16].decode().strip("\x00 ")
    assert len(out) == len(chunk), (len(out), len(chunk))
    return out


def names(rows, pool):
    rows = sorted(set(int(r) for r in rows))
    chunks = [rows[k:k + 200] for k in range(0, len(rows), 200)]
    res = {}
    for d in pool.map(fetch_chunk, chunks):
        res.update(d)
    return res


t0 = time.time()
with ThreadPoolExecutor(16) as pool:
    samp = list(range(0, N, STEP)) + [N - 1]
    known = names(samp, pool)
    print("sampled", len(known), round(time.time() - t0), "s", flush=True)
    pending = [(a, b) for a, b in zip(samp[:-1], samp[1:], strict=True) if known[a] != known[b]]
    bounds = []
    it = 0
    while pending:
        it += 1
        need = [(a + b) // 2 for a, b in pending if b - a > 1 and (a + b) // 2 not in known]
        known.update(names(need, pool))
        nxt = []
        for a, b in pending:
            if b - a <= 1:
                bounds.append(b)
                continue
            m = (a + b) // 2
            if known[m] != known[a]:
                nxt.append((a, m))
            if known[m] != known[b]:
                nxt.append((m, b))
        pending = nxt
        print("iter", it, "pending", len(pending), "bounds", len(bounds), round(time.time() - t0), "s", flush=True)

bounds = sorted(set(bounds))
starts = np.array([0] + bounds, dtype=np.int64)
ends = np.array(bounds + [N], dtype=np.int64)
run_names = [known[s] for s in starts]
tiles = sorted(set(run_names))
tid = {t: i for i, t in enumerate(tiles)}
runs = np.stack([np.array([tid[n] for n in run_names]), starts, ends - starts], axis=1).astype(np.int64)
print("runs", len(runs), "tiles", len(tiles), "rows", int(runs[:, 2].sum()))
out = __import__("pathlib").Path(__file__).resolve().parents[1] / "data" / "kids1000_gold_tile_runs.npz"
assert runs[0, 1] == 0 and (runs[1:, 1] == runs[:-1, 1] + runs[:-1, 2]).all() and runs[:, 2].sum() == N
np.savez_compressed(out, tiles=np.array(tiles), runs=runs.astype(np.int32))
print("done", round(time.time() - t0), "s")
