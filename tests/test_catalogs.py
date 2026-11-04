"""Shear-catalogue tests. Offline tests use a synthetic catalogue with a known shear signal;
network tests read real KiDS-1000 gold-catalogue rows with HTTP Range requests."""

import numpy as np
import pytest
from astropy.table import Table

from tools.catalogs import (
    _tile_index,
    fetch_shear_catalog,
    kernels,
    measure_nz_from_catalog,
    measure_shear_correlation,
    plot_shear_catalog,
)
from tools.registry import DATASETS


def _synthetic_catalog(path, n=60000, seed=1, g=(0.0, 0.0), c=(0.0, 0.0), m=0.0):
    """Random 2x2 deg patch; e = (1+m) g + c + shape noise; tomographic bins 3 (z_B=0.6) and 4 (0.8)."""
    rng = np.random.default_rng(seed)
    ra = rng.uniform(180, 182, n)
    dec = rng.uniform(-1, 1, n)
    zb = rng.choice([0.6, 0.8], n)
    e1 = (1 + m) * g[0] + c[0] + rng.normal(0, 0.27, n)
    e2 = (1 + m) * g[1] + c[1] + rng.normal(0, 0.27, n)
    t = Table({"ra": ra, "dec": dec, "e1": e1, "e2": e2, "weight": rng.uniform(1, 15, n), "z_b": zb,
               "tomo_bin": kernels.tomo_bin_of(zb), "tile": np.full(n, "KIDS_181p0_0p5")})
    t.write(path, overwrite=True)
    return str(path)


def test_tile_index_covers_catalogue():
    idx = _tile_index()
    runs = idx["runs"]
    sc = DATASETS["kids_1000"]["shear_catalog"]
    assert runs[0, 1] == 0 and np.all(runs[1:, 1] == runs[:-1, 1] + runs[:-1, 2])
    assert int(runs[:, 2].sum()) == sc["n_rows"] == kernels.KIDS1000_LAYOUT["n_rows"]
    assert 980 <= len(idx["tiles"]) <= 1006
    assert kernels.tile_centre("KIDS_129p0_m1p5") == (129.0, -1.5)


def test_select_tiles_and_tomo_bins():
    idx = _tile_index()
    sel = kernels.select_tiles(idx["tiles"], ra_range=[179.5, 182.5], dec_range=[-1, 1])
    assert sel and all(178.9 <= kernels.tile_centre(t)[0] <= 183.1 for t in sel)
    south = kernels.select_tiles(idx["tiles"], patch="S", max_tiles=5)
    assert len(south) == 5 and all(kernels.tile_centre(t)[1] < -20 for t in south)
    wrap = kernels.select_tiles(idx["tiles"], ra_range=[355, 5], patch="S")
    assert wrap and all(kernels.tile_centre(t)[0] > 350 or kernels.tile_centre(t)[0] < 10 for t in wrap)
    assert list(kernels.tomo_bin_of([0.1, 0.2, 0.3, 0.31, 0.9, 1.2, 1.3])) == [0, 1, 1, 2, 4, 5, 0]


def test_other_surveys_raise_with_reason(tmp_path):
    for ds in ("des_y3", "hsc_y3"):
        with pytest.raises(RuntimeError, match="not available"):
            fetch_shear_catalog(str(tmp_path), dataset=ds)


def test_correlation_removes_c_and_corrects_m(tmp_path):
    """Constant shear g on a patch: xi+ = |g|^2 at all theta once c is removed and m divided out."""
    g, c, m = (0.02, -0.01), (0.003, 0.002), -0.05
    path = _synthetic_catalog(tmp_path / "cat.fits", g=g, c=c, m=m)
    payload = kernels.shear_correlation_arrays(kernels.load_catalog(path), [[3, 3]], m_bias={"3": m},
                                               min_sep=2, max_sep=60, nbins=4, remove_c=False)
    xip = np.array([r["xip"] for r in payload[1:]])
    expect = ((1 + m) * np.array(g) + np.array(c)) ** 2 / (1 + m) ** 2
    assert np.allclose(xip, expect.sum(), rtol=0.15)
    # with c removal (constant g is itself removed - c-term and constant shear are degenerate)
    pay2 = kernels.shear_correlation_arrays(kernels.load_catalog(path), [[3, 3]], min_sep=2, max_sep=60, nbins=4)
    assert np.all(np.abs([r["xip"] for r in pay2[1:]]) < 5e-5)
    assert abs(pay2[0]["bins"]["3"]["c1"] - ((1 + m) * g[0] + c[0])) < 3e-3


def test_catalog_tools_offline(tmp_path):
    path = _synthetic_catalog(tmp_path / "cat.fits")
    q = plot_shear_catalog(path, str(tmp_path))
    assert q.metadata["data_kind"] == "catalog" and abs(q.metadata["sigma_e"] - 0.27) < 0.01
    r = measure_shear_correlation(path, str(tmp_path), min_sep=2, max_sep=60, nbins=4)
    assert r.metadata["data_kind"] == "catalog" and r.metadata["bin_pairs"] == [[3, 3], [3, 4], [4, 4]]
    assert measure_shear_correlation.weight == "heavy"
    near = r.metadata["xip_near_10arcmin"]["33"]
    assert abs(near["xip"]) < 5 * near["xip_err"]          # pure noise -> consistent with 0
    nz = measure_nz_from_catalog(path, str(tmp_path), compare_release=False)
    assert nz.metadata["mean_z"] == {"3": 0.6, "4": 0.8}


def test_measure_dispatch_branch(tmp_path, monkeypatch):
    import mcp_server.dispatch as d

    path = _synthetic_catalog(tmp_path / "cat.fits", n=3000)
    seen = {}

    def fake_run_kernel(function, args, pip_deps=None, duration=600, **kw):
        seen.update(function=function, args=args, pip_deps=pip_deps)
        rows = [{"bin1": 3, "bin2": 3, "angbin": 0, "theta_arcmin": 10.0, "theta_nominal": 10.0, "xip": 1e-4,
                 "xim": 1e-5, "xip_err": 1e-5, "xim_err": 1e-5, "xip_im": 0.0, "npairs": 10.0, "weight": 1.0}]
        return {"result": {"bins": {"3": {"n": 3000}}, "rows": rows}, "host": "polaris"}

    monkeypatch.setattr(d, "remote_site", lambda: "polaris")
    monkeypatch.setattr(d, "run_kernel", fake_run_kernel)
    r = measure_shear_correlation(path, str(tmp_path), bin_pairs=[[3, 3]])
    assert r.metadata["computed_on"] == "polaris"
    assert seen["function"] == "catalogs.kernels.shear_correlation"
    assert "treecorr>=5.0.1" in seen["pip_deps"] and seen["args"]["catalog_file"] == "cat.fits"
    sub = seen["args"]["subset"]   # node rebuilds the same subset: tiles, bounds, bins, row count
    assert list(sub["tile_ranges"]) == ["KIDS_181p0_0p5"] and sub["max_rows"] == 3000
    assert sub["tomo_bins"] == [3, 4] and 180 <= sub["ra_range"][0] < sub["ra_range"][1] <= 182


@pytest.mark.network
def test_kids1000_range_read_matches_header():
    lay = kernels.fits_table_layout(kernels.KIDS1000_URL)
    assert {k: lay[k] for k in kernels.KIDS1000_LAYOUT} == kernels.KIDS1000_LAYOUT
    cols = {c[0]: (c[2], c[1]) for c in lay["columns"]}
    for slim, src in kernels.KIDS1000_SOURCE_NAMES.items():
        assert cols[src] == kernels.KIDS1000_COLUMNS[slim], slim
    d = kernels.read_row_ranges([[0, 5], [1_000_000, 3]])
    assert np.allclose(d["ra"][:2], [128.614, 128.6121], atol=1e-3)
    assert np.all((d["z_b"] > 0) & (d["z_b"] < 3)) and np.all(d["weight"] > 0)


@pytest.mark.network
def test_fetch_one_kids1000_tile(tmp_path):
    from tools.catalogs import _subset_spec

    r = fetch_shear_catalog(str(tmp_path), tiles=["KIDS_181p0_0p5"], tomo_bins=[3, 4, 5], max_rows=5000)
    assert r.metadata["n_rows"] == 5000 and r.metadata["data_kind"] == "catalog"
    t = Table.read(r.files[0])
    assert 180.3 < np.median(t["ra"]) < 181.7 and -0.2 < np.median(t["dec"]) < 1.2
    assert 0.2 < np.std(t["e1"]) < 0.4 and set(np.unique(t["tomo_bin"])) <= {3, 4, 5}
    # the dispatch subset spec rebuilds the identical catalogue (what a facility node would do)
    cat = kernels.load_catalog(r.files[0])
    kernels.kids1000_subset(out_file=str(tmp_path / "rebuilt.fits"), **_subset_spec(cat))
    t2 = Table.read(tmp_path / "rebuilt.fits")
    assert len(t2) == len(t) and np.array_equal(t2["e1"], t["e1"]) and np.array_equal(t2["ra"], t["ra"])
