"""Weak-lensing tool tests. Offline tests use synthetic files; network tests hit the real portals."""

import json

import numpy as np
import pytest
from astropy.io import fits

from tools.chains import compare_s8_with_published, kernels
from tools.datavectors import io
from tools.registry import DATASETS, PUBLISHED_S8


def _write_2pt(path):
    n = 6
    cov = np.diag(np.full(2 * n, 1e-10)) + 1e-12
    h0 = fits.PrimaryHDU()
    hc = fits.ImageHDU(cov, name="COVMAT")
    hc.header.update({"NAME_0": "xip", "STRT_0": 0, "NAME_1": "xim", "STRT_1": n})
    cols = lambda v: fits.BinTableHDU.from_columns([  # noqa: E731
        fits.Column("BIN1", "K", array=[1, 1, 1, 2, 2, 2]), fits.Column("BIN2", "K", array=[1, 1, 1, 2, 2, 2]),
        fits.Column("ANGBIN", "K", array=[0, 1, 2, 0, 1, 2]), fits.Column("VALUE", "D", array=v),
        fits.Column("ANG", "D", array=[3.0, 10, 30, 3, 10, 30])])
    xp, xm = cols(np.full(n, 1e-4)), cols(np.full(n, 5e-5))
    xp.name, xm.name = "xip", "xim"
    z = np.linspace(0.005, 2.995, 300)
    nz = fits.BinTableHDU.from_columns([fits.Column("Z_LOW", "D", array=z - 0.005),
                                        fits.Column("Z_MID", "D", array=z), fits.Column("Z_HIGH", "D", array=z + 0.005),
                                        fits.Column("BIN1", "D", array=np.exp(-((z - 0.5) / 0.1) ** 2)),
                                        fits.Column("BIN2", "D", array=np.exp(-((z - 1.0) / 0.1) ** 2))])
    nz.name = "nz_source"
    fits.HDUList([h0, hc, xp, xm, nz]).writeto(path)


def test_read_2pt_fits_blocks_and_nz(tmp_path):
    p = tmp_path / "t.fits"
    _write_2pt(p)
    dv = io.read_2pt_fits(p)
    assert [o[0] for o in dv["order"]] == ["xip", "xim"]
    assert io.block_slice(dv, "xim") == slice(6, 12)
    nz = dv["nz"]["nz_source"]
    zbar = np.trapezoid(nz["z"] * nz["bins"][1], nz["z"]) / np.trapezoid(nz["bins"][1], nz["z"])
    assert abs(zbar - 1.0) < 1e-3


def test_chain_kernel_recovers_gaussian(tmp_path):
    rng = np.random.default_rng(1)
    om = rng.normal(0.3, 0.03, 40000)
    S8 = rng.normal(0.78, 0.02, 40000)
    s8 = S8 / np.sqrt(om / 0.3)
    p = tmp_path / "c.txt"
    with open(p, "w") as fh:
        fh.write("#cosmological_parameters--omega_m\tcosmological_parameters--sigma_8\tweight\n")
        np.savetxt(fh, np.c_[om, s8, np.ones_like(om)])
    cols, data = kernels.parse_cosmosis_text(p)
    r = kernels.summarize_columns(kernels.standard_columns(cols, data))
    assert abs(r["params"]["S_8"]["mean"] - 0.78) < 0.002
    assert abs(r["params"]["S_8"]["std"] - 0.02) < 0.002
    assert abs(r["corr_S8_Om"]) < 0.05


def test_compare_s8_published_tensions(tmp_path):
    r = compare_s8_with_published(str(tmp_path))
    t = r.metadata["tension_sigma"]
    assert 2.4 < t["DES Y3 3x2pt"] < 2.8
    assert t["KiDS-Legacy COSEBIs"] < 1.0
    assert set(PUBLISHED_S8) <= set(t)
    assert PUBLISHED_S8["HSC Y3 real space"]["S8"] == 0.769


def test_registry_and_y6_datavector_placeholder(tmp_path):
    from tools import registry
    from tools.datavectors import fetch_wl_datavector
    y6 = DATASETS["des_y6"]
    assert y6["datavector_status"]["available"] is False and registry.DES_Y6_DATAVECTORS_PUBLIC is False
    assert "3x2pt_lcdm" in y6["chains"] and y6["products"]["published_result"]
    with pytest.raises(io.DataNotPublicError, match="not yet public"):
        io.load_statistic("des_y6", "xip")
    with pytest.raises(io.DataNotPublicError, match="releases/y6a2"):
        io.load_nz("des_y6")
    r = fetch_wl_datavector("des_y6", "xip", str(tmp_path))
    assert r.metadata["available"] is False and r.files == [] and "not yet public" in r.message
    assert "DES_Y6_DATAVECTORS_PUBLIC" in r.metadata["enable_hint"]
    for k in ("des_y3", "kids_1000", "kids_legacy", "hsc_y3"):
        assert set(io.STATISTICS[k]) and DATASETS[k]["citation"]


def test_hsc_chain_is_credentialed(tmp_path, monkeypatch):
    monkeypatch.setenv("COSMO_MCP_DATA", str(tmp_path))
    from tools.chains import chain_source
    with pytest.raises(RuntimeError, match="credentials"):
        chain_source("hsc_y3", "real_lcdm")


def test_analyze_dispatch_branch(tmp_path, monkeypatch):
    import mcp_server.dispatch as d
    from tools.chains import analyze_wl_chain
    seen = {}

    def fake_run_kernel(function, args, pip_deps=None, duration=600, **kw):
        seen.update(function=function, args=args)
        fake = {"n_samples": 1, "n_eff": 1.0, "corr_S8_Om": 0.0,
                "params": {p: {"mean": 0.8, "std": 0.02} for p in ("S_8", "omega_m", "sigma_8")}}
        return {"result": fake, "host": "polaris"}

    monkeypatch.setattr(d, "remote_site", lambda: "polaris")
    monkeypatch.setattr(d, "run_kernel", fake_run_kernel)
    r = analyze_wl_chain("kids_1000", "cosebis", str(tmp_path))
    assert r.metadata["computed_on"] == "polaris"
    assert seen["function"] == "chains.kernels.chain_marginals"
    assert seen["args"]["url"].endswith(".tgz") and seen["args"]["member"].endswith("output_multinest_C.txt")


@pytest.mark.network
def test_des_y3_datavector_and_snr(tmp_path):
    from tools.datavectors import compute_signal_to_noise, fetch_wl_datavector
    r = fetch_wl_datavector("des_y3", "xip", str(tmp_path))
    assert r.metadata["n_points"] == 200 and r.metadata["n_bin_pairs"] == 10
    s = compute_signal_to_noise(["des_y3"], ["xip", "xim"], str(tmp_path)).metadata["results"][0]
    assert 40 < s["snr"] < 50


@pytest.mark.network
def test_des_y3_3x2pt_s8(tmp_path):
    from tools.chains import analyze_wl_chain
    r = analyze_wl_chain("des_y3", "3x2pt_lcdm", str(tmp_path))
    p = r.metadata["params"]["S_8"]
    assert abs(p["mean"] - 0.776) < 0.005 and abs(p["std"] - 0.017) < 0.003
    json.loads(open(r.files[0]).read())


@pytest.mark.network
def test_kids_1000_cosebis_s8(tmp_path):
    from tools.chains import analyze_wl_chain
    r = analyze_wl_chain("kids_1000", "cosebis", str(tmp_path))
    m = r.metadata
    # paper quotes MAP 0.759 +0.024/-0.021; mean is lower (skewed posterior)
    assert abs(m["max_post_sample"]["S_8"] - 0.759) < 0.005
    assert abs(m["params"]["S_8"]["mode"] - 0.759) < 0.01


@pytest.mark.network
def test_hsc_y3_real_and_fourier(tmp_path):
    from tools.datavectors import fetch_wl_datavector
    assert fetch_wl_datavector("hsc_y3", "xip", str(tmp_path)).metadata["n_bin_pairs"] == 10
    nz = io.load_nz("hsc_y3")
    assert len(nz["bins"]) == 4


def _write_nautilus(path, n=30000, seed=2, extra_cols=()):
    rng = np.random.default_rng(seed)
    om = rng.normal(0.32, 0.02, n)
    S8 = rng.normal(0.79, 0.011, n)
    s8 = S8 / np.sqrt(om / 0.3)
    lw = np.zeros(n)
    lw[:50] = -np.inf  # non-finite log weights must be dropped
    names = ["cosmological_parameters--omega_m", "cosmological_parameters--h0", *extra_cols,
             "COSMOLOGICAL_PARAMETERS--SIGMA_8", "COSMOLOGICAL_PARAMETERS--S_8", "log_weight", "prior", "post"]
    cols = [om, np.full(n, 0.68)] + [rng.normal(-1, 0.1, n) for _ in extra_cols] + [s8, S8, lw, np.zeros(n),
                                                                                     np.zeros(n)]
    with open(path, "w") as fh:
        fh.write("#" + "\t".join(names) + "\n#sampler=nautilus\n")
        np.savetxt(fh, np.column_stack(cols))
        fh.write("#log_z=1.0\n")


def test_fetch_des_y6_chain_offline(tmp_path, monkeypatch):
    """Nautilus layout (log_weight, -inf rows) and csislike layout ('weights') both parse from the header."""
    monkeypatch.setenv("COSMO_MCP_DATA", str(tmp_path / "cache"))
    from tools.chains import fetch_des_y6_chain
    from tools.common import data_cache_dir
    from tools.registry import DES_Y6_CHAIN_FILES
    d = data_cache_dir() / "des_y6"
    d.mkdir(parents=True)
    _write_nautilus(d / DES_Y6_CHAIN_FILES["3x2pt_wcdm"]["file"], extra_cols=("cosmological_parameters--w",))
    r = fetch_des_y6_chain(str(tmp_path), "3x2pt_wcdm")
    P = r.metadata["params"]
    assert r.metadata["data_kind"] == "published_result" and "comparison_only" in r.metadata
    assert r.metadata["n_samples"] == 30000 - 50
    assert abs(P["S_8"]["mean"] - 0.79) < 0.002 and abs(P["omega_m"]["std"] - 0.02) < 0.002
    assert abs(P["w"]["mean"] + 1) < 0.01
    assert any(f.endswith("_S8_om.png") for f in r.files) and any(f.endswith(".json") for f in r.files)
    rng = np.random.default_rng(3)
    om, s8 = rng.normal(0.31, 0.01, 5000), rng.normal(0.8, 0.01, 5000)
    with open(d / DES_Y6_CHAIN_FILES["3x2pt_cmb_lcdm"]["file"], "w") as fh:
        fh.write("#cosmological_parameters--omega_m\tCOSMOLOGICAL_PARAMETERS--SIGMA_8\tweights\tpost\n")
        np.savetxt(fh, np.c_[om, s8, rng.integers(1, 6, 5000), np.zeros(5000)])
    r2 = fetch_des_y6_chain(str(tmp_path), "3x2pt_cmb_lcdm", plot=False)
    exp = np.average(s8 * np.sqrt(om / 0.3), weights=np.loadtxt(d / DES_Y6_CHAIN_FILES["3x2pt_cmb_lcdm"]["file"])[:, 2])
    assert abs(r2.metadata["params"]["S_8"]["mean"] - exp) < 1e-3
    with pytest.raises(RuntimeError, match="allow_large"):
        fetch_des_y6_chain(str(tmp_path), "3x2pt_nonlin_wcdm")


def test_published_s8_has_des_y6():
    assert PUBLISHED_S8["DES Y6 3x2pt"]["S8"] == 0.789 and PUBLISHED_S8["DES Y6 3x2pt"]["plus"] == 0.012


@pytest.mark.network
def test_des_y6_fiducial_3x2pt_s8(tmp_path):
    from tools.chains import fetch_des_y6_chain
    r = fetch_des_y6_chain(str(tmp_path), "3x2pt_lcdm", plot=False)
    P = r.metadata["params"]
    assert abs(P["S_8"]["mean"] - 0.789) < 0.003 and abs(P["S_8"]["std"] - 0.012) < 0.002
    assert abs(P["omega_m"]["mean"] - 0.325) < 0.01 and abs(P["sigma_8"]["mean"] - 0.759) < 0.01
