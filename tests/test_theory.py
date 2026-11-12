"""predict_shear_correlation tests. Offline: synthetic n(z) files; network: compare to served data vectors."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("camb")

from tools.theory import kernels, predict_shear_correlation  # noqa: E402

TH = [2.0, 5.0, 10.0, 30.0, 100.0]


def _nz_file(path):
    z = np.linspace(0.0, 3.0, 301)
    n1 = np.exp(-0.5 * ((z - 0.5) / 0.15) ** 2)
    n2 = np.exp(-0.5 * ((z - 0.9) / 0.2) ** 2)
    np.savetxt(path, np.c_[z, n1, n2], header="z n1 n2")
    return str(path)


@pytest.fixture(scope="module")
def base(tmp_path_factory):
    d = tmp_path_factory.mktemp("th")
    nz = _nz_file(d / "nz.txt")
    r = predict_shear_correlation(str(d), nz_source=nz, theta_arcmin=TH, S8=0.8, Om=0.3, plot=True)
    return d, nz, r


def test_prediction_shape_and_physics(base):
    d, nz, r = base
    assert r.metadata["data_kind"] == "summary_statistic" and r.metadata["bin_pairs"] == [[1, 1], [1, 2], [2, 2]]
    df = pd.read_csv(r.files[0], comment="#")
    assert len(df) == 3 * len(TH) and r.files[1].endswith(".png")
    for _, g in df.groupby(["bin1", "bin2"]):
        assert np.all(g.xip > 0) and np.all(np.diff(g.xip.values) < 0)   # positive, falling
        assert g.xim.iloc[0] < g.xip.iloc[0]
    x11 = df[(df.bin1 == 1) & (df.bin2 == 1)].xip.values
    x22 = df[(df.bin1 == 2) & (df.bin2 == 2)].xip.values
    assert np.all(x22 > 1.5 * x11)                     # deeper sources lens more
    # order of magnitude: xi+ ~ 1e-5 - 1e-4 at 10 arcmin for z_s ~ 0.9
    assert 5e-6 < x22[2] < 3e-4


def test_amplitude_scales_with_sigma8(base):
    d, nz, r = base
    r2 = predict_shear_correlation(str(d), nz_source=nz, theta_arcmin=TH, S8=0.8 * 1.1, Om=0.3, plot=False)
    a = pd.read_csv(r.files[0], comment="#").xip.values
    b = pd.read_csv(r2.files[0], comment="#").xip.values
    ratio = b / a
    assert np.all((ratio > 1.19) & (ratio < 1.1 ** 3.5))   # ~s8^2 (linear); steeper on small scales
    assert ratio[0] > ratio[4]


def test_baryon_constant_ratio_and_hydroemu_format(base):
    d, nz, _ = base
    p = d / "hydro.csv"
    k = np.geomspace(1e-4, 1e3, 50)
    with open(p, "w") as fh:  # hydroemu-mcp-server layout
        fh.write("# label: Pk_ratio z=0\n# observable: Pk_ratio\n# redshift: 0.0\nx,y_mean,y_std\n")
        for kk in k:
            fh.write(f"{kk},0.9,0.01\n")
    r = predict_shear_correlation(str(d), nz_source=nz, theta_arcmin=TH, S8=0.8, baryon_file=str(p), plot=False)
    df = pd.read_csv(r.files[0], comment="#")
    assert np.allclose(df.ratio_xip, 0.9, atol=2e-3) and np.allclose(df.ratio_xim, 0.9, atol=2e-3)
    assert r.metadata["baryons"]["z_dependent"] is False


def test_baryon_z_dependent_long_format(tmp_path):
    k = np.geomspace(0.01, 50, 40)
    rows = [(kk, z, 1.0 if z < 0.5 else 0.8) for z in (0.0, 0.5, 3.0) for kk in k]
    pd.DataFrame(rows, columns=["k", "z", "ratio"]).to_csv(tmp_path / "b.csv", index=False)
    b = kernels.read_baryon_ratio(str(tmp_path / "b.csv"))
    assert b["ratio"].shape == (3, 40)
    S = kernels.baryon_factor(b, np.array([1.0, 1.0, 1.0]), np.array([0.0, 0.25, 2.0]), h=0.7)
    assert np.allclose(S, [1.0, 0.9, 0.8])
    # below the tabulated k the ratio is 1
    assert kernels.baryon_factor(b, np.array([1e-4]), np.array([2.0]), h=0.7)[0] == 1.0


def test_bad_inputs(tmp_path):
    with pytest.raises(ValueError, match="nz_source"):
        predict_shear_correlation(str(tmp_path), nz_source="nope")


@pytest.mark.network
@pytest.mark.parametrize("dataset,cosmo", [
    ("kids_1000", dict(Om=0.27, S8=0.76, h=0.7, Ob=0.046, ns=0.97)),
    ("des_y3", dict(Om=0.29, S8=0.78, h=0.68, Ob=0.048, ns=0.97)),
])
def test_validation_against_served_xip(tmp_path, dataset, cosmo):
    """Survey-like cosmology: data / GG-only theory within ~20% (IA lowers the data slightly)."""
    r = predict_shear_correlation(str(tmp_path), nz_source=dataset, compare_dataset=dataset, plot=False, **cosmo)
    c = r.metadata["comparison"]
    assert 0.8 < c["A"] < 1.2, c
