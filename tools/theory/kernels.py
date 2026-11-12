"""Cosmic-shear theory: CAMB halofit P(k) -> Limber C_ell -> Hankel xi+/xi- (numpy + scipy + camb).

No intrinsic alignments, shear-calibration or photo-z nuisance terms: this is the pure gravitational-
lensing (GG) prediction for a given n(z) and cosmology, optionally multiplied by a baryonic
suppression ratio P_hydro / P_dmo (k[, z]).

Equations (flat LCDM):
  q_i(chi) = 3/2 Om (H0/c)^2 chi / a(chi) int_chi^inf dchi' n_i(chi') (chi' - chi) / chi'
  C_ell^ij = int dchi q_i q_j / chi^2  P_NL(k = (ell + 1/2) / chi, z(chi)) * S(k, z)
  xi_+/-(theta) = 1 / (2 pi) int d ell  ell C_ell J_0/4(ell theta)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

C_KMS = 299792.458
ARCMIN = np.pi / 180.0 / 60.0


def require_camb():
    try:
        import camb
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise RuntimeError("predict_shear_correlation needs CAMB: pip install camb  (or pip install -e '.[theory]'). "
                           "pyccl is not used because its macOS wheels need libomp.") from exc
    return camb


def camb_power(Om: float, Ob: float, h: float, ns: float, sigma8: float, mnu: float = 0.06, w0: float = -1.0,
               zmax: float = 4.0, halofit: str = "takahashi"):
    """Nonlinear matter P(k, z) interpolator (k in 1/Mpc, P in Mpc^3) normalised to sigma8, plus CAMB results."""
    camb = require_camb()
    omnuh2 = mnu / 93.14
    omch2 = (Om - Ob) * h * h - omnuh2
    if omch2 <= 0:
        raise ValueError(f"Om={Om} too small for Ob={Ob} and mnu={mnu}")

    def _pars(As, nonlinear):
        p = camb.CAMBparams()
        p.set_cosmology(H0=100 * h, ombh2=Ob * h * h, omch2=omch2, mnu=mnu, omk=0.0)
        if w0 != -1.0:
            p.set_dark_energy(w=w0, dark_energy_model="ppf")
        p.InitPower.set_params(As=As, ns=ns)
        p.set_matter_power(redshifts=[0.0], kmax=20.0)
        p.WantCls = False
        if nonlinear:
            p.NonLinear = camb.model.NonLinear_pk
            p.NonLinearModel.set_params(halofit_version=halofit)
        return p

    As0 = 2.1e-9
    r0 = camb.get_results(_pars(As0, False))
    s80 = float(r0.get_sigma8_0())
    As = As0 * (sigma8 / s80) ** 2
    pars = _pars(As, True)
    res = camb.get_background(pars)
    PK = camb.get_matter_power_interpolator(pars, nonlinear=True, hubble_units=False, k_hunit=False, kmax=20.0,
                                            zmax=zmax, extrap_kmax=1e4, var1="delta_tot", var2="delta_tot")
    return PK, res, {"As": As, "omch2": omch2, "sigma8": sigma8, "halofit": halofit}


# ----------------------------------------------------------------------------- baryon suppression
def read_baryon_ratio(path: str) -> dict:
    """Baryonic P(k) ratio file -> {"k_hmpc": array, "z": array | None, "ratio": 2D (nz, nk) or 1D}.

    Accepted formats:
    - hydroemu-mcp-server Pk-ratio CSV: '# ...' comment lines, then header 'x,y_mean[,y_std]' (or 'x,y'),
      x = k [h/Mpc], y = P_hydro / P_dmo (one redshift; the '# redshift:' comment is recorded).
    - plain CSV with columns k (h/Mpc), z, ratio (long format, several z) or k, ratio.
    """
    import pandas as pd

    lines = Path(path).read_text().splitlines()
    zcomment = None
    for ln in lines:
        if ln.startswith("#") and "redshift" in ln.lower() and ":" in ln:
            try:
                zcomment = float(ln.split(":", 1)[1])
            except ValueError:
                pass
    df = pd.read_csv(path, comment="#")
    cols = {c.lower().strip(): c for c in df.columns}
    kcol = next((cols[c] for c in ("k", "k_hmpc", "k_h_mpc", "x") if c in cols), None)
    rcol = next((cols[c] for c in ("ratio", "y_mean", "y", "suppression") if c in cols), None)
    if kcol is None or rcol is None:
        raise ValueError(f"{path}: need a k column (k or x, h/Mpc) and a ratio column (ratio / y_mean / y); "
                         f"got {list(df.columns)}")
    if "z" in cols and df[cols["z"]].nunique() > 1:
        zs = np.sort(df[cols["z"]].unique())
        ks = np.sort(df[kcol].unique())
        grid = np.empty((len(zs), len(ks)))
        for a, z in enumerate(zs):
            sub = df[df[cols["z"]] == z].sort_values(kcol)
            grid[a] = np.interp(np.log(ks), np.log(sub[kcol].values), sub[rcol].values)
        return {"k_hmpc": ks, "z": zs, "ratio": grid, "source": str(path)}
    sub = df.sort_values(kcol)
    return {"k_hmpc": sub[kcol].values.astype(float), "z": None, "ratio": sub[rcol].values.astype(float),
            "z_of_file": zcomment, "source": str(path)}


def baryon_factor(b: dict, k_mpc: np.ndarray, z: np.ndarray, h: float) -> np.ndarray:
    """S(k, z) on matching arrays; ratio -> 1 below the tabulated k, clamped above; constant in z if 1D."""
    kh = np.log(np.clip(k_mpc / h, 1e-12, None))
    lk = np.log(b["k_hmpc"])

    def one(r):
        out = np.interp(kh, lk, r, left=1.0, right=r[-1])
        return out

    if b["z"] is None:
        return one(np.asarray(b["ratio"], float))
    zs = np.asarray(b["z"], float)
    vals = np.stack([one(r) for r in b["ratio"]])  # (nz_tab, ...)
    zc = np.clip(z, zs[0], zs[-1])
    i = np.clip(np.searchsorted(zs, zc) - 1, 0, len(zs) - 2)
    t = (zc - zs[i]) / (zs[i + 1] - zs[i])
    sel = np.indices(kh.shape)
    return (1 - t) * vals[(i, *sel)] + t * vals[(i + 1, *sel)]


# ----------------------------------------------------------------------------- Limber + Hankel
def lensing_kernels(z_nz: np.ndarray, nz_bins: list, chi_of_z, Om: float, h: float, z_grid: np.ndarray):
    """q_i(chi) on z_grid for each n(z) bin (normalised to unit integral)."""
    chi = chi_of_z(z_grid)
    a = 1.0 / (1.0 + z_grid)
    pref = 1.5 * Om * (100 * h / C_KMS) ** 2
    q = []
    for n in nz_bins:
        n = np.clip(np.interp(z_grid, z_nz, np.asarray(n, float), left=0.0, right=0.0), 0, None)
        norm = np.trapezoid(n, z_grid)
        if norm <= 0:
            raise ValueError("an n(z) bin integrates to zero on the redshift grid")
        n = n / norm  # per unit z
        # int_z^inf dz' n(z') (chi' - chi)/chi'  (n dz = n_chi dchi)
        frac = np.where(chi[None, :] > 0, (chi[None, :] - chi[:, None]) / np.where(chi > 0, chi, 1)[None, :], 0.0)
        frac = np.triu(np.clip(frac, 0, None))
        g = np.trapezoid(frac * n[None, :], z_grid, axis=1)
        q.append(pref * chi / a * g)
    return chi, np.array(q)


def limber_cls(ell: np.ndarray, chi: np.ndarray, z_grid: np.ndarray, q: np.ndarray, pairs: list, PK, h: float,
               baryons: dict | None = None) -> np.ndarray:
    """C_ell for each (i, j) pair (0-based) -> array (npairs, nell)."""
    m = chi > 0
    chi_, z_, q_ = chi[m], z_grid[m], q[:, m]
    k = (ell[:, None] + 0.5) / chi_[None, :]
    zz = np.broadcast_to(z_, k.shape)
    P = PK.P(zz.ravel(), k.ravel(), grid=False).reshape(k.shape)
    P = np.where(k < 1e4, P, 0.0)
    if baryons is not None:
        P = P * baryon_factor(baryons, k, zz, h)
    out = np.empty((len(pairs), len(ell)))
    for n, (i, j) in enumerate(pairs):
        integrand = (q_[i] * q_[j] / chi_ ** 2)[None, :] * P
        out[n] = np.trapezoid(integrand, chi_, axis=1)
    return out


def hankel_xi(ell_c: np.ndarray, cls: np.ndarray, theta_arcmin: np.ndarray, lmax: float = 1e5):
    """xi+ / xi- from C_ell by dense quadrature (log grid to ell=1000, then step 1 ... 4 to lmax)."""
    from scipy.special import jv

    L = np.unique(np.concatenate([np.geomspace(1.0, 1000.0, 3000), np.arange(1000.0, 10000.0, 1.0),
                                  np.arange(10000.0, lmax, 4.0)]))
    lc = np.log(ell_c)
    C = np.array([np.exp(np.interp(np.log(L), lc, np.log(np.clip(c, 1e-300, None)))) for c in cls])
    C[:, L > ell_c[-1]] = 0.0
    th = np.asarray(theta_arcmin, float) * ARCMIN
    x = L[:, None] * th[None, :]
    w = np.gradient(L)  # trapezoid-like weights on a non-uniform grid
    taper = np.where(L > 0.7 * lmax, np.cos(0.5 * np.pi * (L - 0.7 * lmax) / (0.3 * lmax)) ** 2, 1.0)
    kern = (w * L * taper)[:, None] / (2 * np.pi)
    xip = C @ (kern * jv(0, x))
    xim = C @ (kern * jv(4, x))
    return xip, xim


def predict_xi(z_nz, nz_bins, theta_arcmin, cosmo: dict, pairs: list | None = None, baryons: dict | None = None,
               halofit: str = "takahashi", nell: int = 160) -> dict:
    """Full chain for a set of bins; pairs 0-based. Returns ell, cls, xip, xim (each (npairs, ...)) and optional
    no-baryon versions when ``baryons`` is given."""
    nb = len(nz_bins)
    pairs = pairs or [(i, j) for i in range(nb) for j in range(i, nb)]
    zmax = float(min(max(np.asarray(z_nz)[np.any(np.asarray(nz_bins) > 0, axis=0)].max(), 1.0) + 0.1, 6.0))
    PK, res, info = camb_power(cosmo["Om"], cosmo["Ob"], cosmo["h"], cosmo["ns"], cosmo["sigma8"],
                               cosmo.get("mnu", 0.06), cosmo.get("w0", -1.0), zmax=zmax + 0.1, halofit=halofit)
    z_grid = np.concatenate([[0.0], np.geomspace(1e-3, zmax, 400)])
    chi, q = lensing_kernels(np.asarray(z_nz, float), nz_bins, res.comoving_radial_distance, cosmo["Om"],
                             cosmo["h"], z_grid)
    ell = np.geomspace(1.0, 1e5, nell)
    out = {"ell": ell, "pairs": pairs, "camb": info}
    cls = limber_cls(ell, chi, z_grid, q, pairs, PK, cosmo["h"])
    out["cls_dmo"] = cls
    out["xip_dmo"], out["xim_dmo"] = hankel_xi(ell, cls, theta_arcmin)
    if baryons is not None:
        clb = limber_cls(ell, chi, z_grid, q, pairs, PK, cosmo["h"], baryons)
        out["cls"] = clb
        out["xip"], out["xim"] = hankel_xi(ell, clb, theta_arcmin)
    else:
        out["cls"], out["xip"], out["xim"] = cls, out["xip_dmo"], out["xim_dmo"]
    return out
