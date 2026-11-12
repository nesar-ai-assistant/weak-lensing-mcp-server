"""Single source of truth for the datasets this server exposes.

Keys are short dataset ids used by every tool's ``dataset`` argument.
"access" is "bundled" (file in data/), "download" (fetched on first use
into the data cache), "credentialed" (needs a login; manual download into
the cache) or "unreleased".

Every URL below was downloaded successfully during development (2026-11-04)
unless the entry says otherwise.
"""

DES_Y3_DV = "https://desdr-server.ncsa.illinois.edu/despublic/y3a2_files/datavectors/"
DES_Y3_CHAINS = "https://desdr-server.ncsa.illinois.edu/despublic/y3a2_files/chains/"
KIDS1000_TGZ = "https://kids.strw.leidenuniv.nl/DR4/data_files/KiDS1000_cosmic_shear_data_release.tgz"
KIDS_LEGACY_TGZ = "https://kids.strw.leidenuniv.nl/sci_data/KiDS_Legacy_cosmic_shear_data_release.tar.gz"
CSL_RAW = "https://raw.githubusercontent.com/cosmosis-developers/cosmosis-standard-library/main/likelihood/"

DES_Y6_CHAINS = "https://desdr-server.ncsa.illinois.edu/despublic/y6a2_files/y6_key/chains_release/"
# chain key -> file in DES_Y6_CHAINS, size and content (sizes from HTTP HEAD, 2026-11-12).
DES_Y6_CHAIN_FILES: dict[str, dict] = {
    "3x2pt_lcdm": {"file": "fiducial32pt_wolensbin2_nonlin_lcdm_DV0728_R0820.txt", "size_mb": 370,
                   "desc": "FIDUCIAL 3x2pt LCDM (nonlinear galaxy bias), Nautilus"},
    "1x2pt_lcdm": {"file": "12pt_wolensbin2_lin_lcdm_DV0723_R0723.txt", "size_mb": 15,
                   "desc": "cosmic shear only (TATT IA), LCDM, Nautilus"},
    "1x2pt_nla_lcdm": {"file": "shear_nla_wolensbin2_lcdm_DV0723_R0723.txt", "size_mb": 10,
                       "desc": "cosmic shear only with NLA intrinsic alignments, LCDM, Nautilus"},
    "2x2pt_lcdm": {"file": "22pt_wolensbin2_lin_lcdm_DV0723_R0723.txt", "size_mb": 128,
                   "desc": "2x2pt (gammat + wtheta, linear bias), LCDM, Nautilus"},
    "3x2pt_linbias_lcdm": {"file": "32pt_wolensbin2_lin_lcdm_DV0723_R0723.txt", "size_mb": 144,
                           "desc": "3x2pt with linear galaxy bias, LCDM, Nautilus"},
    "3x2pt_fixmnu_lcdm": {"file": "fixmnu_32pt_wolensbin2_lin_lcdm_DV0728_R0909.txt", "size_mb": 76,
                          "desc": "3x2pt, fixed neutrino mass, LCDM, Nautilus"},
    "3x2pt_wcdm": {"file": "fiducial32pt_wolensbin2_lin_wcdm_DV0723_R0723.txt", "size_mb": 388,
                   "desc": "3x2pt wCDM (linear bias), Nautilus"},
    "3x2pt_nonlin_wcdm": {"file": "fiducial32pt_wolensbin2_nonlin_wcdm_DV0728_R0820.txt", "size_mb": 4486,
                          "desc": "3x2pt wCDM nonlinear bias - 4.5 GB, only on explicit request"},
    "cmb_lcdm": {"file": "fiducial_cmb_tp_external_R1006.txt", "size_mb": 12,
                 "desc": "external CMB (Planck/ACT/SPT T+P) alone, LCDM, Nautilus"},
    "3x2pt_cmb_lcdm": {"file": "chain_des_3x2pt_cmb_t_p_chain_csislike.txt", "size_mb": 2.5,
                       "desc": "DES Y6 3x2pt + CMB T+P, LCDM (emcee resample, 'csislike')"},
    "3x2pt_bao_sn_clusters_lcdm": {
        "file": "chain_des_3x2pt_des_bao_no_desi_des_sn_desi_dr2_bao_spt_clusters_chain_csislike.txt",
        "size_mb": 6.3, "desc": "3x2pt + DES BAO + DES SN + DESI DR2 BAO + SPT clusters, LCDM (csislike)"},
    "3x2pt_cmb_wcdm": {"file": "chain_des_3x2pt_wcdm_cmb_t_p_wcdm_chain_csislike.txt", "size_mb": 10.8,
                       "desc": "3x2pt + CMB T+P, wCDM (csislike)"},
    "3x2pt_bao_bbn_lcdm": {"file": "chain_des_y6_3x2pt_des_y6_bao_bbn_prior_chain_csislike.txt", "size_mb": 5.8,
                           "desc": "3x2pt + DES Y6 BAO + BBN prior, LCDM (csislike)"},
}

# DES Y6 2pt data vectors / covariance: NOT public yet (2026-11-12). One-line switch: flip to True and
# set DES_Y6_DV_URL to the released CosmoSIS 2pt FITS; fetch_wl_datavector('des_y6', ...) then works like des_y3.
DES_Y6_DATAVECTORS_PUBLIC = False
DES_Y6_DV_URL = ""
DES_Y6_DV_NOTICE = ("DES Y6 2pt data vectors/covariance not yet public; check https://des.ncsa.illinois.edu/releases/y6a2 "
                    "(the Y6 published chains ARE public: fetch_des_y6_chain).")

PLANCK18_S8 = {"mean": 0.832, "sigma": 0.013,
               "citation": "Planck Collaboration 2020, A&A 641, A6 (TT,TE,EE+lowE+lensing)"}

DATASETS: dict[str, dict] = {
    "des_y3": {
        "title": "DES Year 3 cosmic shear + 3x2pt (MagLim lenses; redMaGiC alternative)",
        "products": {"catalog": [],
                     "summary_statistic": ["xip", "xim", "gammat", "wtheta", "nz_source", "nz_lens", "covariance"],
                     "published_result": ["chains (~30: 1x2pt/2x2pt/3x2pt, LCDM/wCDM, +Planck)"]},
        "catalog_note": "The Y3 Metacalibration shape catalogue is public (des.ncsa.illinois.edu/releases/y3a2/"
                      "Y3key-catalogs) but only as DESY3_metacal_v03-004.h5 = 312 GB (plus 27 GB index cat, "
                      "37 GB SOMPZ) at desdr-server.ncsa.illinois.edu/despublic/y3a2_files/y3kp_cats/, HDF5 "
                      "with no subset service (verified listing 2026-11-04). Not fetched; use kids_1000.",
        "access": "download",
        "statistics": {"xip": "xi_+(theta) cosmic shear, 4x4 source bins (10 pairs x 20 angles)",
                       "xim": "xi_-(theta) cosmic shear",
                       "gammat": "galaxy-galaxy lensing gamma_t(theta), 6 lens x 4 source bins",
                       "wtheta": "galaxy clustering w(theta), 6 MagLim lens bins (auto)"},
        "files": {"datavector": DES_Y3_DV + "2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits",
                  "datavector_redmagic": DES_Y3_DV + "2pt_NG_final_2ptunblind_02_24_21_wnz_redmagic_covupdate.fits"},
        "chains": {"3x2pt_lcdm": "chain_3x2pt_lcdm_SR_maglim.txt",
                   "3x2pt_wcdm": "chain_3x2pt_wcdm_SR_maglim.txt",
                   "1x2pt_lcdm": "chain_1x2pt_lcdm_SR_maglim.txt",
                   "2x2pt_lcdm": "chain_2x2pt_lcdm_SR_maglim.txt",
                   "3x2pt_lcdm_redmagic": "chain_3x2pt_lcdm_SR_redmagic.txt",
                   "planck_lcdm": "chain_p-TTTEEE-lowE_lcdm.txt",
                   "planck_des_lcdm": "chain_planck_des_lcdm.txt"},
        "units": "theta in arcmin; xi dimensionless; n(z) per unit z on a 0.01-wide grid",
        "urls": [DES_Y3_DV, DES_Y3_CHAINS],
        "citation": "DES Collaboration, Abbott et al. 2022, PRD 105, 023520 (arXiv:2105.13549); "
                    "Amon et al. 2022 / Secco et al. 2022 (cosmic shear)",
        "reference_S8": {"3x2pt_lcdm": "0.776 +/- 0.017"},
        "tools": ["fetch_wl_datavector", "plot_shear_correlation", "plot_redshift_distributions",
                  "plot_covariance", "compute_signal_to_noise", "fetch_wl_chain", "analyze_wl_chain"],
    },
    "des_y6": {
        "title": "DES Year 6 3x2pt: published chains (2pt data vectors / covariance not yet public)",
        "products": {"catalog": [],
                     "summary_statistic": (["xip", "xim", "gammat", "wtheta", "nz_source", "nz_lens", "covariance"]
                                           if DES_Y6_DATAVECTORS_PUBLIC else []),
                     "published_result": ["chains (Nautilus: 1x2pt/2x2pt/3x2pt LCDM/wCDM; emcee 'csislike' "
                                          "combinations with CMB / BAO / SN / clusters)"]},
        "access": "download (chains, 2.5 MB - 390 MB each); data vectors not yet public",
        "chains_base": DES_Y6_CHAINS,
        "chains": {k: v["file"] for k, v in DES_Y6_CHAIN_FILES.items()},
        "chain_info": DES_Y6_CHAIN_FILES,
        "datavector_status": {
            "product": "DES Y6 2pt data vectors (xi+-, gammat, wtheta), n(z) and covariance",
            "available": DES_Y6_DATAVECTORS_PUBLIC,
            "notice": DES_Y6_DV_NOTICE,
            "enable": "set DES_Y6_DATAVECTORS_PUBLIC = True in tools/registry.py and put the released 2pt "
                      "FITS URL in DES_Y6_DV_URL (CosmoSIS 2pt FITS layout, as DES Y3).",
        },
        "files": {"datavector": DES_Y6_DV_URL},
        "units": "chains: CosmoSIS text; Nautilus files carry log_weight (natural log), csislike files "
                 "carry 'weights' (sample multiplicities)",
        "urls": [DES_Y6_CHAINS, "https://des.ncsa.illinois.edu/releases/y6a2"],
        "citation": "DES Collaboration 2026, 'Dark Energy Survey Year 6 Results: cosmological constraints from "
                    "galaxy clustering and weak lensing' (arXiv:2601.14559); 140M source shapes, 9M lenses, "
                    "~5000 deg^2",
        "reference_S8": {"3x2pt_lcdm": "0.789 +/- 0.012 (paper; Om = 0.333 +0.023/-0.028)",
                         "3x2pt_lcdm (this server, chain mean +/- std)": "0.7888 +/- 0.0112, "
                                                                          "Om = 0.325 +/- 0.023, "
                                                                          "sigma8 = 0.759 +/- 0.033"},
        "tools": ["fetch_des_y6_chain", "fetch_wl_chain", "analyze_wl_chain", "plot_s8_om_contours",
                  "compare_s8_with_published", "fetch_wl_datavector (returns the not-yet-public notice)"],
    },
    "kids_1000": {
        "title": "KiDS-1000 cosmic shear (Asgari et al. 2021): COSEBIs, band powers, xi+-",
        "products": {"catalog": ["gold WL shear catalogue (21,262,011 galaxies, 1006 tiles; subsets by tile / "
                               "patch / RA-Dec box / tomographic bin / row limit via HTTP Range reads)"],
                     "summary_statistic": ["xip", "xim", "cosebis", "bandpowers", "nz_source", "covariance"],
                     "published_result": ["MultiNest chains for each statistic"]},
        "access": "download (release tarball 17 MB; shear catalogue 17.7 GB FITS read as Range subsets)",
        "shear_catalog": {
            "url": "https://kids.strw.leidenuniv.nl/DR4/data_files/KiDS_DR4.1_ugriZYJHKs_SOM_gold_WL_cat.fits",
            "page": "https://kids.strw.leidenuniv.nl/DR4/KiDS-1000_shearcatalogue.php",
            "size_bytes": 17712469440,
            "n_rows": 21262011, "n_tiles": 1006,
            "tile_index": "data/kids1000_gold_tile_runs.npz",
            "zb_edges": [0.1, 0.3, 0.5, 0.7, 0.9, 1.2],
            "m_bias": [-0.009, -0.011, -0.015, 0.002, 0.007],
            "m_sigma": [0.019, 0.020, 0.017, 0.012, 0.010],
            "m_bias_ref": "Asgari et al. 2021, A&A 645, A104, Table 1",
            "citation": "Giblin et al. 2021, A&A 645, A105 (shear catalogue); Kuijken et al. 2019 (DR4)",
        },
        "statistics": {"xip": "xi_+ (5 bins, 15 pairs x 9 angles)", "xim": "xi_-",
                       "cosebis": "COSEBIs E_n, n=1..5 per pair (fiducial)",
                       "bandpowers": "C_EE band powers PeeE (8 ell bins)"},
        "archive": KIDS1000_TGZ,
        "files": {"xipm": "data_fits/xipm_KIDS1000_BlindC_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits",
                  "cosebis": "data_fits/cosebis_KIDS1000_BlindC_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits",
                  "bandpowers": "data_fits/bp_KIDS1000_BlindC_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits"},
        "chains": {"cosebis": "chains_and_config_files/main_chains_iterative_covariance/cosebis/chain/output_multinest_C.txt",
                   "bandpowers": "chains_and_config_files/main_chains_iterative_covariance/bp/chain/output_multinest_C.txt",
                   "xipm": "chains_and_config_files/main_chains_iterative_covariance/xipm/chain/output_multinest_C.txt"},
        "units": "theta in arcmin; COSEBIs dimensionless; band powers ell^2 C_ell/2pi convention of the release",
        "urls": [KIDS1000_TGZ, "https://kids.strw.leidenuniv.nl/sciencedata.php"],
        "citation": "Asgari et al. 2021, A&A 645, A104 (arXiv:2007.15633); Hildebrandt et al. 2021 (n(z))",
        "reference_S8": {"cosebis": "0.759 +0.024 -0.021 (MAP + PJ-HPD)"},
        "tools": ["fetch_wl_datavector", "plot_shear_correlation", "plot_redshift_distributions",
                  "plot_covariance", "compute_signal_to_noise", "fetch_wl_chain", "analyze_wl_chain",
                  "fetch_shear_catalog", "plot_shear_catalog", "measure_shear_correlation",
                  "measure_nz_from_catalog", "compare_measured_with_published_xi"],
    },
    "kids_legacy": {
        "title": "KiDS-Legacy cosmic shear (Wright et al. 2025): COSEBIs, band powers, xi+-",
        "products": {"catalog": [],
                     "summary_statistic": ["xip", "xim", "cosebis", "cosebis_b", "bandpowers", "bandpowers_b", "nz_source",
                               "covariance"],
                     "published_result": ["Nautilus chains (+ cosebis_plus_external)"]},
        "catalog_note": "The KiDS-Legacy (DR5) shear catalogue is not linked for public download on the KiDS "
                      "science-data pages as of 2026-11-04; use kids_1000 for catalogue-level work.",
        "access": "download",
        "statistics": {"xip": "xi_+ (6 bins, 21 pairs x 9 angles)", "xim": "xi_-",
                       "cosebis": "COSEBIs E_n, 20 modes per pair (analysis uses n<=6)",
                       "cosebis_b": "COSEBIs B_n (null test)",
                       "bandpowers": "PeeE band powers", "bandpowers_b": "PeeB (null test)"},
        "archive": KIDS_LEGACY_TGZ,
        "archive_size_mb": 131,
        "files": {"xipm": "data/KiDS_Legacy_xipm.fits", "cosebis": "data/KiDS_Legacy_cosebis.fits",
                  "bandpowers": "data/KiDS_Legacy_bandpowers.fits"},
        "chains": {"cosebis": "chains_and_config_files/cosebis/output_nautilus_cosebis_Fiducial.txt",
                   "bandpowers": "chains_and_config_files/bandpowers/output_nautilus_bandpowers_Fiducial.txt",
                   "xipm": "chains_and_config_files/xipm/output_nautilus_xipm_Fiducial.txt",
                   "cosebis_plus_external": "chains_and_config_files/cosebis_plus_external/output_nautilus_cosebis_plus_external.txt"},
        "units": "as KiDS-1000",
        "urls": [KIDS_LEGACY_TGZ, "https://kids.strw.leidenuniv.nl/DR5/legacy_wl.php"],
        "citation": "Wright et al. 2025, A&A 703, A158 (arXiv:2503.19441); Stolzner et al. 2025",
        "reference_S8": {"cosebis": "0.815 +0.016 -0.021 (Wright+2025)"},
        "tools": ["fetch_wl_datavector", "plot_shear_correlation", "plot_redshift_distributions",
                  "plot_covariance", "compute_signal_to_noise", "fetch_wl_chain", "analyze_wl_chain"],
    },
    "hsc_y3": {
        "title": "HSC Year 3 cosmic shear: real-space xi+- (Li+2023) and Fourier C_ell (Dalal+2023)",
        "products": {"catalog": [],
                     "summary_statistic": ["xip", "xim", "cl_ee", "nz_source", "covariance"],
                     "published_result": ["chains (credentialed)"]},
        "catalog_note": "The HSC Y3 shape catalogue needs an HSC-SSP account (hsc-release.mtk.nao.ac.jp returns "
                      "HTTP 401); not fetched by this server.",
        "access": "download (real-space FITS 0.45 MB, Fourier SACC 22 MB); chains credentialed",
        "statistics": {"xip": "xi_+ (4 bins, scale-cut, 130 points)", "xim": "xi_- (100 points)",
                       "cl_ee": "shear C_ell (4 bins, 10 pairs x 17 ell, 170 points)"},
        "files": {"real": CSL_RAW + "hsc_cosmic_shear/hsc_y3_real.fits",
                  "fourier": CSL_RAW + "hsc_cosmic_shear/hsc_y3_fourier_shear.sacc"},
        "note": "Data vectors come from the public CosmoSIS standard library mirror. The HSC portal "
                "(hsc-release.mtk.nao.ac.jp/archive/filetree/shape_catalog_y3/...) returns HTTP 401 without "
                "HSC-SSP credentials, so the HSC Y3 CHAINS are 'credentialed': download them manually "
                "(dalal23/hsc_y3_fourier_cosmic_shear.txt, li23/outputs/...) into the cache dir hsc_y3/ "
                "and fetch_wl_chain will pick them up. n(z) for both statistics is taken from the SACC "
                "(the real-space FITS carries no n(z) HDU; same source sample).",
        "chains": {"real_lcdm": "credentialed", "fourier_lcdm": "credentialed"},
        "urls": [CSL_RAW + "hsc_cosmic_shear/", "https://hsc-release.mtk.nao.ac.jp/doc/index.php/wly3/"],
        "citation": "Li et al. 2023, PRD 108, 123518 (arXiv:2304.00702); Dalal et al. 2023, PRD 108, 123519 "
                    "(arXiv:2304.00701)",
        "reference_S8": {"real": "0.769 +0.031 -0.034", "fourier": "0.776 +0.032 -0.033"},
        "tools": ["fetch_wl_datavector", "plot_shear_correlation", "plot_redshift_distributions",
                  "plot_covariance", "compute_signal_to_noise", "compare_s8_with_published"],
    },
}

# Published headline S8 values, for comparison only (compare_s8_with_published when a chain is not
# public or not yet fetched). Asymmetric errors as (plus, minus).
PUBLISHED_S8: dict[str, dict] = {
    "DES Y6 3x2pt": {"S8": 0.789, "plus": 0.012, "minus": 0.012,
                     "ref": "DES Collaboration 2026, arXiv:2601.14559 (Om = 0.333 +0.023/-0.028)"},
    "DES Y3 3x2pt": {"S8": 0.776, "plus": 0.017, "minus": 0.017,
                     "ref": "DES Collaboration 2022, PRD 105, 023520"},
    "DES Y3 cosmic shear": {"S8": 0.759, "plus": 0.025, "minus": 0.023,
                            "ref": "Amon+2022 / Secco+2022, PRD 105, 023514/023515"},
    "KiDS-1000 COSEBIs": {"S8": 0.759, "plus": 0.024, "minus": 0.021,
                          "ref": "Asgari+2021, A&A 645, A104"},
    "KiDS-Legacy COSEBIs": {"S8": 0.815, "plus": 0.016, "minus": 0.021,
                            "ref": "Wright+2025, A&A 703, A158"},
    "HSC Y3 real space": {"S8": 0.769, "plus": 0.031, "minus": 0.034,
                          "ref": "Li+2023, PRD 108, 123518"},
    "HSC Y3 Fourier": {"S8": 0.776, "plus": 0.032, "minus": 0.033,
                       "ref": "Dalal+2023, PRD 108, 123519"},
}
