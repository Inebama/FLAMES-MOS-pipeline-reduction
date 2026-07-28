"""
This code was commented and edited by AI
convective_turnover.py
──────────────────────────────────────────────────────────────────────────────
Convective turnover time for RGB/AGB red giant stars.

PHYSICS SUMMARY
───────────────
The convective turnover time is:

    τ_conv = H_p / v_c

where H_p is the local pressure scale height and v_c is the convective
velocity from mixing-length theory (MLT).

Two estimates are computed:

  (1) LOCAL (photospheric): evaluated at optical depth τ_Ross ≈ 2/3.
      Valid for estimating pulsation driving/damping efficiency near the
      stellar surface. Requires ρ_phot, obtained self-consistently from:

        - Eddington grey atmosphere:  P_phot = (2g) / (3 κ)
        - Ideal gas:                  P_phot = ρ k_B T / (μ m_H)
        - H⁻ opacity power law:       κ ∝ Z_rel ρ^{1/2} T^9
              (Unsöld 1955; Mihalas 1978; Bohm-Vitense 1992)

      Combining these three, with A calibrated to the solar photosphere
      (ρ_phot,☉ ≈ 1.25×10⁻⁷ g cm⁻³ at T_☉=5778K, g_☉=27420 cm s⁻²):

        ρ_phot = [ 2g μ m_H / (3 A Z_rel k_B T^10) ]^{2/3}

      The MLT velocity in the EFFICIENT convection limit
      (Böhm-Vitense 1958; Kippenhahn, Weigert & Weiss 2012, §7.3):

        v_c³ = (α/2)² × 2g H_p F_c / (ρ c_p T)

      where F_c = L/(4πR²) is the convective flux.

      NOTE: For red giants the efficient limit is appropriate throughout
      most of the convective envelope. At the solar photosphere it is
      marginal, but our A-coefficient calibration absorbs this.

  (2) GLOBAL (depth-integrated, interior): a rough estimate based on
      the integral:

        τ_global ≈ ∫ dr / v_c(r)

      approximated as:

        τ_global ≈ (R³ / (G M))^{1/2} × (L / L_☉)^{-1/6}
                    × f_cal

      where f_cal is calibrated using stellar evolution models.
      More precisely, this uses the result of Gunn et al. (1998):

        τ_global = τ_☉,global × (L/L_☉)^{-1/3} × (M/M_☉)^{1/3}
                              × (R/R_☉)^2 / (α/α_☉)

      Calibrated to τ_☉,global ≈ 12 days (Corsaro et al. 2021),
      α = 1.6 (solar calibrated), α_☉ = 1.6.

      This gives τ_global that is typically 100–1000× τ_local for
      red giants, consistent with detailed model grids.

REFERENCES
──────────────────────────────────────────────────────────────────────────────
  [1]  Böhm-Vitense, E. (1958), ZAp, 46, 108
       The mixing length theory — foundational.

  [2]  Unsöld, A. (1955), "Physik der Sternatmosphären"
       H⁻ opacity power-law approximation.

  [3]  Mihalas, D. (1978), "Stellar Atmospheres", 2nd ed., W.H. Freeman
       H⁻ opacity and grey atmosphere models.

  [4]  Noyes, R. W., et al. (1984), ApJ, 279, 763
       Observational calibration of convective turnover time and Rossby
       number from chromospheric activity of solar-type stars.

  [5]  Gunn, A. G., Mitrou, C. K. & Doyle, J. G. (1998), MNRAS, 296, 150
       Global convective turnover time scaling relation for active stars.

  [6]  Kippenhahn, R., Weigert, A. & Weiss, A. (2012),
       "Stellar Structure and Evolution", 2nd ed., Springer  §§ 7.1–7.4
       Complete MLT derivation; efficient-limit velocity formula.

  [7]  Corsaro, E., et al. (2021), Nature Astronomy, 5, 272
       τ_global for red giants from asteroseismic data;
       solar reference value τ_☉,global ≈ 12 days.

  [8]  Cranmer, S. R. & Saar, S. H. (2011), ApJ, 741, 54
       MLT-based convective velocity at the photosphere for activity
       modelling; Eddington + H⁻ opacity approach.
──────────────────────────────────────────────────────────────────────────────
"""

import numpy as np


# ── Physical constants (CGS — everything in this module uses CGS) ──────────
_k_B   = 1.380649e-16   # erg K⁻¹
_m_H   = 1.67262192e-24 # g
_G_cgs = 6.67430e-8     # cm³ g⁻¹ s⁻²
_sigma = 5.670374e-5    # erg cm⁻² s⁻¹ K⁻⁴  (Stefan-Boltzmann)

# ── Solar constants (CGS) ──────────────────────────────────────────────────
_T_SUN      = 5778.0        # K
_LOG_G_SUN  = 4.438         # log10(g_☉ / [cm s⁻²]);  g_☉ = 27420 cm s⁻²
_R_SUN      = 6.957e10      # cm
_M_SUN      = 1.989e33      # g
_L_SUN      = 3.828e33      # erg s⁻¹

# Solar photospheric density (from standard solar models e.g. ATLAS9, MARCS):
_RHO_PHOT_SUN = 1.25e-7     # g cm⁻³   at τ_Ross = 2/3

# ── Global turnover time reference ────────────────────────────────────────
_TAU_GLOBAL_SUN = 12.0      # days  (Corsaro et al. 2021)

# ── H⁻ opacity calibration coefficient (CGS: κ in cm² g⁻¹) ───────────────
# κ_H⁻ = A_HM * Z_rel * ρ^{1/2} * T^9
# A_HM calibrated so that ρ_phot(T_☉, g_☉, μ_☉=0.6, Z_rel=1) = 1.25e-7 g/cm³
# Derivation: from ρ^{3/2} = 2 g μ m_H / (3 A_HM Z_rel k_B T^10)
# → A_HM = 2 g_☉ μ_☉ m_H / (3 Z_☉ k_B T_☉^10 × ρ_phot,☉^{3/2})
def _compute_A_HM():
    g_sun   = 10.0 ** _LOG_G_SUN              # cm s⁻²
    mu_sun  = 0.6                              # partially ionised photosphere
    rho_32  = _RHO_PHOT_SUN ** 1.5
    T10     = _T_SUN ** 10
    return 2.0 * g_sun * mu_sun * _m_H / (3.0 * 1.0 * _k_B * T10 * rho_32)

_A_HM = _compute_A_HM()


def photospheric_density(g_cgs, T_eff, mu, Z_rel):
    """
    Photospheric gas density at τ_Ross = 2/3 from H⁻ opacity + Eddington
    grey atmosphere + ideal gas (closed-form, no iteration needed).

    Parameters
    ----------
    g_cgs  : float  — surface gravity   [cm s⁻²]
    T_eff  : float  — effective temp    [K]
    mu     : float  — mean mol. weight  (dimensionless)
    Z_rel  : float  — Z/Z_☉ metallicity ratio

    Returns
    -------
    rho  [g cm⁻³]
    """
    # ρ^{3/2} = 2 g μ m_H / (3 A_HM Z_rel k_B T^10)
    rho_32 = 2.0 * g_cgs * mu * _m_H / (3.0 * _A_HM * Z_rel * _k_B * T_eff**10)
    return rho_32 ** (2.0 / 3.0)


def convective_turnover_local(log_L, T_eff, log_g,
                               alpha_MLT = 1.6,
                               gamma     = 5.0/3.0,
                               mu        = 1.3,
                               Z_rel     = 0.2):
    """
    LOCAL (photospheric) convective turnover time from MLT.

    Uses efficient-convection MLT velocity (Böhm-Vitense 1958 /
    Kippenhahn et al. 2012, §7.3):

        v_c³ = (α/2)² × 2g H_p F_c / (ρ c_p T)

    Parameters
    ----------
    log_L    : float — log10(L/L_☉)
    T_eff    : float — effective temperature [K]
    log_g    : float — log10(g / [cm s⁻²]) standard spectroscopic
    alpha_MLT: float — mixing length parameter (default 1.6, solar calibrated)
    gamma    : float — adiabatic index (default 5/3)
    mu       : float — mean molecular weight of photospheric gas
                       (default 1.3 for a cool, largely neutral RGB atmosphere)
    Z_rel    : float — Z/Z_☉ (default 0.2 for 47 Tuc, [Fe/H]≈−0.7)

    Returns
    -------
    dict with:
        H_p_cm        : pressure scale height at photosphere  [cm]
        rho_phot      : photospheric density                  [g cm⁻³]
        kappa_phot    : derived photospheric opacity          [cm² g⁻¹]
        v_conv_cms    : convective velocity (MLT)             [cm s⁻¹]
        v_conv_kms    : convective velocity                   [km s⁻¹]
        tau_local_s   : local turnover time                   [s]
        tau_local_d   : local turnover time                   [days]
        R_Rsun        : radius                                [R_☉]
        M_Msun        : mass (from g and R)                   [M_☉]
    """
    g_cgs   = 10.0 ** log_g                             # cm s⁻²
    g_sun   = 10.0 ** _LOG_G_SUN                        # cm s⁻²
    L_ratio = 10.0 ** log_L                             # L/L_☉

    # ── Radius and mass ───────────────────────────────────────────────────
    R_ratio  = np.sqrt(L_ratio) * (_T_SUN / T_eff) ** 2  # R/R_☉
    M_ratio  = (g_cgs / g_sun) * R_ratio ** 2             # M/M_☉
    R_cm     = R_ratio * _R_SUN                            # cm

    # ── Pressure scale height ─────────────────────────────────────────────
    H_p = _k_B * T_eff / (mu * _m_H * g_cgs)             # cm

    # ── Convective (radiative) flux ───────────────────────────────────────
    L_cgs = L_ratio * _L_SUN                               # erg s⁻¹
    F_c   = L_cgs / (4.0 * np.pi * R_cm ** 2)             # erg cm⁻² s⁻¹

    # ── Photospheric density from H⁻ opacity + Eddington ─────────────────
    rho_phot = photospheric_density(g_cgs, T_eff, mu, Z_rel)  # g cm⁻³

    # ── Derived opacity (diagnostic) ─────────────────────────────────────
    # From Eddington + ideal gas: κ = 2μ m_H g / (3 ρ k_B T)
    kappa = 2.0 * mu * _m_H * g_cgs / (3.0 * rho_phot * _k_B * T_eff)

    # ── Specific heat at constant pressure (ideal gas) ───────────────────
    c_p = (gamma / (gamma - 1.0)) * _k_B / (mu * _m_H)   # erg g⁻¹ K⁻¹

    # ── Convective velocity (efficient MLT limit) ─────────────────────────
    # v_c³ = (α/2)² × 2 g H_p F_c / (ρ c_p T)
    vc3 = (alpha_MLT / 2.0) ** 2 * 2.0 * g_cgs * H_p * F_c / (rho_phot * c_p * T_eff)
    v_c = vc3 ** (1.0 / 3.0)                              # cm s⁻¹

    # ── Turnover time ─────────────────────────────────────────────────────
    tau_s = H_p / v_c                                     # seconds
    tau_d = tau_s / 86400.0                               # days

    return {
        "H_p_cm":       H_p,
        "rho_phot":     rho_phot,
        "kappa_phot":   kappa,
        "v_conv_cms":   v_c,
        "v_conv_kms":   v_c / 1e5,
        "tau_local_s":  tau_s,
        "tau_local_d":  tau_d,
        "R_Rsun":       R_ratio,
        "M_Msun":       M_ratio,
    }


def convective_turnover_global(log_L, T_eff, log_g,
                                alpha_MLT     = 1.6,
                                tau_sun_global = _TAU_GLOBAL_SUN):
    """
    GLOBAL convective turnover time (approximate).

    Based on the scaling of Gunn et al. (1998), calibrated to the solar
    global turnover time from Corsaro et al. (2021):

        τ_global / τ_☉ = (L/L_☉)^{-1/3} × (M/M_☉)^{1/3} × (R/R_☉)^2

    This represents the integral ∫ dr/v_c(r) over the convection zone,
    approximated assuming that L_conv = L throughout (fully convective
    envelope) and using the mean convective velocity from the stellar
    luminosity and mass.

    NOTE: Stellar structure models are required for an accurate τ_global.
    This scaling can be off by factor 2–5× depending on the assumed
    M_core (core mass excluded from the convective envelope).

    Parameters
    ----------
    Same as convective_turnover_local (alpha_MLT and Z_rel not used here).
    tau_sun_global : float — solar global τ_conv in days (default 12 d,
                              Corsaro et al. 2021)

    Returns
    -------
    dict with:
        tau_global_d : global turnover time  [days]
        R_Rsun, M_Msun
    """
    g_cgs  = 10.0 ** log_g
    g_sun  = 10.0 ** _LOG_G_SUN
    L_ratio = 10.0 ** log_L
    R_ratio = np.sqrt(L_ratio) * (_T_SUN / T_eff) ** 2
    M_ratio = (g_cgs / g_sun) * R_ratio ** 2

    # Gunn et al. 1998 scaling
    tau_ratio = (L_ratio ** (-1.0/3.0) *
                 M_ratio ** ( 1.0/3.0) *
                 R_ratio **   2.0)
    tau_global = tau_sun_global * tau_ratio

    return {
        "tau_global_d": tau_global,
        "R_Rsun":       R_ratio,
        "M_Msun":       M_ratio,
    }


def convective_turnover_full(log_L, T_eff, log_g, **kwargs):
    """
    Convenience wrapper — returns both local and global estimates plus
    a ratio useful for diagnosing the pulsation regime.

    Also returns:
        P_max_days     : period at νmax   [days]  (for comparison)
        ratio_Pmax_tau : P_max / τ_local  (>1 suggests stochastic driving)
    """
    loc = convective_turnover_local(log_L, T_eff, log_g, **kwargs)
    glb = convective_turnover_global(log_L, T_eff, log_g)

    # ν_max period (recomputed here for convenience)
    g_cgs  = 10.0 ** log_g
    g_sun  = 10.0 ** _LOG_G_SUN
    nu_max = 3090.0 * (g_cgs / g_sun) * (_T_SUN / T_eff) ** 0.5  # μHz
    P_max  = 1.0 / (nu_max * 1e-6) / 86400.0                       # days

    return {
        **loc,
        "tau_global_d":   glb["tau_global_d"],
        "P_max_days":     P_max,
        "ratio_Pmax_tau": P_max / loc["tau_local_d"],
    }


def convective_turnover_table(log_L_arr, T_eff_arr, log_g_arr, **kwargs):
    """Vectorised wrapper — accepts arrays, returns dict of arrays."""
    log_L_arr = np.asarray(log_L_arr, dtype=float)
    T_eff_arr = np.asarray(T_eff_arr, dtype=float)
    log_g_arr = np.asarray(log_g_arr, dtype=float)
    results = [
        convective_turnover_full(lL, T, lg, **kwargs)
        for lL, T, lg in zip(log_L_arr, T_eff_arr, log_g_arr)
    ]
    keys = results[0].keys()
    return {k: np.array([r[k] for r in results]) for k in keys}


# ── Demo ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    print("=" * 72)
    print("SOLAR SELF-CHECK")
    print("=" * 72)
    # For the Sun we need log_g = log(g_sun) = 4.438, log_L=0, T=5778K, mu=0.6
    # (partially ionised photosphere), Z_rel=1
    sol = convective_turnover_local(0.0, 5778.0, 4.438, mu=0.6, Z_rel=1.0)
    print(f"  ρ_phot,☉  = {sol['rho_phot']:.3e} g/cm³  (expected ~1.25e-7)")
    print(f"  κ_phot,☉  = {sol['kappa_phot']:.3f} cm²/g  (expected ~0.18)")
    print(f"  v_conv,☉  = {sol['v_conv_kms']:.3f} km/s  (expected ~0.5-2 km/s)")
    print(f"  τ_local,☉ = {sol['tau_local_s']/60:.1f} min  (expected ~5-15 min for granulation)")

    sol_g = convective_turnover_global(0.0, 5778.0, 4.438)
    print(f"  τ_global,☉= {sol_g['tau_global_d']:.1f} days (by definition = 12 d)")

    print()
    print("=" * 72)
    print("47 TUC RGB SAMPLE  (Z_rel=0.2 for [Fe/H]≈-0.7, mu=1.3)")
    print("=" * 72)

    sample = [
        ("MED_RGB_3872",  2.12, 4483, 1.81),
        ("MED_RGB_557",   2.14, 4500, 1.80),
        ("MED_RGB_1807",  2.70, 4082, 1.05),
        ("MED_RGB_85",    3.05, 3864, 0.59),
        ("MED_VAR_618",   3.27, 3545, 0.24),
        ("MED_VAR_1399",  3.27, 3532, 0.23),
        ("MED_VAR_1373",  3.29, 3493, 0.19),
    ]

    hdr  = f"{'Star':<20} {'R/R☉':>6} {'M/M☉':>6} "
    hdr += f"{'v_c [km/s]':>11} {'H_p [km]':>10} "
    hdr += f"{'τ_local [d]':>12} {'τ_global [d]':>13} "
    hdr += f"{'P_max [d]':>10} {'P/τ_loc':>8}"
    print(hdr)
    print("-" * len(hdr))

    for name, lL, Teff, lg in sample:
        r = convective_turnover_full(lL, Teff, lg)
        print(
            f"{name:<20} "
            f"{r['R_Rsun']:6.1f} "
            f"{r['M_Msun']:6.3f} "
            f"{r['v_conv_kms']:11.2f} "
            f"{r['H_p_cm']/1e5:10.0f} "
            f"{r['tau_local_d']:12.3f} "
            f"{r['tau_global_d']:13.1f} "
            f"{r['P_max_days']:10.2f} "
            f"{r['ratio_Pmax_tau']:8.1f}"
        )

    print()
    print("P/τ_local >> 1  →  stochastic (solar-like) oscillations")
    print("P/τ_local  ~ 1  →  transition regime (semi-regular pulsators)")
    print("P/τ_local  < 1  →  convective damping dominates")
