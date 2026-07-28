"""
This code was commented and edited by AI
stellar_oscillation_params.py
──────────────────────────────────────────────────────────────────────────────
Functions to compute asteroseismic and structural parameters for red giant
stars from their fundamental photometric/spectroscopic parameters.

INPUT CONVENTION
----------------
  log_L  : log10(L / L_sun)
  T_eff  : effective temperature in Kelvin
  log_g  : standard spectroscopic surface gravity in CGS units, i.e. log10(g)
            where g is in cm s⁻². Values ~0.2–2.0 are typical for RGB/AGB.
            NOT log(g/g_sun) in a dimensionless sense.

KEY SCALING RELATIONS USED
---------------------------
  ν_max   ∝  g / sqrt(T_eff)            [Brown et al. 1991; Kjeldsen & Bedding 1995]
  ν_ac    ∝  g / sqrt(T_eff)            [Chaplin & Miglio 2013; Belkacem et al. 2011]
  Δν      ∝  sqrt(ρ̄)                    [Ulrich 1986]

  ρ̄ from Stefan-Boltzmann + g = GM/R²:
    R / R_sun  =  sqrt(L/L_sun) * (T_sun / T_eff)^2
    M / M_sun  =  (g/g_sun) * (R/R_sun)^2
    ρ̄ / ρ̄_sun  =  (M/M_sun) / (R/R_sun)^3  =  (g/g_sun) / (R/R_sun)

REFERENCES
──────────────────────────────────────────────────────────────────────────────
  [1]  Ulrich, R. K. (1986), ApJ, 306, L37
       Δν ∝ sqrt(mean density) — the first Δν scaling relation.

  [2]  Brown, T. M., et al. (1991), ApJ, 368, 599
       Showed ν_max is proportional to the acoustic cut-off frequency,
       hence ν_max ∝ g / sqrt(T_eff).

  [3]  Kjeldsen, H. & Bedding, T. R. (1995), A&A, 293, 87
       Practical scaling relations: ν_max ∝ g / sqrt(T_eff),
       Δν ∝ sqrt(M/R³). Solar reference values established.

  [4]  Chaplin, W. J. & Miglio, A. (2013), ARA&A, 51, 353
       Comprehensive review of all asteroseismic scaling relations,
       theoretical underpinning, and observational validation.

  [5]  Belkacem, K., et al. (2011), A&A, 530, A142
       Physical explanation for why ν_max traces the acoustic cut-off:
       the turbulent driving rate peaks at the cut-off frequency.

  [6]  Balmforth, N. J. & Gough, D. O. (1990), ApJ, 362, 256
       Theoretical derivation of the acoustic cut-off in a stratified
       stellar atmosphere: ω_ac = γ g / (2 c_s).

  Solar reference values used here:
    T_sun    = 5778  K
    g_sun    = 27420 cm s⁻²  →  log g_sun = 4.438
    ρ̄_sun   = 1408  kg m⁻³
    ν_max,sun = 3090  μHz      (Kjeldsen & Bedding 1995)
    ν_ac,sun  = 5300  μHz      (Chaplin & Miglio 2013 Table 1)
    Δν_sun   = 134.9 μHz      (Chaplin & Miglio 2013 Table 1)
──────────────────────────────────────────────────────────────────────────────
"""

import numpy as np


# ── Physical constants (SI) ────────────────────────────────────────────────
_k_B  = 1.380649e-23   # Boltzmann constant   [J K⁻¹]
_m_H  = 1.67262192e-27 # proton mass          [kg]
_G    = 6.67430e-11    # Gravitational const  [m³ kg⁻¹ s⁻²]

# ── Solar constants ────────────────────────────────────────────────────────
_T_SUN       = 5778.0   # K
_LOG_G_SUN   = 4.438    # log10(g_sun / [cm s⁻²]);  g_sun = 27420 cm s⁻²
_RHO_SUN     = 1408.0   # kg m⁻³  (mean density)
_NU_MAX_SUN  = 3090.0   # μHz
_NU_AC_SUN   = 5300.0   # μHz
_DELTA_NU_SUN = 134.9   # μHz
_R_SUN       = 6.957e8  # m
_M_SUN       = 1.989e30 # kg

# ── Photosphere defaults ───────────────────────────────────────────────────
_GAMMA = 5.0 / 3.0  # adiabatic index (mono-atomic ideal gas)
_MU    = 1.3        # mean molecular weight for a cool (~4000 K),
                    # largely neutral, H+He photosphere


def stellar_oscillation_params(
    log_L: float,
    T_eff: float,
    log_g: float,
    gamma: float = _GAMMA,
    mu:    float = _MU,
    T_sun: float = _T_SUN,
    log_g_sun:    float = _LOG_G_SUN,
    rho_sun:      float = _RHO_SUN,
    nu_max_sun:   float = _NU_MAX_SUN,
    nu_ac_sun:    float = _NU_AC_SUN,
    delta_nu_sun: float = _DELTA_NU_SUN,
) -> dict:
    """
    Compute asteroseismic and structural parameters for a red giant star.

    Parameters
    ----------
    log_L : float
        log10(L / L_sun)
    T_eff : float
        Effective temperature  [K]
    log_g : float
        Spectroscopic surface gravity, log10(g / [cm s⁻²]).
        Typical range for RGB/AGB giants: 0.2 – 2.0.
    gamma : float, optional
        Adiabatic index. Default 5/3 (ideal mono-atomic gas).
    mu : float, optional
        Mean molecular weight of the photospheric gas.
        Default 1.3 (cool, largely neutral H + He atmosphere).
    T_sun, log_g_sun, rho_sun, nu_max_sun, nu_ac_sun, delta_nu_sun
        Solar reference values — see module docstring.

    Returns
    -------
    dict
        R_Rsun           : radius                          [R_sun]
        M_Msun           : mass inferred from g and R      [M_sun]
        rho_mean_gcm3    : mean density                    [g cm⁻³]
        rho_mean_solar   : mean density                    [ρ_sun]
        nu_ac_uHz        : acoustic cut-off freq (scaling) [μHz]
        nu_ac_phys_uHz   : acoustic cut-off freq (ab initio, Balmforth & Gough 1990) [μHz]
        P_ac_days        : period at ν_ac (scaling)        [days]
        nu_max_uHz       : ν_max from scaling              [μHz]
        P_max_days       : period at ν_max                 [days]
        delta_nu_uHz     : large frequency separation      [μHz]
        log_g_sun_used   : solar log g used as reference

    Notes
    -----
    The two acoustic cut-off estimates (scaling vs ab initio) should agree
    to within ~10 %.  The scaling version is the observationally calibrated
    one recommended for comparison with data.
    """

    # ── Dimensionless ratios ──────────────────────────────────────────────
    L_ratio = 10.0 ** log_L                   # L / L_sun
    g_cgs   = 10.0 ** log_g                   # g  [cm s⁻²]
    g_sun_cgs = 10.0 ** log_g_sun             # g_sun [cm s⁻²]
    g_ratio = g_cgs / g_sun_cgs               # g / g_sun  (dimensionless)
    T_ratio = T_eff / T_sun                   # T_eff / T_sun

    # ── Radius  (Stefan-Boltzmann) ────────────────────────────────────────
    # L = 4π R² σ T⁴  →  R/R_sun = sqrt(L/L_sun) * (T_sun/T_eff)²
    R_ratio = np.sqrt(L_ratio) * (T_sun / T_eff) ** 2   # R / R_sun

    # ── Mass  (from g = GM/R²) ────────────────────────────────────────────
    # M/M_sun = (g/g_sun) * (R/R_sun)²
    M_ratio = g_ratio * R_ratio ** 2          # M / M_sun

    # ── Mean density ──────────────────────────────────────────────────────
    # ρ̄/ρ̄_sun = (M/M_sun) / (R/R_sun)³ = (g/g_sun) / (R/R_sun)
    rho_ratio = g_ratio / R_ratio             # ρ̄ / ρ̄_sun
    rho_SI    = rho_ratio * rho_sun           # kg m⁻³
    rho_gcm3  = rho_SI * 1e-3                 # g cm⁻³

    # ── ν_max  (scaling relation) ─────────────────────────────────────────
    # ν_max / ν_max_sun = (g/g_sun) * (T_eff/T_sun)^{-1/2}
    # Refs: Brown+1991; Kjeldsen & Bedding 1995
    nu_max   = nu_max_sun * g_ratio * T_ratio ** (-0.5)   # μHz
    P_max_s  = 1.0 / (nu_max * 1e-6)                      # seconds
    P_max_d  = P_max_s / 86400.0                          # days

    # ── ν_ac  (scaling relation, same proportionality as ν_max) ──────────
    # ν_ac / ν_ac_sun = (g/g_sun) * (T_eff/T_sun)^{-1/2}
    # Refs: Chaplin & Miglio 2013; Belkacem+2011
    nu_ac   = nu_ac_sun * g_ratio * T_ratio ** (-0.5)     # μHz
    P_ac_s  = 1.0 / (nu_ac * 1e-6)                        # seconds
    P_ac_d  = P_ac_s / 86400.0                            # days

    # ── ν_ac  (ab initio, from first principles) ──────────────────────────
    # ω_ac = γ g / (2 c_s)  with  c_s = sqrt(γ k_B T / (μ m_H))
    # → ν_ac = γ g / (4π c_s)      [Balmforth & Gough 1990]
    g_SI = g_cgs * 1e-2                                   # m s⁻²
    c_s  = np.sqrt(gamma * _k_B * T_eff / (mu * _m_H))   # m s⁻¹
    nu_ac_phys = (gamma * g_SI / (4.0 * np.pi * c_s)) * 1e6   # μHz

    # ── Δν  (large frequency separation) ─────────────────────────────────
    # Δν / Δν_sun = sqrt(ρ̄ / ρ̄_sun)      [Ulrich 1986]
    delta_nu = delta_nu_sun * np.sqrt(rho_ratio)          # μHz

    return {
        "R_Rsun":          R_ratio,
        "M_Msun":          M_ratio,
        "rho_mean_gcm3":   rho_gcm3,
        "rho_mean_solar":  rho_ratio,
        "nu_ac_uHz":       nu_ac,
        "nu_ac_phys_uHz":  nu_ac_phys,
        "P_ac_days":       P_ac_d,
        "nu_max_uHz":      nu_max,
        "P_max_days":      P_max_d,
        "delta_nu_uHz":    delta_nu,
        "log_g_sun_used":  log_g_sun,
    }


# ── Convenience: vectorised wrapper ───────────────────────────────────────
def stellar_oscillation_table(log_L_arr, T_eff_arr, log_g_arr, **kwargs) -> dict:
    """
    Vectorised wrapper — accepts array inputs, returns dict of arrays.

    Parameters
    ----------
    log_L_arr, T_eff_arr, log_g_arr : array-like
        Same meaning as in stellar_oscillation_params().
    **kwargs : passed through to stellar_oscillation_params()

    Returns
    -------
    dict of numpy arrays, one entry per star.
    """
    log_L_arr = np.asarray(log_L_arr, dtype=float)
    T_eff_arr = np.asarray(T_eff_arr, dtype=float)
    log_g_arr = np.asarray(log_g_arr, dtype=float)

    results = [
        stellar_oscillation_params(lL, T, lg, **kwargs)
        for lL, T, lg in zip(log_L_arr, T_eff_arr, log_g_arr)
    ]
    keys = results[0].keys()
    return {k: np.array([r[k] for r in results]) for k in keys}


# ── Quick demo ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 62)
    print("Hypothetical case 2 stars different masses and teff but similar luminosity")
    print("=" * 62)

    cases = [
        ("AGB  0.6 M☉", np.log10(150), 4300, None),
        ("RGB  0.8 M☉", np.log10(150), 4100, None),
    ]

    # For the hypothetical case we don't have log_g directly,
    # so we derive it from the stated masses and L/T.
    for label, lL, Teff, _ in cases:
        L_ratio  = 10 ** lL
        R_ratio  = np.sqrt(L_ratio) * (5778 / Teff) ** 2
        M_Msun   = 0.6 if "AGB" in label else 0.8
        g_cgs    = M_Msun * 10 ** 4.438 / R_ratio ** 2   # g in cm/s²
        lg       = np.log10(g_cgs)
        res = stellar_oscillation_params(lL, Teff, lg)
        print(f"\n{label}  |  log_g = {lg:.3f}")
        print(f"  R = {res['R_Rsun']:.1f} R_sun     M = {res['M_Msun']:.3f} M_sun")
        print(f"  ρ̄ = {res['rho_mean_gcm3']:.4e} g/cm³  ({res['rho_mean_solar']:.4e} ρ_sun)")
        print(f"  ν_ac  = {res['nu_ac_uHz']:.3f} μHz   →  P_ac  = {res['P_ac_days']:.2f} days")
        print(f"  ν_max = {res['nu_max_uHz']:.3f} μHz   →  P_max = {res['P_max_days']:.2f} days")
        print(f"  Δν    = {res['delta_nu_uHz']:.4f} μHz")

    print("\n" + "=" * 62)
    print("47 Tuc sample — first and last 3 stars from table")
    print("=" * 62)

    # Subset of the table  [log_L, T_eff, log_g]
    sample = [
        ("MED_RGB_3872",  2.12, 4483, 1.81),
        ("MED_RGB_557",   2.14, 4500, 1.80),
        ("MED_RGB_459",   2.14, 4469, 1.79),
        ("MED_VAR_1399",  3.27, 3532, 0.23),
        ("MED_VAR_3405",  3.28, 3596, 0.24),
        ("MED_VAR_1373",  3.29, 3493, 0.19),
    ]

    hdr = f"{'Star':<20} {'R/R☉':>6} {'M/M☉':>6} {'ρ [g/cm³]':>12} "
    hdr += f"{'ν_ac [μHz]':>11} {'P_ac [d]':>9} {'ν_max [μHz]':>12} {'Δν [μHz]':>9}"
    print(hdr)
    print("-" * len(hdr))

    for name, lL, Teff, lg in sample:
        r = stellar_oscillation_params(lL, Teff, lg)
        print(
            f"{name:<20} "
            f"{r['R_Rsun']:6.1f} "
            f"{r['M_Msun']:6.3f} "
            f"{r['rho_mean_gcm3']:12.3e} "
            f"{r['nu_ac_uHz']:11.3f} "
            f"{r['P_ac_days']:9.2f} "
            f"{r['nu_max_uHz']:12.3f} "
            f"{r['delta_nu_uHz']:9.4f}"
        )
