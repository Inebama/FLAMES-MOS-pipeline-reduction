# 11 — Manual mode: the `Fast_usage/` tools

[← back to index](README.md)

The `Fast_usage/` folder contains five stand-alone scripts. The first two show how to drive
the pipeline's two engines (**cross-correlation** and **continuum fitting**) on arbitrary
spectra, with no `INPUT.py`, no `MAIN.py` and no Observing-Block structure. The other three
are independent physics/astrometry utilities that were developed alongside the pipeline.

All of them start with the same two lines that make the package importable from the
sub-folder:

```python
import sys
path_package = __file__.replace('Fast_usage/<script name>.py', '')
sys.path.append(path_package)
```

> The example paths inside these scripts point at the author's own folders
> (`TEST_OUTPUT/Test_data_ian/`, `GIRAFFE_TEMPLATE_HR14b_R37000.dat`, `/home/ian/...`) which do
> not all exist in the distributed repository. **Edit the paths before running.** The bundled
> equivalents are `TESTING_FOLDER/TEST_OUTPUT/Test_example_try_to_reproduce_it/` and
> `GIRAFFE_TEMPLATE_HR14b_R36900.dat`.

---

## 1. `Fast_RV.py` — cross-correlation on arbitrary spectra

A ~70-line template showing the minimal RV measurement.

```python
from CrossCorrelation import CCF
from User_Functions import TEMPLATE_READER

Template_wave, Template_flux = TEMPLATE_READER(complete_path_template)
RADVEL_RANGE = (-120, 120)

def get_w_f(path, MIN=6391.4, MAX=6613):       # ← adapt to YOUR file layout
    with fits.open(path) as hdul:
        OBJ  = hdul[0].header['OBJECT']
        RVC  = hdul[0].header['RVBARY']
        wave = hdul[1].data['wave']
        mask = (MIN < wave) & (wave < MAX)
        flux = hdul[1].data['flux'][mask]
    return wave[mask], flux, OBJ, MAG, RVC, RA, DEC

def get_rv(path):
    wave, flux, OBJ, MAG, RVC, RA, DEC = get_w_f(path)
    myCCF = CCF(template_wavelength = Template_wave,
                template_flux       = Template_flux,
                observed_wavelength = wave,
                observed_flux       = flux,
                IGNORE   = False,
                PRINT    = True,
                RV_range = RADVEL_RANGE)
    RV_sample, Corr, RV, ERR = myCCF.find_rv(steps3=500)
    myCCF.get_errs_info(PRINT=True)     # full error breakdown on the terminal
    myCCF.PLOT()                        # five-panel diagnostic figure

for path in path_stars:
    get_rv(path)
```

To use it on your own data you only need to rewrite `get_w_f` so that it returns the observed
wavelength and flux; everything else is generic. The class handles the flux scaling, the
coverage checks and the trimming itself. See
[08 — Stage 3](08_stage3_radial_velocities.md) for the meaning of every argument.

**Minimal version**, if you already have arrays in memory:

```python
from CrossCorrelation import CCF
myCCF = CCF(tw, tf, ow, of, RV_range=(-200, 200))
rv_grid, ccf, RV, ERR = myCCF.find_rv()
print(f'RV = {RV:.3f} ± {ERR:.3f} km/s')
```

## 2. `Fast_continuum.py` — the continuum GUI on arbitrary spectra

```python
from Continuum_fitting import cont_fit

PATH  = '/where/my/spectra/live/'
Files = ['star1.fits', 'star2.fits', ...]        # names only, not full paths

WAVE_RANGE = (6391.4, 6613)      # explicit window, or a float in (0,1) to trim the edges

def read_fits(hdul):             # receives an OPEN HDUList, returns (wave, flux)
    try:
        return hdul[-1].data.wave, hdul[-1].data.flux
    except:
        return hdul[-2].data.wave, hdul[-2].data.flux

CONT = cont_fit(Files, PATH,
                wave_range    = WAVE_RANGE,
                Continue      = True,        # resume from Stored_params.npy in PATH
                READ_FUNCTION = read_fits)
CONT.store_cont_sub(Force=False)             # writes PATH/Continuum_sub/
```

Everything about the GUI, the four sliders and the iteration semantics is documented in
[09 — Stage 4](09_stage4_continuum.md). Remember: **FITS input only**, and closing the window
is how you end the session.

For the unsupervised variant on an arbitrary folder:

```python
from automated_cont_fit import multiple_spectra_cont
multiple_spectra_cont(spectra=Files, path=PATH, wave_range=0.05)   # → PATH/CONT_SUB/CS-*.fits
```

---

## 3. `barycentric_vels.py` — barycentric velocity explorer

Independent of the pipeline. Resolves a target by name in SIMBAD (needs `astroquery`) and
computes the barycentric velocity correction with `astropy`.

### `plot_barycentric_velocity(...)`

```python
plot_barycentric_velocity(target_name = "47 Tuc",
                          mjd_range   = (60000, 60365),   # or mjd_array = [...]
                          n_points    = 1200,
                          observatory = "paranal",        # paranal | lasilla | ctio | campana
                          savefile    = "bary.txt",       # optional ASCII output
                          return_values = True)           # → (mjd, velocities)
```

Plots the correction over the requested time span and optionally writes a two-column ASCII
file with a self-describing header (target, ICRS coordinates, observatory).

Use it to **plan observations** (when is the correction stationary? when does it change
fastest?) or to **check** the values the pipeline wrote into `BARY0…BARYn`.

### `compute_components(...)` and `components_plot(...)`

The interesting part. They decompose the correction into

* the **Earth's orbital** motion projected on the line of sight (~±30 km s⁻¹ over a year), and
* the **observatory's rotational** motion (~±0.46 km s⁻¹ over a day),

and compare three treatments of the same quantity:

| Curve | What it is |
|---|---|
| classical, ICRS direction | velocities projected on a *fixed* unit vector towards the target |
| classical, GCRS direction | projected on the *time-dependent* apparent direction |
| `astropy` full correction | `SkyCoord.radial_velocity_correction`, including aberration, precession/nutation and relativistic terms |

The residual panels quantify what the classical approximations miss (typically tens of m s⁻¹),
which is exactly the level at which you need to care for precise RV work. A peak-to-peak
summary is printed at the end.

```python
components_plot(target_name="47 Tuc", observatory="Paranal",
                mjd_start=60000.0, mjd_end=60365, n_points=4000)
```

Note that the observatory list is hard-coded (Paranal, La Silla, CTIO, Las Campanas); add
entries to the `observatories` dictionary using `EarthLocation.of_site` or explicit
coordinates.

---

## 4. `stellar_oscillation_params.py` — asteroseismic scaling relations

Pure-numpy, no I/O, fully referenced in its module docstring. Given
`log L/L☉`, `T_eff` [K] and `log g` [cgs], it returns:

| Output | Relation |
|---|---|
| `R_Rsun` | Stefan–Boltzmann: `R/R☉ = √(L/L☉)·(T☉/T_eff)²` |
| `M_Msun` | from `g = GM/R²`: `M/M☉ = (g/g☉)·(R/R☉)²` |
| `rho_mean_gcm3`, `rho_mean_solar` | `ρ̄/ρ̄☉ = (g/g☉)/(R/R☉)` |
| `nu_max_uHz`, `P_max_days` | `ν_max ∝ g/√T_eff`, scaled to `ν_max,☉ = 3090 μHz` (Brown et al. 1991; Kjeldsen & Bedding 1995) |
| `nu_ac_uHz`, `P_ac_days` | acoustic cut-off with the same scaling, `ν_ac,☉ = 5300 μHz` (Chaplin & Miglio 2013; Belkacem et al. 2011) |
| `nu_ac_phys_uHz` | *ab initio* cut-off `ν_ac = γg/(4πc_s)` with `c_s = √(γk_BT/μm_H)` (Balmforth & Gough 1990) — a consistency check on the scaling version (agreement expected to ~10 %) |
| `delta_nu_uHz` | `Δν ∝ √ρ̄`, `Δν☉ = 134.9 μHz` (Ulrich 1986) |

```python
from stellar_oscillation_params import stellar_oscillation_params, stellar_oscillation_table

r = stellar_oscillation_params(log_L=2.12, T_eff=4483, log_g=1.81)
print(r['P_max_days'], r['R_Rsun'], r['M_Msun'])

tab = stellar_oscillation_table(logL_array, teff_array, logg_array)   # dict of arrays
```

Defaults: `gamma = 5/3`, `mu = 1.3` (cool, largely neutral H+He photosphere). Solar reference
values are module constants and can be overridden per call.

**Important convention:** `log_g` is the *spectroscopic* `log₁₀(g / cm s⁻²)` (≈ 0.2–2.0 for
RGB/AGB giants), **not** `log(g/g☉)`.

Running the module directly prints a demo: two hypothetical giants of equal luminosity but
different mass/T_eff, and a table for a sample of 47 Tuc stars.

## 5. `convective_turnover.py` — convective turnover times

Companion to the previous module; same input convention. Computes two quantities:

**Local (photospheric)** `τ_local = H_p / v_c`, where

* `H_p = k_B T/(μ m_H g)` is the pressure scale height;
* `ρ_phot` comes from a closed-form combination of an **Eddington grey atmosphere**
  (`P = 2g/3κ`), the **ideal gas law** and an **H⁻ opacity power law** (`κ ∝ Z_rel ρ^½ T⁹`),
  with the coefficient calibrated so that the solar photosphere reproduces
  `ρ_phot,☉ = 1.25×10⁻⁷ g cm⁻³`;
* `v_c` is the mixing-length velocity in the **efficient-convection limit**,
  `v_c³ = (α/2)²·2g H_p F_c/(ρ c_p T)` (Böhm-Vitense 1958; Kippenhahn, Weigert & Weiss 2012 §7.3).

**Global (envelope-integrated)** via the Gunn et al. (1998) scaling calibrated to
`τ_☉,global ≈ 12 d` (Corsaro et al. 2021):

```
τ_global/τ_☉ = (L/L☉)^(−1/3) · (M/M☉)^(1/3) · (R/R☉)²
```

```python
from convective_turnover import convective_turnover_full, convective_turnover_table

r = convective_turnover_full(log_L=3.27, T_eff=3545, log_g=0.24)
r['tau_local_d'], r['tau_global_d'], r['v_conv_kms'], r['ratio_Pmax_tau']
```

The convenience wrapper also returns `P_max_days` and the diagnostic ratio
`P_max/τ_local`:

| `P_max/τ_local` | Regime |
|---|---|
| ≫ 1 | stochastic, solar-like oscillations |
| ~ 1 | transition regime (semi-regular pulsators) |
| < 1 | convective damping dominates |

Default parameters: `alpha_MLT = 1.6`, `gamma = 5/3`, `mu = 1.3`, `Z_rel = 0.2` (i.e.
`[Fe/H] ≈ −0.7`, tuned for 47 Tuc). Override them per call for other systems.

Running the module directly prints a **solar self-check** (ρ_phot, κ, v_conv, τ_local against
their expected values) followed by a table for the 47 Tuc sample — run it once after any edit
to confirm the calibration still holds.

> Both modules carry the note *"This code was commented and edited by AI"* at the top of their
> docstrings, and `convective_turnover.py` explicitly warns that the global scaling can be off
> by a factor 2–5 depending on the assumed core mass. They are order-of-magnitude tools for
> interpretation, not substitutes for stellar-structure models.
