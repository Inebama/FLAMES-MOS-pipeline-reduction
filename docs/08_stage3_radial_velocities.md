# 08 — Stage 3: radial velocities by cross-correlation

[← back to index](README.md)

**Modules:** `CrossCorrelation.py` (the `CCF` class — the engine) and `Getting_RV.py` (the
driver script)
**Driver:** the `if RADIAL_VEL …` branches of `MAIN.py`
**Products:** `OB-n/RESTFRAME/*.fits`, `OUTPUT_RV.txt`, `CCFS.png`, and the `RV`/`RVERR`/
`RVTRUE` keywords added to the input spectra

---

## 1. How the driver and the engine fit together

`Getting_RV.py` is **not imported** — `MAIN.py` reads it and `exec`s it inside its own
namespace:

```python
with open(path_package+'Getting_RV.py') as f:
    exec(f.read())
```

so `Getting_RV.py` sees every variable from `INPUT.py` plus `path_folder` and `spectrums`
(the OB folder and its file list). This is why the script uses bare globals rather than
function arguments, and why you can freely edit its top block (`method`, `SAVE_NAME`,
`steps3=500`, …) without touching anything else.

Per spectrum the driver:

1. reads `wave`, `flux` and the metadata (`get_w_f`),
2. builds a `CCF` object,
3. calls `find_rv(steps3=500)`,
4. prints one line, writes one line to `OUTPUT_RV.txt`,
5. stamps `RV`, `RVERR`, `RVTRUE` into the parent file (`UPDATE_header`),
6. writes a rest-frame copy into `RESTFRAME/` (`new_fits`).

## 2. Preparing the observed spectrum (`get_w_f`)

```python
wave = hdul[1].data['wave']
mask = (MIN < wave) & (wave < MAX)
flux = hdul[1].data['flux'][mask]
if not (0.5 < median(flux) < 1.5):
    flux /= median(flux)
```

* `MIN`/`MAX` come from `MIN_WAVE_CCF`/`MAX_WAVE_CCF` and are in **Ångström**. They remove the
  low-efficiency CCD edges and the zero-padded pixels Stage 2 leaves at both ends.
  **Set both or neither:** if either is `None` the code substitutes an automatic trim of 1/20 of
  the wavelength range at each end, and a lone value is not used. On a merged HR14b spectrum
  the automatic window is 6394.200 – 6614.840 Å; the HR14b preset gives 6391.400 – 6612.990 Å.
* **The trim affects only the CCF.** `new_fits` re-reads the parent file, so the rest-frame
  product keeps the full, untrimmed wavelength range.
* The flux is rescaled to a median near 1 if it is not already, because the CCF is a plain
  (un-normalised) product sum and therefore sensitive to the flux scale.
* `OBJECT`, `RVBARY`, `RA`, `DEC`, `MAG` are read from the primary header at the same time.

## 3. The `CCF` class

```python
CCF(template_wavelength, template_flux,
    observed_wavelength, observed_flux,
    RV_range = None,        # (min, max) km/s
    IGNORE   = False,       # skip all consistency checks (requires RV_range)
    PRINT    = False,       # verbose checks
    Cosmic_rays = True,     # clip the top 1 % of the flux
    Check_template_relative_to_obs = False)   # plot max blue/red-shifted template vs spectrum
```

### 3.1 Consistency checks (`IGNORE = False`)

**Flux scale.** Template and observed flux are each divided by their median if that median is
outside `[0.5, 1.5]`, with a message on the terminal.

**Wavelength coverage.** The template must cover the observed range *after* being shifted by
the extreme RVs to be tested. Two cases:

* **`RV_range = None`** → the class *derives* the largest testable range from the actual
  coverage,

  ```
  RV_min = c·(λ_max,obs/λ_max,temp − 1)      # blue shift needed to reach the red edge
  RV_max = c·(λ_min,obs/λ_min,temp − 1)      # red  shift needed to reach the blue edge
  ```

  and refuses to continue if the template is narrower than the spectrum. It also sets
  `self.Computed = True`, which forces the iterative search (§4).

* **`RV_range` given** → the class computes the wavelengths the shifted template will span,
  ```
  λ_min,eff = λ_min,temp·(1 + RV_max/c)      λ_max,eff = λ_max,temp·(1 + RV_min/c)
  ```
  and, if the observed spectrum sticks out, **trims the observed spectrum** with the message
  ```
  In order to test RV=[…] I'm gona cut the obs spectra since your template is to short in range of wavelength
  ```
  Conversely, if the template is much wider than needed it is trimmed for speed.

  This is deliberate: it guarantees that **the same set of pixels is compared at every trial
  velocity**, so the CCF is not biased by pixels entering and leaving the sum.

**Ordering.** `RV_range` must be `[min, max]` with `min < max`, else an exception is raised.

### 3.2 Cosmic-ray / emission clipping

```python
p99 = np.percentile(obs_flux, 99)
obs_flux[obs_flux > p99] = p99
```

Simple and effective for absorption-line work: it removes cosmic rays and hot pixels that
would otherwise dominate a product-sum correlation. **It also flattens genuine emission
lines** — if you are correlating on emission features, pass `Cosmic_rays=False`.

### 3.3 The correlation itself

```python
def cross_correlation(self, RV):
    shifted = template_wavelength * (1 + RV/299792.458)
    synth   = interp1d(shifted, template_flux, bounds_error=False, fill_value=0)(observed_wavelength)
    return np.sum(observed_flux * synth)
```

i.e.

```math
\mathrm{CCF}(v)=\sum_i S_{\rm obs}(\lambda_i)\; S_{\rm tmpl}\!\left(\lambda_i - \Delta\lambda(v)\right)
```

a straight product sum with the template interpolated onto the observed pixel grid — no
normalisation by the template norm, no continuum removal (the code contains a commented-out
"template-match" variant that divides by `Σ|obs − synth|`).

Because there is no normalisation, the **flux scaling of §2/§3.1 matters**, and the absolute
CCF value is meaningless; only the position of the maximum is used, and the CCF is rescaled to
`[0, 1]` before any fitting.

## 4. The three-pass search (`find_rv`)

```python
find_rv(steps3=300, steps1=100, steps2=100, Tolerance=90, Manual=False)
```

`Getting_RV.py` calls it with `steps3=500`.

| Pass | Range | Samples | Purpose |
|---|---|---|---|
| 1 | the full `RV_range` | `steps1` = 100 | locate the global maximum coarsely |
| 2 | the width of the pass-1 peak | `steps2` = 100 | resolve the peak, expose nearby secondary peaks or asymmetries |
| 3 | the width of the pass-2 peak | `steps3` = 500 | final, high-resolution sampling — the reported RV is `argmax` of this pass |

The three samplings are concatenated into a single, globally ascending `RV_sample`/
`Correlation` pair (coarse in the wings, dense in the core), which is what gets plotted and
returned. If anything goes wrong in the nested passes, the code falls back to a single uniform
scan with `steps1+steps2+steps3` points.

Setting `Manual=True` **and** providing an explicit `RV_range` performs a single uniform scan
with `steps3` points instead.

### The peak finder

`find_peaks(Y, npeak=1, peak='max', tolerance=60, SHOW=False)` is a small hand-written routine:
it takes the global extremum, then walks outwards on both sides until the curve drops below

```
LIM = Y_min + (Y_max − Y_min) · tolerance/100
```

or the array ends, returning the peak index and its `[left, right]` bounds. `find_rv` uses
`tolerance=90` for the final peak, i.e. only the top 10 % of the CCF height is used for the
analytic fits. Passing `SHOW=True` draws a small explanatory figure — useful for
understanding or debugging the routine.

## 5. The error budget

The CCF is first rescaled to `[0, 1]`:

```python
Correlation = (Correlation − min) / (max − min)
```

Four analytic functions are then fitted (`scipy.optimize.curve_fit`) to the **peak region only**:

| Fit | Model | Free parameters |
|---|---|---|
| Fixed Gaussian | `exp(−(x−x₀)²/2σ²)` | `x₀, σ` (height fixed at 1) |
| Fixed parabola | `−a(x−x₀)² + 1` | `x₀, a` (height fixed at 1) |
| Gaussian | `A·exp(−(x−x₀)²/2σ²)` | `x₀, σ, A` |
| Parabola | `−a(x−x₀)² + c` | `x₀, a, c` |

Fits that do not converge are recorded as `nan` (and a `WARNING: Some of the fits didn't
converge` is printed) but do not stop the run.

The final uncertainty combines three independent terms in quadrature:

```math
\sigma_{\rm RV}=\sqrt{\;\sigma^2_{\rm centres}\;+\;\widetilde{\sigma^2_{\rm fit}}\;+\;\sigma^2_{\rm pix}\;}
```

| Term | Code | Meaning |
|---|---|---|
| `σ_centres` | `np.nanstd([RV, c₁, c₂, c₃, c₄])` | scatter of the five estimates of the peak position (the sampled maximum plus the four fitted centres) — captures the systematic effect of the assumed peak shape |
| `σ_fit` | `np.nanmedian([e₁, e₂, e₃, e₄])` | median of the formal errors on the centre from the covariance matrices |
| `σ_pix` | `ERR_init` | the sampling floor, computed at construction: `median(Δλ)/median(λ) · c/2` — the RV corresponding to half a pixel at the central wavelength |

`get_errs_info(PRINT=True)` prints the full breakdown; `get_errs_info(ALL_info=True)` returns
`(RV, ERR, median_of_centres, std_of_centres, median_of_fit_errors)`, which `find_rv` also
stores as `self.RV_error`, `self.RV_median`, `self.RV_std`, `self.median_err_fits`.

> **In practice `σ_pix` dominates for high-S/N spectra.** For the bundled GIRAFFE test data
> `σ_pix = 0.3840 km s⁻¹`, and the reported `RVERR` is `0.38404` — the fit terms contribute at
> the 10⁻³ level.

### Why `σ_pix` is evaluated on the combined grid

`σ_pix` is computed from the spectrum the CCF actually sees, i.e. the **Stage-2 output grid**,
which is ~3× finer than the native detector pixels. That is deliberate, and it reflects real
information rather than an accounting trick: the exposures merged into one OB are taken at
different times, so each is shifted by a *different* barycentric velocity before co-addition.
Their pixel grids therefore land at different sub-pixel positions on the stellar line profile,
and the combined spectrum samples that profile more finely than any single exposure — the
classical dithered/differential-sampling gain. Resampling onto a ~3× grid preserves that gain
instead of discarding it at the coarsest common step.

Concretely, for the bundled GIRAFFE data the combined grid gives `σ_pix = 0.384 km s⁻¹`, where
a raw 0.05 Å pixel at 6500 Å would imply ≈ 1.15 km s⁻¹.

Two practical consequences:

* The gain is real only to the extent that the epochs in the OB genuinely span different
  barycentric velocities. A block of consecutive exposures taken within a few minutes shares
  essentially the same correction, so the sub-pixel diversity — and hence the justification for
  the finer grid — is much weaker. Short-cadence OBs are the case to be careful about.
* If your data do not justify it, the oversampling factor is a single number in
  `Merging_and_substract_sky.MERGING_MULTIPLE_SPECTRA`:

  ```python
  samp_wave = extend_linear_range(WMIN*0.9999999, WMAX*1.000000001, 3 * N_pixels)
  #                                                                 ^ change this factor
  ```

  Setting it to `1 * N_pixels` keeps the native sampling and raises `RVERR` accordingly. The
  measured RV itself is essentially unchanged either way; only the quoted floor moves.

Note also that adjacent pixels of the oversampled grid are correlated by construction, which
matters if you compute χ², S/N per pixel or equivalent widths from these products.

## 6. Sign conventions

This trips everyone up once, so here it is explicitly.

* Stage 2 already multiplied the wavelengths by `(1 + v_bary/c)`. The spectra therefore live
  in the **barycentric frame**.
* The CCF shifts the *template* by `+RV` until it matches the observed spectrum, so a positive
  `RV` means the star's lines are redshifted with respect to the template.
* Consequently, `RV` as returned by `find_rv` **is the barycentric radial velocity**.

What gets written where:

| Destination | Value | Interpretation |
|---|---|---|
| `OUTPUT_RV.txt` column `RV_corrected` (7th) | `RV` | ✅ the physical (barycentric) radial velocity — **use this one** |
| FITS header `RVTRUE` | `RV` | same |
| `OUTPUT_RV.txt` column `RV` (4th) | `RV − RVBARY` | the value before removing the Earth's motion |
| FITS header `RV` | `RV − RVBARY` | same |
| FITS header `RVERR` | `ERR` | the uncertainty, identical for both |

The rest-frame spectra are built with the barycentric value:
`wave_restframe = wave / (1 + RV/c)`.

`plotingRVS.py` reads column −2 of `OUTPUT_RV.txt`, i.e. `RV_corrected` — consistent.

> The header comment on `RVTRUE` reads *"(Measured RV+RVBARY)"*, which does not match the
> arithmetic actually performed (`RVTRUE = RV`, `RV = RVTRUE − RVBARY`). Trust the table above.

## 7. The `'Template'` method

At the top of `Getting_RV.py`:

```python
method = 'Synthetic'      # or 'Template'
```

* **`'Synthetic'`** (default, and the only mode wired into `INPUT.py`): every star is
  correlated against the synthetic template.
* **`'Template'`**: one real star (`path_temp_star`, which you must define in the script) is
  first measured against the synthetic template and shifted to rest; that observed, rest-framed
  spectrum then becomes the template for all the others, and `MIN`/`MAX` are tightened by
  10 Å at each end to keep the interpolation valid. This can reduce template-mismatch
  systematics (the relative velocities become differential measurements) at the cost of tying
  the zero point to one star. `OUTPUT_RV.txt` records the template star's own RV in a comment
  line and labels each measurement `relative_to_Temp_Star_<name>`.

## 8. Diagnostic plots

### `CCF.PLOT(...)` — one spectrum, full diagnostics (`PLOT_EACH_ONE = True`)

Five panels:

| Panel | Content |
|---|---|
| top left | the full CCF vs trial RV; the selected peak region marked, the adopted RV as a dotted line |
| top right | zoom on the peak with the **four analytic fits** overplotted, plus vertical lines at the best RV (grey) and the median of the centres (orange) |
| middle | full observed spectrum (red, shifted to rest) over the template (blue-grey) |
| bottom left / bottom right | two zoom windows on the same comparison, connected to the middle panel with inset-zoom indicators |

Arguments: `RV`, `window`, `window1..3` (axis ranges), `zoom`, `PAUSE`/`time` (flash the figure
for *n* seconds — used when looping), `Save`, `Save_name`.

`plot_fit_peak(ax)` draws just the top-right panel and can be reused on your own axes.

### `CCFS.png` — one figure per OB (`PLOT_EACH_ONE = False`)

Left: every CCF of the OB. Right: all rest-framed spectra over the template (black), full
range plus two zooms. This is the figure to inspect routinely: misaligned lines or double
peaks jump out immediately.

## 9. Using the engine on your own spectra

`CrossCorrelation.py` has no pipeline dependencies. A complete example is in
`Fast_usage/Fast_RV.py` ([11 — Manual mode](11_manual_mode.md)):

```python
from CrossCorrelation import CCF

myCCF = CCF(template_wavelength = tw, template_flux = tf,
            observed_wavelength = ow, observed_flux = of,
            RV_range = (-120, 120))

RV_sample, Corr, RV, ERR = myCCF.find_rv(steps3=500)
myCCF.get_errs_info(PRINT=True)     # full error breakdown on the terminal
myCCF.PLOT()                        # the five-panel diagnostic
```

## 10. Checklist

* [ ] The template covers the observed range **after** shifting by both ends of `RADVEL_RANGE`
      (set `Check_template_relative_to_obs=True` once to see it).
* [ ] The template's parameters (Teff, log g, [M/H], resolution) are representative of your
      stars.
* [ ] `RADVEL_RANGE` is wide enough to contain the true velocity but narrow enough to avoid
      aliases from repeated line patterns.
* [ ] `CCFS.png` shows single, clean peaks.
* [ ] Outliers in `OUTPUT_RV.txt` have been inspected individually with `PLOT_EACH_ONE=True`.
* [ ] You are aware that `RVERR` is the oversampled-grid floor (§5).
