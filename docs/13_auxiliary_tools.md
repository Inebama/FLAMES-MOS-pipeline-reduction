# 13 — Auxiliary tools: RV review, multi-epoch viewer, bookkeeping

[← back to index](README.md)

Three tools that sit *after* the pipeline. None is called by `MAIN.py`; each is configured by a
block at the top of its own file.

| Tool | Purpose |
|---|---|
| [`plotingRVS.py`](#1-plotingrvspy--multi-epoch-radial-velocity-review) | combine the `OUTPUT_RV.txt` of several OBs, flag RV variables, periodograms, sample RV distribution |
| [`Review_Multiple_Epochs.py`](#2-review_multiple_epochspy--interactive-multi-epoch-viewer) | interactive star-by-star review of every epoch's spectrum against a reference |
| [`cheating.py`](#3-cheatingpy--scripted-bookkeeping-recipes) | scripted bookkeeping recipes (share continuum parameters, index targets, combine epochs) |

The first two were generalised from the author's original science-case-specific scripts and
carry a header noting that the rewrite was done by AI; the astrophysical intent, statistics and
figure design are the author's.

---

# 1. `plotingRVS.py` — multi-epoch radial-velocity review

Reads the Stage-3 RV catalogues of several Observing Blocks, matches the stars across epochs,
and produces the mean RV per star, a list of RV-variable candidates with periodograms, and the
velocity distribution of the whole sample.

## 1.1 Configuration

Everything lives in the `CONFIGURATION` block; nothing else needs editing.

| Variable | Default | Meaning |
|---|---|---|
| `PATH` | test-output folder | Root of the pipeline output (must end with `/`). |
| `LIST` | `None` | Observing blocks to combine. `None` reads them from `Observing_blocks.npy`; otherwise give a list, e.g. `['OB-0','OB-3','OB-7']`. |
| `SUBFOLDER` | `'/RESTFRAME/Continuum_sub/'` | Where the spectra live. Used **only** to recover each epoch's MJD from a FITS header; `''` looks in the OB folder itself. |
| `RV_FILE` | `'OUTPUT_RV.txt'` | Name of the RV catalogue inside each OB. |
| `RV_COLUMN` | `'RV_corrected'` | Which column to analyse. `RV_corrected` is the barycentric (physical) velocity — normally what you want. |
| `MATCH_MODE` | `'name'` | `'name'` matches stars by their `OBJECT` name (robust for this pipeline's own output); `'coords'` cross-matches on the sky with `MATCH_RADIUS_ARCSEC`. |
| `DISCARD_STARS` | `[]` | Stars to drop entirely, by name. |
| `DISCARD_OBS` | `{}` | Per-epoch rejections: `{'OB-2': ['STAR_A','STAR_B']}`, or `{'OB-5': '*'}` to drop a whole epoch. |
| `MAX_RV_ERR` | `None` | Also drop individual measurements whose formal error exceeds this (km s⁻¹). |
| `VAR_THRESHOLD` | `None` | Combined-σ threshold for flagging a variable. `None` → automatic (see below). |
| `N_SIGMA_AUTO` | `5.0` | Multiplier used by the automatic threshold. |
| `MIN_EPOCHS` | `2` | Minimum surviving epochs for a star to be analysed. |
| `MAX_CANDIDATES_PLOTTED` | `12` | Keeps the time-series figure readable. |
| `GLS_PATH` | `''` | Folder containing `gls.py`, if it is not installed. |
| `PERIOD_RANGE` | `(2.0, 1000.0)` | Period search window in days. Automatically capped at 3× the time span. |
| `OVERSAMPLING` | `2000` | Periodogram `ofac`. Lower it if the search feels slow. |
| `SAVE_SUMMARY` | `'RV_SUMMARY.txt'` | Written into `PATH`. `None` to skip. |
| `SAVE_FIGURES` | `None` | A folder path saves the figures instead of showing them. |

### Recording bad observations

This is the mechanism that replaced the hard-coded special cases of the original script:

```python
DISCARD_STARS = ['MED_RGB_370']                          # drop the star everywhere
DISCARD_OBS   = {'OB-2': ['MED_RGB_370','MED_RGB_402'],  # these two were bad that night
                 'OB-5': '*'}                            # this whole epoch was unusable
MAX_RV_ERR    = 5.0                                      # and drop anything with sigma > 5 km/s
```

Discarded measurements become `NaN` and are ignored by every statistic downstream, so a star
dropped in one epoch still contributes its other epochs.

## 1.2 What it computes

1. **Per-epoch ingest.** The column names are read from the commented header of
   `OUTPUT_RV.txt`, so the reader keeps working if the column order changes; it falls back to
   the documented order if no header is found.
2. **Matching.** The first surviving epoch defines the star list and their order.
3. **Epoch times.** The MJD of each epoch is taken from the header of *any* spectrum in it. If
   no MJD can be read anywhere, the OB order is used as the time axis and the period search is
   skipped — the rest of the analysis still runs.
4. **Combination.**
   ```python
   rv_mean = nanmean(rv, axis=0)
   rv_err  = sqrt(nanmean(rve, axis=0)**2 + nanstd(rv, axis=0)**2)
   ```
   Measurement error and epoch-to-epoch scatter in quadrature, so a genuinely variable star
   ends up with a large σ — which is exactly what the selection keys on.
5. **Variability selection.** With `VAR_THRESHOLD = None` the threshold is
   `median(σ) + N_SIGMA_AUTO × MAD-based scatter of σ`, i.e. it adapts to your dataset's own
   precision instead of a hard-coded number. Candidates are sorted by decreasing σ.
6. **Period search.** Per candidate, on the mean-subtracted time series.

## 1.3 Periodogram backends

The script tries, in order:

1. **GLS** — the generalised Lomb–Scargle of **Zechmeister & Kürster (2009), A&A 496, 577**
   ([2009A&A...496..577Z](https://ui.adsabs.harvard.edu/abs/2009A%26A...496..577Z)).
   Reference implementation: <https://github.com/mzechmeister/GLS>, Python package directly at
   <https://github.com/mzechmeister/GLS/tree/master/python>.
   Either `pip install gls`, or download it and set `GLS_PATH` to the folder holding `gls.py`.
2. **`astropy.timeseries.LombScargle`** — automatic fallback, with a message on the terminal.
3. **Nothing** — the period analysis is skipped and the phase panels say `no period`.

## 1.4 Figures

**RV variables** — one row per candidate: the RV range with error bars (left), the time series
with the best-fit sinusoid (middle), and the phase-folded curve labelled with the period
(right). Rows share axes and only the bottom row carries tick labels.

**RV distribution** — a histogram of the mean RVs plus:

| Curve | Meaning |
|---|---|
| Gaussian fit to the histogram | the classical answer, printed to the terminal as well |
| variable-bandwidth KDE | **each star's own uncertainty is its kernel bandwidth**, so precise measurements produce narrow kernels — an honest way to show a distribution built from heterogeneous errors |
| Epanechnikov KDE | fixed bandwidth = 5 × mean error, for comparison (needs `scikit-learn`; silently skipped if absent) |

plus vertical lines at the mean and median and a rug plot of the individual measurements.

## 1.5 Output

`RV_SUMMARY.txt` in `PATH`:

```
# RV summary from RV_corrected over 2 epochs: OB-0, OB-1
# variability threshold = 0.4680 km/s
#OBJ	RA	DEC	RV	RV_ERR	N_EPOCHS	VARIABLE
```

## 1.6 Reusable pieces

| Function | Purpose |
|---|---|
| `Folding(t, width, ref=None)` | phase-fold a time array, returning phases in `[0,1)` |
| `variable_bandwidth_kde(X, bandwidths, query_points, kernel=...)` | KDE with a per-point bandwidth |
| `robust_scatter(x)` | MAD scaled to a Gaussian σ |
| `reorderLegend(ax, order)` | force a readable legend order |

## 1.7 Example on the bundled test data

With the shipped defaults (2 epochs, 116 stars):

```
116 stars, 2 epochs (116 stars with >= 2 epochs)
   mean       median     std        median(sigma)
   -17.138    -17.714    11.756     0.411

Automatic variability threshold: 0.468 km/s (median sigma + 5 x robust scatter)
  idx Object                        RA         DEC    sigma  N_ep
   99 MED_RGB_1365            6.017417  -72.115056    0.502     2
  110 MED_RGB_752             5.872042  -72.080056    0.473     2
   64 MED_VAR_1142            6.265333  -72.158806    0.469     2
```

The mean of −17.1 km s⁻¹ is a useful sanity check: it is the systemic velocity of 47 Tuc. With
only two epochs no period can be derived, and the phase panels correctly say so.

## 1.8 Using it with a single epoch

This script is about *multi-epoch* scatter, so with one epoch there is no variability
information to extract, and it says so rather than producing something misleading:

```
116 stars, 1 epochs (0 stars with >= 2 epochs)

No star has 2 or more surviving epochs -- nothing to analyse.
Lower MIN_EPOCHS, add more OBs to LIST, or relax DISCARD_OBS / MAX_RV_ERR.
```

If you set `MIN_EPOCHS = 1` it will run, and two things are worth knowing:

* **The RV distribution figure is still fully meaningful** — histogram, Gaussian fit and both
  kernel densities of the sample's velocities. This is the useful output for a single epoch.
* **The "variable candidates" are not variables.** With one epoch `nanstd` over epochs is 0, so
  the combined σ collapses to the formal CCF error. On the bundled data every star then sits at
  the sampling floor of ≈0.385 km s⁻¹ and the flagged stars are merely those with the largest
  formal error — worth checking as possible bad CCFs, but not variability. The period search
  correctly returns nothing, since the time span is zero.

So for one epoch: use `plotingRVS.py` for the distribution, and use
`Review_Multiple_Epochs.py` with `REFERENCE='template'` to inspect individual spectra.

---

# 2. `Review_Multiple_Epochs.py` — interactive multi-epoch viewer

An eyeball tool: step through your stars and see, for each one, **every epoch overplotted plus
a residual against a reference of your choice**.

```
+--------------------------------------------------+---------------------+
|  spectrum: all epochs + reference (black)        |  parameter panel    |
|                                                  |  e.g. Teff vs logL  |
+--------------------------------------------------+  (or RA vs DEC)     |
|  residual: flux / reference, with sigma bands    |  click to jump      |
+--------------------------------------------------+---------------------+
|  [<] [ star slider ] [>]     (arrow keys work too)                     |
+------------------------------------------------------------------------+
```

## 2.1 Controls

| Action | Effect |
|---|---|
| slider, `[<]` `[>]`, ←/→ (or ↑/↓) | change star |
| click a point in the parameter panel | jump to that star |
| `f` | toggle the line-profile fit |
| `b` | toggle the line bisector |
| `s` | save the current figure to `SAVE_PREFIX<star>.png` |
| `h` | print the help on the terminal |

## 2.2 Configuration

| Variable | Default | Meaning |
|---|---|---|
| `PATH` | test-output folder | Output root (must end with `/`). |
| `LIST` | `None` | Epochs to show; `None` reads `Observing_blocks.npy`. |
| `SUBFOLDER` | `'/RESTFRAME/Continuum_sub/'` | Sub-folder holding the spectra (`''` = the OB folder). |
| `FILE_PREFIX` | `''` | Set to `'CS-'` to view the automatic-continuum products. |
| `W_MIN`, `W_MAX` | `None` | Display range; `None` = no limit. |
| `N_SAMPLE` | `None` | Resample every epoch onto a common grid of this many points; `None` uses the first epoch's own grid. |
| `REFERENCE` | `'mean'` | What to divide by: `'mean'`, `'median'`, an epoch name such as `'OB-0'`, or `'template'`. |
| `TEMPLATE_PATH` | `''` | Used only with `REFERENCE='template'`. |
| `NORMALISE` | `lambda w,f: np.median(f)` | Re-scale each spectrum before comparing; `None` leaves it as read. |
| `MATCH_MODE` | `'name'` | `'name'` (the `OBJECT` part of the filename) or `'coords'`. |
| `DISCARD_STARS`, `DISCARD_OBS` | `[]`, `{}` | Same semantics as in `plotingRVS.py`. |
| `SORT_BY` | `'name'` | `'name'`, `'ra'`, `'param_x'`, `'param_y'`, `'n_epochs'` or `None`. |
| `PARAM_X`, `PARAM_Y` | `None` | Dicts `{star_name: value}` for the parameter panel — see below. |
| `PARAM_X_LABEL`, `PARAM_Y_LABEL` | Teff, log L | Axis labels. |
| `PARAM_X_INVERT`, `PARAM_Y_INVERT` | `True`, `False` | Invert axes (Teff wants `True`). |
| `PARAM_BACKGROUND` | `None` | Optional `(x, y)` arrays drawn in grey behind, e.g. the full cluster. |
| `SHOW_SIGMA_BANDS` | `True` | ±1/3/5σ envelope on the residual panel. |
| `SMOOTH_N` | `0` | Boxcar width in pixels, **display only**. |
| `MARK_LINES` | `[(6562.819, 'Hα')]` | Vertical reference lines `(wavelength, label)`. |
| `HIGHLIGHT` | `[]` | Shaded regions `(w1, w2, label, colour)` — for sky residuals, telluric bands, anything you want to keep an eye on. |
| `FIT_CENTER`, `FIT_WINDOW` | 6562.819, 5.0 Å | Window for the profile fit and the bisector. |
| `FIT_ON_START`, `BISECTOR_ON_START` | `False` | Enable a diagnostic for every star instead of on keypress. |
| `CMAP`, `FIGSIZE`, `SAVE_PREFIX`, `VERBOSE` | | Cosmetics. |

## 2.3 The parameter panel

**Reading your own auxiliary data is your job** — the tool does not care where it comes from.
Build two plain dicts keyed by star name and assign them:

```python
import numpy as np
t = np.genfromtxt('my_params.txt', names=True, dtype=None, encoding=None)
PARAM_X = dict(zip(t['OBJ'], t['Teff']))
PARAM_Y = dict(zip(t['OBJ'], t['logL']))
PARAM_X_LABEL, PARAM_Y_LABEL = r'$T_{\rm eff}$ [K]', r'$\log(L/L_\odot)$'
PARAM_X_INVERT = True
```

Any pair of variables works — a Kiel diagram, a CMD, mass vs radius, EW vs luminosity. Stars
missing from the dicts are simply not plotted (you get a count on the terminal). **If you leave
both as `None` the panel shows RA vs DEC**, which the tool always knows from the filenames, so
it is useful out of the box with no extra input.

The red crosshair marks the current star; clicking any point jumps to it.

## 2.4 Reference and residual

The residual panel shows `flux / reference` for each epoch. The reference is chosen with
`REFERENCE`:

| Value | Behaviour |
|---|---|
| `'mean'` | mean of the epochs shown for that star — always available, the default |
| `'median'` | median of the epochs; more robust with many epochs |
| `'OB-3'` | that epoch becomes the reference — the natural choice for differential work |
| `'template'` | the spectrum in `TEMPLATE_PATH`, resampled and re-normalised onto the star's grid |

Robustness is deliberate: a template path that does not exist, or a `REFERENCE` epoch that a
given star lacks, prints one warning and falls back to the mean.

### Using it with a single epoch

It works — nothing crashes, and the spectrum and parameter panels are fully usable. But be
aware of what the residual means:

| `REFERENCE` | 1 epoch | Result |
|---|---|---|
| `'mean'` / `'median'` / that epoch's own name | reference **is** the spectrum | residual ≡ 1 exactly; the panel says so and the σ bands are skipped |
| `'template'` | reference is the model | residual is informative — this is the configuration to use |

Because a flat line at exactly 1.0 could easily be misread as "the epochs agree perfectly",
the viewer detects the case and prints the reason in the residual panel:

> *single epoch compared with itself: residual is 1 by construction*
> *set REFERENCE = 'template' (or another epoch) for a meaningful residual*

with a matching one-time hint on the terminal. So for single-epoch review, set

```python
REFERENCE     = 'template'
TEMPLATE_PATH = './TESTING_FOLDER/SYNTH_TEMPLATE/GIRAFFE_TEMPLATE_HR14b_R36900.dat'
```

and you get a genuine data-minus-model comparison. (Remember that a mismatch then reflects the
template's parameters as much as the star's — on the bundled test data a generic RGB template
gives residuals spanning ≈0.33–2.3 in the deep lines.)

The σ bands are intentionally skipped with one epoch: they are measured *from* the
epoch-to-epoch scatter, so with a single spectrum there is nothing to estimate.

A single **star** also works (the slider simply has nowhere to go).

### The sigma bands

With more than one epoch, the scatter of the residuals is measured after scaling by
`√reference` and then propagated back as `±kσ/√reference`. The point is that the same photon
noise produces a *larger relative* deviation where the reference flux is low — i.e. inside line
cores — so a flat band would badly misrepresent the significance of a residual in a deep line.
The three shaded levels are 1, 3 and 5σ.

## 2.5 Line diagnostics

**Profile fit (`f`).** Fits `offset − depth·exp(−(λ−λ₀)²/2σ²)` inside `FIT_CENTER ±
FIT_WINDOW` on the reference spectrum and reports, in the figure title, the fitted centre, its
velocity shift, the FWHM and an equivalent width:

```
fit: centre 6562.842 A (+1.0 km/s), FWHM 1.226 A, EW 1.131 A
```

**Bisector (`b`).** The midpoint of the two line wings at each of `BISECTOR_LEVELS` flux
levels. A vertical bisector means a symmetric line; curvature indicates asymmetry, which for a
resolved line traces velocity fields — granulation, expansion, or a blend.

Both are exposed as plain functions, so you can use them outside the GUI:

```python
from Review_Multiple_Epochs import fit_line_profile, line_bisector
fit = fit_line_profile(wave, flux, center=6562.819, halfwidth=5.0)
bx, by = line_bisector(wave, flux, center=6562.819, halfwidth=5.0, n_levels=25)
```

## 2.6 Performance

Spectra are loaded **lazily, per star, and cached** — the same philosophy as the interactive
continuum stage. Opening a folder with thousands of files is immediate; only the stars you
actually visit are read from disk, and revisiting one is instant.

## 2.7 Adapting the reader

If your FITS layout differs, edit the one function `READ_SPECTRUM(hdul)` — it receives an open
`HDUList` and returns `(wave, flux)`. The default tries the `nflux`/`flux`/`Flux` columns of the
last and first extensions, then the older `wavelength`/`Flux` naming, then a 1-D IRAF image with
`CRVAL1`/`CDELT1`, so it already covers every product this pipeline writes.
`READ_TEMPLATE(path)` does the same for the template (ASCII by default, FITS if the extension
says so).

---

# 3. `cheating.py` — scripted bookkeeping recipes

Four independent blocks, each guarded by `if 1:` (run) or `if 0:` (skip). Only the first is
enabled in the distributed file. Every block re-imports what it needs, so they are genuinely
independent — enable one at a time.

> The file writes `.npy` files and new FITS products. Read the paths at the top of the block
> you enable **before running it**: block 1 in particular overwrites `Stored_params.npy` in the
> target OBs, and the code says so (`# BEWARE of overwriting this`).

## Recipe 1 — share continuum parameters across OBs *(enabled by default)*

**Problem.** If every epoch of a star is normalised with different parameters, part of the
epoch-to-epoch difference you measure is normalisation noise.

**Solution.** Normalise **one** OB interactively, then copy those parameters to all the others,
matching spectra by target name:

```python
data_or = np.load(path+'/OB-0/RESTFRAME/Stored_params.npy', allow_pickle=True).item()

names = {k.split('=')[0]: k for k in data_or}     # target name → filename in OB-0

for i in ['OB-1', 'OB-2', …]:
    spectrums = np.genfromtxt(path+i+'/Oficial_list.list', dtype='U100')
    params = {j: data_or[names[j.split('=')[0]]] for j in spectrums}
    np.save(path+i+'/RESTFRAME/Stored_params.npy', params)      # ← overwrites
```

This is where the `<OBJ>=<RA>=<DEC>=<DATE>.fits` naming convention pays off: `split('=')[0]`
recovers the target name, which is the only part stable across OBs (the date differs, and the
coordinates can differ in the last digits).

Afterwards, run the continuum stage on the other OBs with `Continue=True` — every spectrum
already has parameters, so you can go straight to `store_cont_sub()`, or step through the GUI
to spot-check.

A `KeyError` here means a target present in one OB is missing from the reference OB.

## Recipe 2 — build a per-target index across all OBs *(`if 0`)*

Cross-matches a master target catalogue against every `INFO_FILES.txt` with a 1.75″ tolerance
and builds

```python
{ 'MED_RGB_107': {'OB-1': 'MED_RGB_107=…fits', 'OB-2': None, …}, … }
```

saved as `<k>Targets_files.npy`. The `SUPER_LIST` dictionary maps a *fine* OB grouping onto
coarser ones ("these 5 fine OBs form one 48-hour epoch"), which is how the author explored
several time-binning choices without re-reducing; three alternative groupings (≈1 h, 48 h,
all-merged) are kept as commented blocks. This index is the input of recipe 4.

## Recipe 3 — rotational broadening test *(`if 0`)*

Implements a limb-darkened rotational kernel and applies it in **log-wavelength space**
(constant velocity steps), which is the correct way to do it:

```python
rotational_profile(v, vsini, epsilon=0.6)                        # limb-darkened profile
apply_rotational_broadening_loglambda(wave, flux, vsini, ...)    # FFT convolution
```

The block compares a synthetic spectrum broadened by this routine against one broadened
natively by AUTOKUR at the same `v sin i` — a validation exercise. Useful if you need to match
a template to rotating stars.

## Recipe 4 — median-combine all epochs of each star *(`if 0`)*

Using the index from recipe 2, for every star it:

* opens the spectrum of each contributing OB (skipping `OB-0` and missing entries),
* tries in turn the column pairs `wavelength/Flux`, `wave/nflux`, `wave/flux` so it works on
  merged, rest-frame or continuum-subtracted products alike,
* stores each epoch as its own `ImageHDU` (resampled ×3, with `CRVAL1`/`CDELT1` and
  `EXTNAME` = the OB name),
* interpolates all epochs onto the **common overlap range** and appends a
  `SPECTRUM merged` binary table with the median.

The result is one multi-extension FITS per star: all epochs plus their median in a single file.

---

## Adapting these scripts safely

1. For `cheating.py`, copy the block you want into your own file rather than editing the original.
2. Replace every hard-coded path with your own.
3. Check `LIST` / `SUPER_LIST` / `Folders` against your actual OB names — a wrong name gives a
   `FileNotFoundError` or, worse, a silent partial result.
4. Back up `Stored_params.npy` before running recipe 1.
5. Remember that `cheating.py` still encodes dataset-specific decisions. Delete what does not
   apply to you — leaving them in is how you get results that are quietly wrong.
   `plotingRVS.py` and `Review_Multiple_Epochs.py` no longer have that problem: their
   equivalents are the explicit `DISCARD_STARS` / `DISCARD_OBS` settings.
