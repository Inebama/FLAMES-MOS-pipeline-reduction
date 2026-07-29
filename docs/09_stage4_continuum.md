# 09 — Stage 4: continuum estimation and normalisation

[← back to index](README.md)

**Modules:** `Continuum_fitting.py` (class `cont_fit`, interactive) and
`automated_cont_fit.py` (class `automated_continuum` + helper `multiple_spectra_cont`,
unsupervised)
**Driver:** the `if CONTINUUM …` branches of `MAIN.py`, switched by `AUTOMATED_CONT`
**Products:** `Continuum_sub/*.fits` (interactive) or `CONT_SUB/CS-*.fits` (automatic),
`Stored_params.npy`, and a cached `CONTINUUM GUESS` extension inside the input spectra

---

## 1. The idea behind the initial guess

Classical continuum placement fits a low-order function to points selected by sigma-clipping.
That works when the continuum is flat, but fails when the *shape* of the continuum varies —
instrumental response, order blaze, molecular bands — because then a global flux distribution
no longer isolates the continuum.

The approach used here builds a **local** flux distribution at every pixel:

1. Draw `Nsamples = 2500` velocity shifts from a Gaussian, `v ~ N(0, shift)` with
   `shift = 180 km/s` (interactive) or `200 km/s` (automatic).
2. For each shift, Doppler-shift the wavelength axis, `λ' = λ/(1+v/c)`, and re-interpolate the
   flux back onto the **original** pixel wavelengths. Out-of-range pixels are filled with a
   representative edge value (the median of the corresponding end of the spectrum).
3. At every pixel you now have 2500 flux values sampled from a small wavelength neighbourhood
   — a *local* flux distribution whose upper envelope traces the continuum.
4. Take a **high percentile** of that distribution (86th interactive, 84→92 ramp automatic) as
   the continuum estimate at that pixel.

Because absorption lines are narrow and sparse compared with the shift width, they occupy the
low tail of each local distribution, while the continuum dominates the high end. The method is
insensitive to the global continuum shape by construction, since each pixel is compared only
with its own neighbourhood.

Then a mask of *continuum candidate* pixels is defined from the normalised flux
`NORM = flux / cont_guess`:

```python
b = percentile(NORM, PerL)      # 80 (interactive) / 75 (automatic)
a = percentile(NORM, PerU)      # 95 (interactive) / 95 (automatic)
init_mask = (b < NORM) & (NORM <= a)
```

i.e. pixels lying in a band just below and around unity — high enough to exclude line cores,
low enough to exclude noise spikes and cosmic rays.

**Known limitation of the shifting technique.** Where the spectral flux changes rapidly with
wavelength (steep blaze, strong broad features, the spectrum edges), the shifted copies no
longer overlap in continuum level and the percentile is biased. Residuals grow there. Two
mitigations are built in: the automatic mode uses a *relaxed* percentile band at the edges
(`percentile_edges`), and iterating the procedure on the successively normalised spectrum
makes the copies overlap better with each pass.

Everything is driven by `np.random.seed(170499)`, so the guess is **exactly reproducible**.

## 2. Interactive mode — `Continuum_fitting.cont_fit`

### 2.1 Launching

`MAIN.py` (with `CONTINUUM=True`, `AUTOMATED_CONT=False`) asks which OB to work on, then:

```python
CONT = cont_fit(spectrums,                       # list of filenames (from Oficial_list.list)
                path_folder + EXTRA_FOLDER_CONT, # where they live
                wave_range    = WAVE_RANGE,
                Continue      = True,            # resume from Stored_params.npy
                READ_FUNCTION = read_fits)       # None ⇒ built-in reader
CONT.store_cont_sub(Force=False)
```

Full constructor:

```python
cont_fit(spectra, path='./', RUN=True, Continue=False,
         wave_range=(-1, np.inf),
         sigma_range=(0.1,5), iter_max=50, smooth_range=(0,5),
         iter_default=0, sigma_L_default=1.25, sigma_U_default=2.5, smooth_default=0.03,
         INP_OUT_FILE='Stored_params.npy', READ_FUNCTION=None, wave_units='Angstrom')
```

> **FITS only.** The interactive tool reads spectra through `astropy.io.fits`; ASCII spectra
> must be converted first. `READ_FUNCTION` receives an already-open `HDUList` and must return
> `(wave, flux)`. The built-in reader tries `hdul[1].data['wave']/['flux']`, then a 1-D IRAF
> layout (`hdul[0].data` + `CRVAL1`/`CDELT1`).

### 2.2 The window

```
┌───────────────────┬──────────────────────────────────────────────────────┐
│                   │   ── flux (blue)                                     │
│  Current Spectra  │   ── init guess (green, always shown for reference)  │
│  ▲ vertical       │   ── current continuum (red)                         │
│  │ slider         │                                                      │
│  │ (list of       │                                                      │
│  │  spectra)      │                                                      │
│                   ├──────────────────────────────────────────────────────┤
│                   │  Sigma below      [────●────]        [-] [+]         │
│                   │  Sigma above      [──●──────]        [-] [+]         │
│ [Previous][Next]  │  Iterations       [●────────]        [-] [+]         │
│                   │  Smoother factor  [─●───────]        [-] [+]         │
└───────────────────┴──────────────────────────────────────────────────────┘
```

* The **figure title** is the current filename.
* The vertical slider on the left selects the spectrum; *Previous*/*Next* step through the list.
* Each horizontal slider has `+`/`−` buttons for fine steps (0.01 for the sigmas and the
  smoothing factor, 1 for the iterations).
* **Closing the window ends the program.** Parameters are saved on every exit path, including
  exceptions (`finally: self.save(...)`).

### 2.3 The four parameters

| Slider | Range | Meaning |
|---|---|---|
| *Sigma below* | 0.1–5 | Lower clipping threshold: pixels with `residual < −σ_low·std` are rejected (line cores). |
| *Sigma above* | 0.1–5 | Upper clipping threshold: pixels with `residual > +σ_high·std` are rejected (spikes). Asymmetric on purpose — you want to keep the upper envelope. |
| *Iterations* | −1 … 50 | See the table below. |
| *Smoother factor* | 0–5 | The `s` parameter of `scipy.interpolate.UnivariateSpline` (cubic). `s=0` interpolates every point; larger values smooth more. |

### 2.4 Iteration semantics

This is the heart of the tool (`cont_fit.CONTINUUM`):

| `Iterations` | Behaviour | Which sliders matter |
|---|---|---|
| **−1** | Return the statistical initial guess unchanged. | none |
| **0** | `UnivariateSpline(wave, cont_guess, s=smooth)` — a smoothed version of the guess. | *Smoother factor* |
| **1** | `UnivariateSpline(wave[init_mask], flux[init_mask], s=smooth)` — one spline through the points inside the initial acceptance band. | *Smoother factor* |
| **> 1** | Classical asymmetric sigma-clipping, `n` passes, starting from `init_mask`: fit spline → residuals → `mask = (r ≤ σ_high·std) & (r ≥ −σ_low·std)` → refit. | all four |

A practical recipe: start at `Iterations = 0` (fast, well-behaved), raise the *Smoother
factor* until the red curve stops chasing the noise; if line wings are pulling the continuum
down, go to `Iterations = 3–5` and tune *Sigma below* first.

### 2.5 Performance and caching

* Only **one spectrum is held in RAM at a time** — each is re-read from disk when selected.
* The expensive initial guess is computed once per spectrum and **appended to the spectrum
  file itself** as the `CONTINUUM GUESS` extension (wave, `init_mask`, `Cont_guess`, plus all
  the parameters used). Revisiting a star is then instantaneous, and the guess is identical
  across sessions.
* Parameters live in `Stored_params.npy` as `{filename: [σ_low, σ_high, iterations, smooth]}`.
  With `Continue=True` the GUI reopens at the first spectrum whose entry is `None`; if the file
  does not exist you get a clear message and a fresh start. If the stored keys do not match the
  current spectrum list, the class refuses to continue (set `Continue=False` to restart).

### 2.6 Writing the output

```python
CONT.store_cont_sub(key='Continuum_sub/', path=None, Force=False)
```

Refuses to write while any spectrum still has `None` parameters, unless `Force=True`. For each
spectrum it re-reads the file, re-applies `wave_range`, recomputes the continuum from the
stored parameters and writes the 6-column product described in
[05 §11](05_file_formats.md#11-normalised-continuum-subtracted-spectrum), with the exact
parameters recorded in the extension header.

## 3. Automatic mode — `automated_cont_fit.py`

```python
multiple_spectra_cont(spectra, path='./',
                      wave_range       = (-1, np.inf),
                      percentile_cont  = (84, 93),
                      percentile_range = (75, 95),
                      percentile_edges = (65, 98),
                      iter_max         = 10,
                      wave_units       = 'Angstrom',
                      READ_FUNCTION    = None,
                      subfolder        = 'CONT_SUB/')
```

`MAIN.py` calls it per OB with only `spectra`, `path` and `wave_range`, so the other defaults
apply. `READ_FUNCTION` here takes a **path** (not an `HDUList`) and returns
`(wave, flux, header)`; the default `default_read` handles the pipeline's own products and the
1-D IRAF layout, and passes the original primary header through to the output.

### The convergence loop

```
for i in 0 … iter_max−1:
    cont      ← percentile(shifted copies of current_flux, percentile_cont[i])
    NORM      ← current_flux / cont
    mask      ← percentile band percentile_range   (relaxed to percentile_edges near the edges)
    n_points  ← number of accepted pixels
    current_flux ← NORM                   (from the 2nd iteration on)
    stop when n_points changes by < 0.5 % of the initial count for two consecutive iterations
final continuum ← UnivariateSpline(wave, flux/current_flux, s=0.04)
```

Key points:

* `percentile_cont` as a **tuple** ramps the continuum percentile linearly with iteration
  (84 → 93 as called by `MAIN.py`; 84 → 92 if you instantiate `automated_continuum` directly).
  Starting low and finishing high lets the estimate approach the true envelope from below
  without locking onto noise.
* The spectrum is *successively normalised*: each iteration works on the previous iteration's
  `NORM`, so the copies overlap better and better (this is exactly the mitigation described in
  §1).
* The two edge regions — defined as the pixels within `3·shift` in velocity of each end — get
  the relaxed band `percentile_edges`, except in the last iteration.
* `iter_max` is forced to be at least 3. `percentile_cont[-1]` must not exceed
  `max(percentile_range)`, else a `ValueError` is raised.
* The final continuum is `flux / current_flux` (the cumulative product of all normalisations)
  smoothed with a cubic spline, `s = 0.04`.

Set `save_plots=True` on `automated_continuum` to get three diagnostic PNGs per spectrum: the
continuum guess at each iteration, the accepted continuum points, and a bar chart of the
convergence of the point count.

### Class signature

```python
automated_continuum(wave, flux, name,
                    path='', subfolder='CONT_SUB/',
                    wave_range=(-1, np.inf),
                    percentile_cont=(84,92), percentile_range=(75,93), percentile_edges=(65,95),
                    iter_max=10, RUN=True, wave_units='Angstrom',
                    save_plots=False, primary_header=None)
```

With `RUN=True` (default) the constructor immediately runs the convergence, creates the output
folder and writes the FITS. `RUN=False` lets you build the object, inspect/patch it and call
`search_convergence()` / `final_cont()` / `UPDATE_fits()` yourself.

## 4. Interactive or automatic?

| | Interactive | Automatic |
|---|---|---|
| Effort | ~5–20 s per spectrum of human attention | none |
| Reproducibility | full (parameters stored) | full (deterministic, seeded) |
| Output folder | `Continuum_sub/` | `CONT_SUB/` |
| Filename | unchanged | prefixed `CS-` |
| Robustness near strong/broad lines | you can see and fix it | can bite the line wings |
| Robustness at the spectrum edges | you can see and fix it | relaxed percentiles help, but check |
| Recommended for | final science normalisation | large samples, first pass, homogeneity tests |

The author's recommendation is the interactive mode for anything that feeds a measurement, with
the automatic mode as a fast first look or a consistency check.

> A practical hybrid, implemented in `cheating.py` recipe 1: normalise one OB interactively,
> then **copy `Stored_params.npy` to the other OBs** (matching spectra by target name), so
> every epoch of a given star is normalised with identical parameters. That removes
> epoch-to-epoch normalisation systematics — important when the science is differential.

## 5. `WAVE_RANGE` in detail

The same parameter is interpreted two ways:

| Value | Meaning |
|---|---|
| `(6391.4, 6613)` (tuple) | explicit wavelength window in Å |
| `0.04` (float in 0–1) | trim that fraction of the *total range* at **each** end |

The trimming is applied identically when the guess is computed and when the output is written,
so `wave`, `nflux`, `cont`, `flux`, `Cont_guess` and `init_mask` all have the same length. This
is why the rest-frame files (14 713 rows) yield 13 535-row normalised products with the default
`0.04`.

## 6. Checklist

* [ ] The spectra are in the frame you want (rest frame if you ran Stage 3 first).
* [ ] `WAVE_RANGE` removes the zero-padded pixels left by Stage 2 at the edges.
* [ ] The green *init guess* line follows the upper envelope of the spectrum — if it does not,
      the problem is upstream (bad merge, strong residuals), not in the sliders.
* [ ] The red *continuum* does not dip into broad line wings (Hα, Ca II, molecular bands).
* [ ] `Stored_params.npy` has no `None` left before calling `store_cont_sub`.
* [ ] If you compare epochs, all epochs of a star were normalised with the same parameters.
