# 15 — API reference

[← back to index](README.md)

Signatures, arguments and return values of everything importable in the package. Modules are
listed in pipeline order, then the stand-alone tools.

---

## `User_Functions.py`

### `KEYS_FOR_TABLES_READER_FUNCTIONS`
`list[str]` — the contract for reader functions:
`['MJD','RAS','DECS','GLAT','GLON','HTEL','WAVE','FLUX','NOT_SKY','OBJS','SHAPE','MAGS']`.
See [10 — Extending](10_extending_user_functions.md) for the meaning of each key.

### `TABLES_GIRAFFE_MOS_REDUCED(file) -> dict`
Reader for one GIRAFFE MOSSPECTRA FITS file. Returns `FLUX`, `WAVE` (from `CRVAL2`/`CDELT2`/
`NAXIS2`, **nanometres**), `OBJS`, `RAS`, `DECS`, `MAGS`, `SHAPE`, `FIBER`, `NOT_SKY`
(`TYPE not in {'S','s'}`), `MJD`, `GLAT`, `GLON`, `HTEL`.

### `TABLES_UVES_MOS_REDUCED(epoch, ARM='U') -> dict`
Reader for one EsoReflex UVES-fibre exposure **folder**. `ARM` is `'U'` or `'L'`.
Returns `NOT_SKY`, `RAS`, `DECS`, `OBJS`, `MJD`, `GLAT`, `GLON`, `HTEL`, `SHAPE`
(= number of fibre files), `FLUX` `(N_fibres, N_pixels)`, `WAVE` (**Ångström**).
No `MAGS`.
*Note:* `MAIN.py` calls readers with one argument, so bind `ARM` with a `lambda`/`partial`.

### `READ_FITS_CONTINUUM(hdul) -> (wave, flux)`
Example Stage-4 reader. Tries `hdul[1].data['wave']/['flux']`, then
`hdul[1].data['wavelength']/['Flux']`, then a 1-D IRAF layout.

### `TEMPLATE_READER(path) -> (wave, flux)`
Reads AUTOKUR/ATLAS9 ASCII templates: `np.genfromtxt(path, skip_header=47)`, columns 0 and 1.
Set `skip_header` to match your template's header length — see
[10 §4](10_extending_user_functions.md#4-custom-template-reader-for-stage-3-template_reader).

---

## `MAIN.py`

### `merge_dicts(dicts: list[dict]) -> dict`
Stacks the per-exposure dictionaries produced by `read_table` into one dictionary of object
arrays, following the fixed key order
`['MJD','RAS','DECS','NOT_SKY','GLAT','GLON','HTEL','WAVE','FLUX','OBJS','SHAPE','MAGS']`.
Missing keys become `None`. The object-array construction (`np.empty(n, dtype=object)`) is
deliberate: it prevents numpy from trying to broadcast ragged per-exposure arrays.

---

## `List_files_related.py`

### `OBSERVING_BLOCK(path, time, epochs, RA, DEC, not_sky=None, OBS_BLOCK_TIME=12, shapes=None, Fname='Observing_blocks.npy', SAVE_TXT=True) -> None`

| Argument | Meaning |
|---|---|
| `path` | output folder (must end with `/`) |
| `time` | MJD per exposure |
| `epochs` | file/folder names |
| `RA`, `DEC` | list of per-fibre coordinate arrays (deg), one entry per exposure |
| `not_sky` | list of boolean masks; `None` ⇒ all fibres are science |
| `OBS_BLOCK_TIME` | grouping window half-width, **hours** |
| `shapes` | optional, printed in the table |
| `Fname` | output `.npy` name (the `.txt` uses the same stem) |
| `SAVE_TXT` | also write the PrettyTable text version |

Writes `Observing_blocks.npy` / `.txt`, prints one table per block, returns `None`.
Algorithm: [06](06_stage1_observing_blocks.md).

---

## `Merging_and_substract_sky.py`

### `MERGING_MULTIPLE_SPECTRA(MJD, path, Wave, Flux, OFiles, RAS, DECS, GLAT, GLON, HTEL, function_used_to_merge=np.median, NOT_SKY=None, OBJS=None, MAGS=None, SAVE=True, SHOW=True, SHOW_SKIES=True, window=[None], AXIS_TO_MERGE=0, wave_to_angstrom=1, coord_sys='fk5')`

Processes **one Observing Block**: sky subtraction, flux scaling, barycentric correction,
resampling, co-addition, writing. All list-like arguments hold one entry per exposure of the
block. Returns `None`. Details: [07](07_stage2_sky_and_merging.md).

### `MISMATCH_SECOND_DIMENSION(fluxes, waves, tolerance=0.02) -> list`
Interpolates every exposure onto `waves[0]` when the pixel grids differ. Raises `ValueError`
if the edges disagree by more than `tolerance` (fraction of the median spectral width).

### `new_fits(Fname, OBJ, RA, DEC, MJD, Ofile, MJDFile, Flux, Wave, GLAT, GLON, HTEL, CRVAL=None, CDELT=None, MAG=None, path='./', count=1, wave_to_angstrom=1, RVCORR=[0])`
Writes one merged-spectrum FITS. Supplying `CRVAL`/`CDELT` also stores the flux as a 1-D
primary image (IRAF-style). Layout: [05 §3](05_file_formats.md#3-merged-sky-subtracted-spectrum).

### `save_sky(skies_obs, Fname, OBJ, RA, DEC, MJD, Ofile, flux, wave, GLAT, GLON, HTEL, path='./', function_used_to_merge=np.median, AXIS_TO_MERGE=0, wave_to_angstrom=1)`
Writes one `SKY_<MJD>.fits` (individual fibres as `ImageHDU`s + the master sky as a binary
table). Creates `path` if needed.

### `extend_linear_range(min_val, max_val, num_points, rv_kms=50.0) -> ndarray`
Linear grid over `[min_val, max_val]` extended by ±`rv_kms` while keeping the pixel step
implied by `num_points`.

### `malla(x, PRINT=0) -> (rows, cols)`
Near-square subplot layout for `x` panels. `PRINT=1` draws an ASCII sketch.

---

## `CrossCorrelation.py`

### `class CCF`

```python
CCF(template_wavelength, template_flux, observed_wavelength, observed_flux,
    RV_range=None, IGNORE=False, PRINT=False,
    Cosmic_rays=True, Check_template_relative_to_obs=False)
```

| Argument | Meaning |
|---|---|
| `template_*`, `observed_*` | equal-length wavelength/flux pairs, same units |
| `RV_range` | `[min, max]` km s⁻¹. `None` ⇒ derived from the coverage (and forces the iterative search) |
| `IGNORE` | skip all consistency checks (then `RV_range` is mandatory) |
| `PRINT` | verbose checks |
| `Cosmic_rays` | clip flux above the 99th percentile (also flattens emission lines) |
| `Check_template_relative_to_obs` | plot the maximally blue/red-shifted template over the spectrum |

**Attributes set at construction:** `t_wave`, `t_flux`, `obs_wave`, `obs_flux`, `RV_range`,
`ERR_init` (the sampling floor), `Computed`.

#### `find_rv(steps3=300, steps1=100, steps2=100, Tolerance=90, Manual=False) -> (RV_sample, Correlation, RV, ERR)`
The three-pass search plus the four analytic peak fits. Also sets `RV`, `RV_error`,
`RV_median`, `RV_std`, `median_err_fits`, `RV_sample`, `Correlation`, `idx_peak`,
`params_F_gauss`, `params_F_parab`, `params_parabola`, `params_gauss`.
Each `params_*` is `[centre, σ_centre, popt, pcov]`, or four `nan`s if the fit failed.

#### `cross_correlation(RV) -> float`
One CCF evaluation: shift the template by `RV`, interpolate onto the observed grid
(`fill_value=0`), return `Σ obs·template`.

#### `get_errs_info(PRINT=False, ALL_info=False)`
Returns `ERR`, or `(RV, ERR, median_centres, std_centres, median_fit_errors)` with
`ALL_info=True`. `PRINT=True` prints the full breakdown.

#### `PLOT(RV=None, window=None, window1=None, window2=None, window3=None, zoom=True, PAUSE=False, time=1, Save=False, Save_name='Output_CCF.png') -> None`
The five-panel diagnostic figure.

#### `plot_fit_peak(AX0=None, RV=None)`
Draws only the zoomed peak with the four fits; pass your own axes to embed it.

#### `find_peaks(Y, npeak=1, peak='max', tolerance=60, SHOW=False) -> (peak_idx, width_peak)`
Standalone peak finder: extremum plus the `[left, right]` indices where the curve falls below
`Y_min + (Y_max−Y_min)·tolerance/100`. `SHOW=True` draws an explanatory figure.

#### Static-like helpers
`gaussian(x, x0, sigma, a=1)`, `parabola(x, x_center, a, c=1)`, `scale_to_1(x)`,
`INSET_ZOOM(AX1, AX2, color='fuchsia', lw=5, alpha=0.5, c1=True, c2=False, c3=True, c4=False)`.

---

## `Getting_RV.py` (script, `exec`'d by `MAIN.py`)

Consumes from the namespace: `path_folder`, `spectrums`, `EXTRA_FOLDER_RV`, `WAVE_TEMP`,
`FLUX_TEMP`, `RADVEL_RANGE`, `MIN_WAVE_CCF`, `MAX_WAVE_CCF`, `OBJ_Template`, `PLOT_EACH_ONE`,
`SAVE_SUMMARY_PLOT`, `WINDOW_1`, `WINDOW_2`, `complete_path_template`.

| Function | Purpose |
|---|---|
| `corr_rv(wave, rv)` | `wave/(1+rv/c)` — shift to rest frame |
| `get_w_f(path)` | → `(wave, flux, OBJ, MAG, RVC, RA, DEC)`; trims the edges and rescales the flux |
| `MYFigure()` | → `(ax0, ax1, ax2, ax3)`, the summary-figure layout |
| `UPDATE_header(path, RV, RVERR, RVTRUE)` | writes the three keywords into an existing file |
| `new_fits(path, name, RV)` | writes the rest-frame copy into `path+'RESTFRAME/'` |
| `get_rv(path)` | ties it together; returns `(RV_sample, Corr, RV, ERR, wave, flux, OBJ, MAG, RVC, RA, DEC)` |

Script-level switches: `method` (`'Synthetic'`/`'Template'`), `path_temp_star`, `SAVE`,
`SAVE_NAME`.

---

## `Continuum_fitting.py`

### `class cont_fit`

```python
cont_fit(spectra, path='./', RUN=True, Continue=False,
         wave_range=(-1, np.inf), sigma_range=(0.1,5), iter_max=50, smooth_range=(0,5),
         iter_default=0, sigma_L_default=1.25, sigma_U_default=2.5, smooth_default=0.03,
         INP_OUT_FILE='Stored_params.npy', READ_FUNCTION=None, wave_units='Angstrom')
```

FITS input only. `READ_FUNCTION(hdul) -> (wave, flux)`; `None` uses `default_read`.
Parameters are saved on **every** exit path (including exceptions).

| Method | Purpose |
|---|---|
| `default_read(hdul)` | built-in reader (binary table, then 1-D IRAF) |
| `get_spectra(name)` | read one spectrum, apply `wave_range`, load or compute + cache the guess |
| `CONT_GUESS()` | the shift-and-percentile initial guess; sets `init_mask`, `cont_guess` |
| `CONTINUUM(wavelength, flux, sigma_clip=(1.25,2.5), max_iter=15, smoothing_factor=1)` | the continuum for a given parameter set; implements the `−1 / 0 / 1 / >1` semantics |
| `Figure_setup(N=0)` | build the GUI |
| `update`, `change_spectrum`, `prev_spectrum`, `next_spectrum`, `increment_slider` | widget callbacks |
| `save(Output=…)`, `read(Input=…)` | persist / restore `stored_params` |
| `RUN(n=0)` | open the GUI at spectrum `n` |
| `UPDATE_fits(path, name, key)` | write one normalised spectrum |
| `store_cont_sub(key='Continuum_sub/', path=None, Force=False)` | write them all |

Guess parameters (attributes, edit before `RUN`): `seed=170499`, `shift=180` (km s⁻¹),
`Nsamples=2500`, `percentile_of_cont=86`, `percentile_range=(80,95)`.

---

## `automated_cont_fit.py`

### `class automated_continuum`

```python
automated_continuum(wave, flux, name, path='', subfolder='CONT_SUB/',
                    wave_range=(-1, np.inf),
                    percentile_cont=(84,92), percentile_range=(75,93), percentile_edges=(65,95),
                    iter_max=10, RUN=True, wave_units='Angstrom',
                    save_plots=False, primary_header=None)
```

With `RUN=True` the constructor runs the convergence, creates `path+subfolder` and writes the
FITS.

| Method | Purpose |
|---|---|
| `CONT_GUESS(iter)` | one iteration of the shift-and-percentile estimate → `(mask, cont)` |
| `search_convergence(plot=False, save=False, key='', percentage_of_convergence=0.005)` | iterate until the accepted-point count stabilises |
| `final_cont(smoothing_factor=0.04)` | cubic-spline smoothing of the converged continuum |
| `RUN_app()` | convergence + folder creation + write |
| `UPDATE_fits(key='')` | write the 6-column product |

Guess parameters: `seed=170499`, `shift=200` km s⁻¹, `Nsamples=2500`.

### `default_read(full_path, table=1, return_header=True) -> (wave, flux[, header])`
Path-based reader for the automatic mode.

### `multiple_spectra_cont(spectra, path='./', wave_range=(-1,np.inf), percentile_cont=(84,93), percentile_range=(75,95), percentile_edges=(65,98), iter_max=10, wave_units='Angstrom', READ_FUNCTION=None, subfolder='CONT_SUB/')`
Loops `automated_continuum` over a list of filenames, prefixing outputs with `CS-`.

---

## `single_isochrone_fitting.py`

### `class ISO_FIT`

```python
ISO_FIT(MAG, COLOR, iso_path=None, OBJ=None, RA=None, DEC=None,
        photometry='Gaia', Header=13, Stages=[1,2,3,4,7,8])
```

Opens the GUI immediately (`run_app`). `iso_path=None` → a file dialog.

| Method | Purpose |
|---|---|
| `photometric_system(system)` | set `self.system` (Mag1/Mag2/Mag3) and `self.absorption_coeff`; derives `abs_mag`, `abs_col` |
| `isochrone_read(path=None, update=True)` | read the PARSEC file, filter `Stages`, list `MH`/`AGE`, set initial `AV`/`UM` |
| `change_iso(update=True)` | select the current (MH, logAge) isochrone |
| `fitted_isochrone()` | apply `UM` and `AV` to the isochrone |
| `Figure_setup(id=170499)` | build the GUI |
| `update`, `update_mh_age`, `renew_points_iso`, `increment_slider`, `increment_slider_mh_age` | widget callbacks |
| `increment_points(n_points=10000, max_jump_allowed_sigmas=4)` | densify the isochrone, splitting at large jumps |
| `match_nearest(Output_file='Isochrone_values.txt')` | nearest-point assignment + write the table + Kiel diagram |
| `waiting_figure(id=1704)`, `kiel_diagram(indexes, id=1704)` | progress notice and result plot |
| `ask_path()` | file dialog |

---

## `multiple_isochrone_fitting.py`

### Module-level functions

#### `log_likelihood_cmd(delta_mag, delta_col, inv_cov_matrix, N=2) -> (weights, D_min, N_w)`
Normalised weights `exp(logL − max logL)/Σ`, the minimum Mahalanobis distance
`(max logL/(−N/2))^½`, and the number of contributing points.

#### `weighted_mean_and_std(variable, weights) -> (mean, -std_lower, std_upper)`
Weighted mean and **asymmetric** weighted dispersions (points above / below the mean treated
separately). The lower error is returned **negative**.

#### `line_of_sight_of_point(delta_y, delta_x, nbins=360) -> indices`
One point per angular bin: the closest in each 1° cone.

#### `init_worker(...)`, `process_single_point(i)`, `process_chunk(chunk_indices)`
Multiprocessing plumbing. `init_worker` publishes the shared arrays as globals in each worker;
`process_single_point` returns `(i, result_dict)` with every isochrone key, its `_err_l` /
`_err_u`, plus `good`, `d_min`, `NP_mod`.

### `class ISO_FIT`

```python
ISO_FIT(MAG, COLOR,
        err_mag1=0.05, err_mag2=0.05, err_mag3=0.05, err_dist=0.03, err_redd=0.05,
        iso_path=None, degree=2.0, N_bins=360,
        OBJ=None, RA=None, DEC=None,
        photometry='Gaia', Header=13, Stages=[1,2,3,4,7,8])
```

`OBJ`, `RA`, `DEC` are effectively required ([14 §A.12](14_troubleshooting.md#a12-typeerror-nonetype-object-is-not-subscriptable-multiple-isochrone-fitting)).

| Method | Purpose |
|---|---|
| `photometric_system(system)` | as above |
| `inverse_covariance_matrix_cmd(err_mag1, err_mag2, err_mag3, err_dist, err_redd)` | build and invert Σ (one matrix, or one per star if arrays are given) |
| `isochrone_read(path=None)` | load the whole grid; set `all_ages`, `all_mh`, `Iso_mask` |
| `Figure_setup(id=170499)` | GUI with two `RangeSlider`s (MH, logAge) and two `Slider`s (A_V, distance modulus) |
| `update`, `update_mh_age`, `fitted_isochrone`, `increment_slider`, `string_age` | widget callbacks and helpers |
| `match_nearest(Output_file='Isochrone_values', min_mag_dist=0.075, min_color_dist=0.035, max_mag_dist=1, max_color_dist=0.3)` | the full Bayesian inference over all cores |
| `save_fits(params, Fname, path='./')` | write the FITS binary table |
| `save_plain_txt(params, Fname, path='./')` | write the PrettyTable text version |
| `kiel_diagram(teff, logg, id=1704)`, `waiting_figure(id=1704)` | plots |

---

## `Fast_usage/` modules

### `barycentric_vels.py`
* `get_target_coords(target_name) -> SkyCoord` — SIMBAD resolution (needs `astroquery`).
* `barycentric_velocity(target_coord, times, location) -> Quantity`
* `plot_barycentric_velocity(target_name, mjd_range=None, n_points=100, mjd_array=None, observatory='Paranal', savefile=None, return_values=False)`
* `compute_components(target_coord, times, location) -> (v_orb, v_rot, v_total)`
* `components_plot(target_name='Alnitak', observatory='Paranal', mjd_start=60000.0, mjd_end=60001.0, n_points=400) -> dict`

### `stellar_oscillation_params.py`
* `stellar_oscillation_params(log_L, T_eff, log_g, gamma=5/3, mu=1.3, …) -> dict`
  with `R_Rsun`, `M_Msun`, `rho_mean_gcm3`, `rho_mean_solar`, `nu_ac_uHz`, `nu_ac_phys_uHz`,
  `P_ac_days`, `nu_max_uHz`, `P_max_days`, `delta_nu_uHz`, `log_g_sun_used`.
* `stellar_oscillation_table(log_L_arr, T_eff_arr, log_g_arr, **kw) -> dict[str, ndarray]`

### `convective_turnover.py`
* `photospheric_density(g_cgs, T_eff, mu, Z_rel) -> float`
* `convective_turnover_local(log_L, T_eff, log_g, alpha_MLT=1.6, gamma=5/3, mu=1.3, Z_rel=0.2) -> dict`
  with `H_p_cm`, `rho_phot`, `kappa_phot`, `v_conv_cms`, `v_conv_kms`, `tau_local_s`,
  `tau_local_d`, `R_Rsun`, `M_Msun`.
* `convective_turnover_global(log_L, T_eff, log_g, alpha_MLT=1.6, tau_sun_global=12.0) -> dict`
* `convective_turnover_full(log_L, T_eff, log_g, **kw) -> dict` — local + global +
  `P_max_days` + `ratio_Pmax_tau`.
* `convective_turnover_table(log_L_arr, T_eff_arr, log_g_arr, **kw) -> dict[str, ndarray]`

### `Fast_RV.py`, `Fast_continuum.py`
Example scripts, not libraries. See [11 — Manual mode](11_manual_mode.md).

---

## `plotingRVS.py`

Configured by its own `CONFIGURATION` block and run as a script (`main()`).
See [13 §1](13_auxiliary_tools.md#1-plotingrvspy--multi-epoch-radial-velocity-review).

| Function | Purpose |
|---|---|
| `gaussian(x, x0, sigma, a=1)` | model for the RV histogram fit |
| `Folding(t, width, ref=None)` | phase-fold a time array → phases in `[0,1)` |
| `reorderLegend(ax=None, order=None)` | force a legend order |
| `variable_bandwidth_kde(X, bandwidths, query_points, kernel='gaussian'|'tophat')` | KDE with a per-point bandwidth |
| `robust_scatter(x)` | MAD scaled to a Gaussian σ |
| `discover_obs(path, requested)` | OB list, from `LIST` or `Observing_blocks.npy` |
| `read_rv_table(filename)` | parse `OUTPUT_RV.txt` → `dict` of arrays, columns taken from its header |
| `epoch_time(path, obs, subfolder)` | MJD of an epoch from any spectrum header, or `nan` |
| `build_catalogue()` | → `(names, coords, mjd, obs_ids, rv, rve)` with rejections applied |
| `get_periodogram_backend()` | → `('gls', Gls)` \| `('astropy', LombScargle)` \| `(None, None)` |
| `best_period(t, y, dy, backend, cls)` | → `(period, T0, amplitude)` or `nan`s |
| `figure_variables(...)`, `figure_distribution(...)` | the two figures |
| `main()` | the whole analysis |

## `Review_Multiple_Epochs.py`

Interactive viewer, configured by its own `CONFIGURATION` block.
See [13 §2](13_auxiliary_tools.md#2-review_multiple_epochspy--interactive-multi-epoch-viewer).

| Function | Purpose |
|---|---|
| `READ_SPECTRUM(hdul)` | `(wave, flux)` from an open `HDUList` — **the hook to edit for your own layout** |
| `READ_TEMPLATE(path)` | `(wave, flux)` for the reference template (ASCII or FITS) |
| `boxcar(y, n)` | running mean, display only |
| `gauss_profile(x, center, sigma, depth, offset)` | absorption-line model |
| `fit_line_profile(wave, flux, center, halfwidth)` | → dict with `center`, `center_err`, `sigma`, `fwhm`, `depth`, `offset`, `ew`, `rv_shift`, `x`, `model`; `None` on failure |
| `line_bisector(wave, flux, center, halfwidth, n_levels=25)` | → `(bisector_wavelengths, levels)` or `(None, None)` |
| `discover_obs(path, requested)` | OB list, honouring `DISCARD_OBS['…'] == '*'` |
| `star_name(filename)`, `coords_from_name(filename)` | parse `<OBJ>=<RA>=<DEC>=<DATE>.fits` |
| `build_index()` | → `(obs_ids, stars)`, `stars` = list of `{'name','ra','dec','files'}` |
| `main()` | build the index, open the viewer |

### `class Viewer`

```python
Viewer(obs_ids, stars)
```

| Method | Purpose |
|---|---|
| `load(i)` | read + resample one star's epochs onto a common grid; cached |
| `reference_of(grid, fluxes)` | → `(reference, label)` honouring `REFERENCE`, with fallbacks |
| `update(i)` | redraw everything for star `i` |
| `step(delta)` | move the slider |
| `on_pick(event)`, `on_key(event)` | click-to-select and keyboard shortcuts |

## `cheating.py`
No importable API — four `if 1:`/`if 0:` recipe blocks. See [13 §3](13_auxiliary_tools.md#3-cheatingpy--scripted-bookkeeping-recipes).

## `cheating.py`
No importable API — four `if 1:`/`if 0:` recipe blocks. See [13](13_auxiliary_tools.md).
