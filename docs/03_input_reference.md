# 03 — `INPUT.py` configuration reference

[← back to index](README.md)

`INPUT.py` is the single configuration file of the semi-automatic pipeline. `MAIN.py` starts
with `from INPUT import *`, and every stage module reads its settings from the resulting
namespace.

**Rules**

1. Never change the *name* of a variable or of an imported function — other modules look them
   up by name. Change only their *contents*.
2. All paths must end with `/`. The code concatenates strings; it does not use `os.path.join`.
3. Prefer absolute paths for real work. The `path_package` trick is only there to make the
   bundled test data portable.
4. Run `python INPUT.py` after every edit — it validates the configuration without processing
   anything.

---

## General

| Variable | Type | Default | Meaning |
|---|---|---|---|
| `path_package` | `str` | `__file__.replace('INPUT.py','')` | Absolute path of `Package_Data_Reduction/`, derived automatically. Used to locate `Getting_RV.py` (which `MAIN.py` executes with `exec`) and the bundled test files. |
| `path_input` | `str` | `…/TESTING_FOLDER/DATA_GIRAFFE/` | Folder containing the input exposures (files **or** subfolders — see `FITS_FILE`). Must end with `/`. |
| `path_output` | `str` | `…/TESTING_FOLDER/TEST_OUTPUT/Test_your_attemp/` | Root of the output tree. Created automatically if missing. Must end with `/`. |

### Stage switches

| Variable | Type | Stage | Effect when `True` |
|---|---|---|---|
| `LIST` | `bool` | 1 | Reads all input files and writes `Observing_blocks.npy` (+`.txt`). |
| `MERGE` | `bool` | 2 | Sky-subtracts, barycentric-corrects, resamples and co-adds each OB. Can run without `LIST` (the input metadata is re-read either way, inside `if LIST or MERGE:`), but it needs `Observing_blocks.npy` to already exist from a previous `LIST` run. |
| `RADIAL_VEL` | `bool` | 3 | Runs the cross-correlation and writes rest-frame spectra. |
| `CONTINUUM` | `bool` | 4 | Runs the continuum normalisation (interactive or automatic, see `AUTOMATED_CONT`). |

You can run several stages in one go, but the recommended workflow is one stage at a time,
checking the diagnostics in between.

> `MAIN.py` also contains a commented-out block that overrides these flags from the command
> line (`python MAIN.py -list -merge …`, keywords `-list`, `-merge`, `-continuum`, `-ravel`).

---

## Stage 1 — Observing Blocks

| Variable | Type | Default | Meaning |
|---|---|---|---|
| `FITS_FILE` | `bool` | `True` | `True` → the pipeline searches for **files** matching `*.fit*` inside `path_input`. `False` → it searches for **subfolders** (`*/`), which is the EsoReflex/UVES layout where each exposure is a folder. Sets `key` at the bottom of the file. |
| `Range_in_hours_to_be_consider_same_Observing_Block` | `float` | `24` | Half-width, in hours, of the time window used to group exposures. See the algorithm note below. `12` ≈ "one night"; `< 0.25` essentially isolates individual exposures. |
| `read_table` | callable | `TABLES_GIRAFFE_MOS_REDUCED` | The instrument reader. Receives the full path of one file (or folder) and returns a dictionary of metadata + data. Swap for `TABLES_UVES_MOS_REDUCED`, or your own — see [10 — Extending](10_extending_user_functions.md). |
| `SAVE_TXT` | `bool` | `True` | Also write the human-readable `Observing_blocks.txt`. Recommended. |
| `key` | `str` | `'*.fit*' if FITS_FILE else '*/'` | The glob pattern actually used. Defined in the *fine adjustments* section near the bottom; edit it directly if you need a more specific pattern (e.g. `'*SCI*.fits'`). |

> **How the window is applied.** The grouping is *not* a chained/gap-based clustering. Files
> are sorted by MJD; the first unassigned exposure becomes the seed, and every unassigned
> exposure within `±Range/24` **days of that seed** joins the OB. Then the next unassigned
> exposure becomes the new seed, and so on. Details and the positional consistency test:
> [06 — Stage 1](06_stage1_observing_blocks.md).

---

## Stage 2 — sky subtraction and merging

| Variable | Type | Default | Meaning |
|---|---|---|---|
| `function_used_to_merge` | callable | `np.nanmedian` | Combination function used **three times**: to build the master sky from the sky fibres, to co-add exposures, and to combine the per-exposure median counts. `np.nanmedian` is the most robust choice; `np.nanmean` converges faster to the true value with few, clean exposures. Any function with the signature `f(array, axis=…)` works. |
| `min_wave`, `max_wave` | `float\|None` | `None` | X-axis limits of the QA figure `Substract_SKY.png`, **expressed in the native wavelength units of your data** (nm for GIRAFFE, Å for UVES) — *not* Ångström. `None` → full range. Commented examples in the file: `655 / 657.5` for GIRAFFE HR14b, `5150 / 5300` and `6550 / 6575` for the two UVES 580 arms. Purely cosmetic. |
| `wave_to_angstrom` | `float` | `10` | Multiplicative factor applied to the wavelength array **when writing FITS binary tables**, so that the `wave` column is always in Ångström. `10` for GIRAFFE (pipeline delivers nm), `1` for UVES (already Å). Get this wrong and every downstream stage sees the wrong wavelengths. |
| `AXIS_TO_MERGE` | `int` | `0` | Axis passed to `function_used_to_merge`. `0` means "combine across exposures / across sky fibres", which is correct for the standard `(N_fibres, N_pixels)` layout. The `INPUT.py` self-test inspects your flux shape and refuses to continue if this looks wrong. |

---

## Stage 3 — radial velocities

| Variable | Type | Default | Meaning |
|---|---|---|---|
| `EXTRA_FOLDER_RV` | `str` | `''` | Sub-folder **inside each `OB-n/`** where the spectra to correlate live. `''` → the OB root (i.e. run RV directly on the Stage-2 products). `'Continuum_sub/'` → run RV after the interactive continuum stage. Must end with `/` when non-empty. |
| `complete_path_template` | `str` | `…/GIRAFFE_TEMPLATE_HR14b_R36900.dat` | Full path to the template spectrum. Four templates ship with the package (GIRAFFE HR13 R27000, GIRAFFE HR14b R36900, UVES 580 L-arm R47000, UVES 580 U-arm R47000). |
| `reads_template` | callable | `TEMPLATE_READER` | Function that turns that path into `(wave, flux)`. The bundled one reads AUTOKUR/ATLAS9 `.dat` output. ⚠️ Note that `MAIN.py` calls `TEMPLATE_READER` directly rather than `reads_template`; if you write your own reader, either name it `TEMPLATE_READER` in `User_Functions.py` or set `WAVE_TEMP`/`FLUX_TEMP` explicitly (below). |
| `RADVEL_RANGE` | `tuple[float,float]` | `(-250, 250)` | Radial-velocity search range in km s⁻¹, `(min, max)`, strictly ordered. Narrow it if you know the systemic velocity; it makes the search both faster and safer against spurious peaks. |
| `MIN_WAVE_CCF`, `MAX_WAVE_CCF` | `float\|None` | `None` | Wavelength limits (**Å**) applied to the observed spectrum before correlating, to remove low-efficiency CCD edges and the zero-padded pixels left by Stage 2. Affects **only** the CCF — the rest-frame file is written from the full, untrimmed spectrum. **Set both or neither:** if either is `None` the code falls back to an automatic trim of 1/20 of the range at each end, and a lone value is not used. Ready-made presets for GIRAFFE HR14b and the two UVES 580 arms are commented in the file. |
| `OBJ_Template` | `str` | `'Synth_Spectra'` | Label written in the last column of `OUTPUT_RV.txt` to record what the RVs are relative to. |
| `WAVE_TEMP`, `FLUX_TEMP` | `array\|None` | `None` | Optional: supply the template arrays directly instead of a path. If both are given and have the same length, the reader is skipped. |
| `PLOT_EACH_ONE` | `bool` | `False` | `True` → the full four-panel CCF diagnostic (`CCF.PLOT`) is displayed for **every** spectrum for 1 s. Use for debugging a handful of stars. `False` → a single summary figure per OB instead. |
| `SAVE_SUMMARY_PLOT` | `bool` | `True` | With `PLOT_EACH_ONE = False`: save the summary as `CCFS.png` instead of showing it interactively. |
| `WINDOW_1`, `WINDOW_2` | `tuple` | commented out | Wavelength zoom windows (Å) for the two lower panels of the summary plot `CCFS.png`. Pick regions with strong, well-known lines so that template/spectrum alignment is easy to judge. If a variable is left undefined, set to `None`, or contains a `None`, a generic window derived from the spectrum length is used instead and a message says so. Presets for HR14b and the two UVES arms are commented in the file. |
| `Across_OB` | `bool` | `True` | `True` → loop over **all** OBs automatically. `False` → ask which single OB to process (Tk dialog, or numbered prompt). |

Additional knobs that are *not* in `INPUT.py` but are worth knowing (they live at the top of
`Getting_RV.py`):

| Variable | Default | Meaning |
|---|---|---|
| `method` | `'Synthetic'` | `'Synthetic'`: correlate every star against the synthetic template. `'Template'`: first measure one real star against the synthetic template, shift it to rest, and use *that star* as the template for all the others (requires setting `path_temp_star`). |
| `SAVE_NAME` | `'OUTPUT_RV.txt'` | Name of the RV catalogue written in each OB folder. |
| `steps3=500` | in the `find_rv` call | Number of samples of the final, finest CCF pass. |

---

## Stage 4 — continuum

| Variable | Type | Default | Meaning |
|---|---|---|---|
| `EXTRA_FOLDER_CONT` | `str` | `'RESTFRAME/'` | Sub-folder inside `OB-n/` holding the spectra to normalise. `'RESTFRAME/'` = normalise after the RV stage (default). `''` = normalise the Stage-2 products directly (continuum-first workflow). |
| `read_fits` | callable\|`None` | `None` | Custom reader for your FITS layout. It receives an **already-open `HDUList`** and must return `(wave, flux)`. `None` → use the built-in reader, which tries `hdul[1].data['wave']/['flux']`, then `hdul[1].data['wavelength']/['Flux']`, then a 1-D IRAF-style `hdul[0].data` + `CRVAL1`/`CDELT1`. `User_Functions.READ_FITS_CONTINUUM` is a ready-made example. |
| `WAVE_RANGE` | `float\|tuple` | `0.04` | Wavelength range used for the fit. A **tuple** `(λmin, λmax)` is an explicit range in Å. A **float between 0 and 1** trims that fraction of the total range at *each* end (`0.04` → 4 % left + 4 % right ≈ 8 % of the pixels removed). Use it to get rid of the messy CCD edges across a heterogeneous set of spectra. |
| `AUTOMATED_CONT` | `bool` | `False` | `False` → interactive GUI (`Continuum_fitting.cont_fit`), one OB at a time, output in `Continuum_sub/`. `True` → unsupervised convergence loop (`automated_cont_fit.multiple_spectra_cont`) over all OBs, output in `CONT_SUB/` with a `CS-` filename prefix. |

Interactive-mode defaults that are set in the `cont_fit` constructor (not in `INPUT.py`) and
that you may want to change in `MAIN.py`:

| Argument | Default | Meaning |
|---|---|---|
| `Continue` | `True` (as called by `MAIN.py`) | Resume from `Stored_params.npy` if present. |
| `sigma_range` | `(0.1, 5)` | Slider limits for the two sigma-clipping thresholds. |
| `iter_max` | `50` | Upper limit of the *Iterations* slider. |
| `smooth_range` | `(0, 5)` | Slider limits for the spline smoothing factor. |
| `iter_default`, `sigma_L_default`, `sigma_U_default`, `smooth_default` | `0`, `1.25`, `2.5`, `0.03` | Starting values for a spectrum with no stored parameters. |
| `INP_OUT_FILE` | `'Stored_params.npy'` | Name of the parameter store, written in the same folder as the spectra. |

Automatic-mode parameters. `MAIN.py` calls `multiple_spectra_cont` with only `spectra`, `path`
and `wave_range`, so the **function** defaults below are what the pipeline actually uses; the
`automated_continuum` class defaults differ slightly and are given in brackets.

| Argument | Default (function / *class*) | Meaning |
|---|---|---|
| `percentile_cont` | `(84, 93)` / *(84, 92)* | Percentile of the shifted-flux distribution taken as the continuum. A tuple ramps linearly from the first to the second value across iterations (recommended); an int keeps it fixed. |
| `percentile_range` | `(75, 95)` / *(75, 93)* | Percentile band of the normalised flux used to select "continuum points" each iteration. |
| `percentile_edges` | `(65, 98)` / *(65, 95)* | Relaxed band applied near the two spectrum edges, which converge more slowly because the velocity shifts have less information there. |
| `iter_max` | `10` | Maximum number of iterations (forced to ≥ 3). Convergence stops early when the number of accepted continuum points changes by < 0.5 % over two consecutive iterations. |
| `subfolder` | `'CONT_SUB/'` | Output sub-folder. |

---

## The `__main__` self-test block

Everything below the `# HERE AFTER NOT TOUCH` banner runs **only** when you execute
`python INPUT.py`. It imports every dependency (with the author's tested versions listed in
comments), then performs the checks described in
[01 — Installation §5](01_installation.md#5-the-self-test-python-inputpy).

The final line

```python
THIS_SHOULD_NOT_BE_DEFINED_ON_RUNNING_MAIN_FILE = 170499
```

is a marker: if that name exists in your namespace while `MAIN.py` runs, then `INPUT.py` was
executed as a script rather than imported, which would mean something is wrong with your entry
point. (`170499` is also the random seed used throughout the package, so results are
reproducible.)

---

## Worked configurations

### GIRAFFE HR14b (the bundled test data)

```python
path_input  = '/data/47Tuc/GIRAFFE_HR14b/'
path_output = '/data/47Tuc/OUT_HR14b/'
LIST, MERGE, RADIAL_VEL, CONTINUUM = True, True, True, True

FITS_FILE  = True
Range_in_hours_to_be_consider_same_Observing_Block = 24
read_table = TABLES_GIRAFFE_MOS_REDUCED

function_used_to_merge = np.nanmedian
min_wave, max_wave     = 655, 657.5      # nm — Hα region
wave_to_angstrom       = 10

EXTRA_FOLDER_RV        = ''
complete_path_template = path_package+'TESTING_FOLDER/SYNTH_TEMPLATE/GIRAFFE_TEMPLATE_HR14b_R36900.dat'
RADVEL_RANGE           = (-250, 250)
MIN_WAVE_CCF, MAX_WAVE_CCF = None, None

EXTRA_FOLDER_CONT = 'RESTFRAME/'
WAVE_RANGE        = 0.04
AUTOMATED_CONT    = False
```

### UVES fibre mode, 580 nm setting, lower arm (EsoReflex output)

```python
path_input  = '/data/47Tuc/UVES_reflex/'   # contains one folder per exposure
path_output = '/data/47Tuc/OUT_UVES_L/'

FITS_FILE  = False                          # → key = '*/'
read_table = TABLES_UVES_MOS_REDUCED        # defaults to ARM='U'; see note

min_wave, max_wave = 5150, 5300             # Å
wave_to_angstrom   = 1

complete_path_template = path_package+'TESTING_FOLDER/SYNTH_TEMPLATE/UVES_L_ARM_580_R47000.dat'
MIN_WAVE_CCF, MAX_WAVE_CCF = 4800, None     # keep one as None
```

> **UVES arm selection.** `TABLES_UVES_MOS_REDUCED(epoch, ARM='U')` takes the arm as a second
> argument, but `MAIN.py` calls `read_table(path)` with a single argument, so you always get
> the default `'U'`. To reduce the L arm, wrap it:
> ```python
> read_table = lambda p: TABLES_UVES_MOS_REDUCED(p, ARM='L')
> ```
> (a `lambda` is fine — the name `read_table` is what matters), and point
> `complete_path_template` at the matching arm template.
