# 14 — Troubleshooting

[← back to index](README.md)

Grouped by symptom. Section [B](#b-behaviours-worth-knowing) lists behaviours that are *not*
errors but will surprise you if you do not expect them.

---

# A. Error messages

## A.1 `ImportError` on `python INPUT.py`

A dependency is missing. See [01 — Installation](01_installation.md):

```bash
python -m pip install numpy scipy astropy matplotlib prettytable
```

`prettytable` is the one people usually lack (`pip install prettytable`).

## A.2 `IndexError: list index out of range` in the `INPUT.py` self-test

```python
files = glob.glob(path_input+key); files.sort(); OUTPUT = read_table(files[0])
```

No file matched. Check that

* `path_input` **ends with `/`**;
* `FITS_FILE` matches your layout — `True` searches `*.fit*` (files), `False` searches `*/`
  (folders, the EsoReflex/UVES case);
* the files really are where you think (`ls $path_input`).

## A.3 `WARNING:; I think that AXIS_TO_MERGE must be 1 for your data`

The self-test compared the two dimensions of your flux array and disagrees with your setting.
For the standard `(N_fibres, N_pixels)` layout the correct value is `0`. Follow the message.

## A.4 `Exception: The Values of RA DEC does not match between files …`

Raised by Stage 2. Some exposure in the OB does not share the science-fibre configuration of
the first exposure — different targets, or a different number of them (matching tolerance
0.5″). Causes and fixes:

| Cause | Fix |
|---|---|
| Two plate configurations on the same night | reduce `Range_in_hours_to_be_consider_same_Observing_Block` so they fall into different OBs, or split the input folders |
| One exposure with a dead/missing fibre | remove that exposure, or relax the check in `List_files_related.py` / `Merging_and_substract_sky.py` |
| Coordinates written with different precision by different pipeline versions | raise the 0.5″ threshold in both modules |

The message prints the MJDs of the conforming and non-conforming files, so you can identify
the culprit immediately.

## A.5 `ValueError: The ends and/or the begining of the spectra that you wanted to force to merge differ more than the tolerance accepted of 2.0%`

From `MISMATCH_SECOND_DIMENSION`: the exposures have different pixel grids **and** their
wavelength ranges differ by more than 2 % of the spectral width. The message tells you the
exact tolerance that would be needed. Either raise `tolerance` in that function (and distrust
the edges), or split those exposures into different OBs.

## A.6 `Exception: The min or max of your Obs spectra exceeds you template …`

From `CCF.__init__` when `RV_range is None`: the template does not cover the observed range at
all. Provide a wider template, or give an explicit `RV_range` — the class will then trim the
observed spectrum for you and tell you it did:

```
In order to test RV=[-250.00,250.00] I'm gona cut the obs spectra since your template is to short in range of wavelength
```

That message is informational, not an error.

## A.7 `WARNING: Some of the fits didn't converge` (Stage 3)

One or more of the four analytic peak fits failed. Its centre and error become `nan` and are
skipped by the `nanstd`/`nanmedian` in the error budget. Occasional occurrences are harmless;
if it happens for most spectra, your CCF peak is malformed — inspect it with
`PLOT_EACH_ONE = True`.

## A.8 The continuum GUI does not open, or the sliders do nothing

Backend problem. Use `TkAgg` from a real terminal:

```bash
MPLBACKEND=TkAgg python MAIN.py
```

Not from an IDE console, not from Jupyter, not over an X-less SSH session. See
[01 §3](01_installation.md#3-matplotlib-backend).

## A.9 `SystemExit: The Stored params are not complete please fill them up …`

`store_cont_sub(Force=False)` refuses to write while some spectrum still has `None` in
`Stored_params.npy`. Either finish the OB (touch **every** spectrum with at least one slider so
its parameters get stored), or change the call in `MAIN.py` to `Force=True` to write only the
spectra that are done.

## A.10 `The spectras given does not match the spectras loaded …` then exit

`cont_fit(..., Continue=True)` found a `Stored_params.npy` whose keys do not match the current
file list (different number of spectra, or a renamed file). Delete/rename that file, or call
with `Continue=False` to start over.

## A.11 `KeyError: 'MH'` / `KeyError: 'label'` (isochrone tools)

The `Header` argument is wrong, so `np.genfromtxt(..., names=True, skip_header=Header)` did not
land on the column-name line. Use `Header=14` for PARSEC v2 / CMD 3.8, `Header=13` for the
older format — and count the leading `#` lines in your actual file if in doubt.

## A.12 `TypeError: 'NoneType' object is not subscriptable` (multiple isochrone fitting)

`OBJ`, `RA` and `DEC` are sliced in the constructor (`OBJ[mask][:nnn]`), so they cannot be
`None`. Pass arrays — dummy ones if you have no identifiers.

## A.13 `! gls not found -- falling back to astropy...` (`plotingRVS.py`)

Informational, not an error: the periodogram is being computed with
`astropy.timeseries.LombScargle` instead of the generalised Lomb–Scargle of
**Zechmeister & Kürster (2009), A&A 496, 577** (<https://github.com/mzechmeister/GLS>, Python
package at <https://github.com/mzechmeister/GLS/tree/master/python>). To use GLS itself, either
`pip install gls` or download it and set

```python
GLS_PATH = '/where/you/put/gls/'
```

If neither backend is available the period search is skipped and the phase panels say
`no period`; every other part of the analysis still runs.

## A.14 `No star has N or more surviving epochs -- nothing to analyse` (`plotingRVS.py`)

Your rejections removed too much, or you have fewer epochs than `MIN_EPOCHS`. Lower
`MIN_EPOCHS` (2 is the minimum for any scatter to exist), add OBs to `LIST`, or relax
`DISCARD_OBS` / `MAX_RV_ERR`.

## A.15 `No spectra found. Check PATH, SUBFOLDER and FILE_PREFIX` (`Review_Multiple_Epochs.py`)

The viewer globbed `PATH + <OB> + SUBFOLDER + FILE_PREFIX + '*.fits'` and found nothing. The
usual causes: `SUBFOLDER` does not match where you actually wrote the products (`''` looks in
the OB folder itself), or you are pointing at the automatic-continuum output and forgot
`FILE_PREFIX = 'CS-'`.

---

# B. Behaviours worth knowing

Not errors — design decisions and consequences of the data flow.

## B.1 The CCF trim needs **both** limits to take effect

`MIN_WAVE_CCF` / `MAX_WAVE_CCF` are applied only when *both* are numbers. If either is `None`
the code falls back to the automatic trim (1/20 of the wavelength range at each end) and any
single value you did supply is not used. On a merged HR14b spectrum
(6381.93–6627.11 Å, 14 713 px):

| `MIN_WAVE_CCF` | `MAX_WAVE_CCF` | window used |
|---|---|---|
| `None` | `None` | 6394.200 – 6614.840 Å (automatic) |
| `6391.4` | `None` | 6394.200 – 6614.840 Å (automatic) |
| `6391.4` | `6613` | 6391.400 – 6612.990 Å (yours) |

So set both, or neither. See [08 §2](08_stage3_radial_velocities.md#2-preparing-the-observed-spectrum-get_w_f).

## B.2 `RVERR` is the sampling floor of the *combined* grid

The dominant term of the RV error budget is `ERR_init = median(Δλ)/median(λ)·c/2`, evaluated on
the Stage-2 output grid. That grid is intentionally ~3× finer than the native pixels, because
the exposures of an OB are shifted by *different* barycentric velocities before co-addition, so
the combined spectrum genuinely samples the line profile more finely than any single exposure.
For the bundled GIRAFFE data this gives 0.384 km s⁻¹ rather than the ≈1.15 km s⁻¹ implied by a
raw 0.05 Å pixel. This is a deliberate choice, not an artefact —
see [08 §5](08_stage3_radial_velocities.md#5-the-error-budget) for the full reasoning and for
where to change the oversampling factor if your data do not justify it.

## B.3 The `RVTRUE` header comment does not match the arithmetic

The comment reads `'True Radial Velocity (Measured RV+RVBARY)'`, but the code stores
`RVTRUE = RV` (the CCF result, already barycentric because Stage 2 shifted the wavelengths) and
`RV = RVTRUE − RVBARY`. The **values** are as described in
[08 §6](08_stage3_radial_velocities.md#6-sign-conventions); only the comment string is
misleading.

## B.4 Exposures with no sky fibres are merged without sky subtraction

You get a warning on the terminal (`# WARNING: This file:"…" was encountered to not have
sky...`) and a zero sky array. Watch for it: it is easy to miss in a long log.

## B.5 Zero-padded pixels at the spectrum edges

Pixels not covered by any exposure after the barycentric shifts are set to `0` in the merged
product. Always trim the edges downstream (`MIN_WAVE_CCF`/`MAX_WAVE_CCF`, `WAVE_RANGE`).

## B.6 Stages 3 and 4 modify their input files

Stage 3 adds `RV`/`RVERR`/`RVTRUE` to the primary header of the Stage-2 products; Stage 4
appends a `CONTINUUM GUESS` extension to the spectra it opens. Data arrays are untouched, but
if you want pristine copies, snapshot the OB folder first.

## B.7 `MERGE` re-reads every input file, and needs `Observing_blocks.npy`

`MAIN.py` reads the input metadata inside `if LIST or MERGE:`, so `MERGE` alone works — but it
still opens and parses **every** input exposure through `read_table`, which is the slow part
for large programmes. It then loads `Observing_blocks.npy`, so that file must already exist
from a previous `LIST` run (`FileNotFoundError` otherwise). Leaving `LIST = True` as well is
harmless: it simply rewrites the identical OB files.

## B.8 The `CONTINUUM` self-test opens a directory

In `INPUT.py`, the `if CONTINUUM:` check does `fits.open(path_output)` — a folder. It only
executes when `read_fits is not None`, so with the default it never runs. If you supply a
custom `read_fits` and run `python INPUT.py`, expect an `IsADirectoryError`; point it at an
actual file to test your reader, or ignore the check.

## B.9 `Oficial_list.list` is read with `dtype='U100'`

Filenames longer than 100 characters are silently truncated, which then fails to open. The
generated names (`<OBJ>=<RA>=<DEC>=<DATE>.fits`) run to ~80 characters, so this only bites with
very long target names. Keep target names short and `path_output` shallow.

## B.10 `Fast_usage/` example paths point at the author's own folders

`Fast_RV.py` and `Fast_continuum.py` reference `TEST_OUTPUT/Test_data_ian/` and
`GIRAFFE_TEMPLATE_HR14b_R37000.dat`. The distributed equivalents are
`TEST_OUTPUT/Test_example_try_to_reproduce_it/` and `GIRAFFE_TEMPLATE_HR14b_R36900.dat`.
Edit before running.

## B.11 `single_isochrone_fitting.py` exits before fitting when run directly

Its `__main__` block is the author's 47 Tuc reddening analysis and calls `sys.exit()` before
reaching `ISO_FIT(...)`. Import the class instead of running the module.

## B.12 Products written before v1.0 use the old `GEOLON`/`GEOLAT` convention

Earlier versions of `new_fits`/`save_sky` wrote the observatory latitude under `GEOLON` and the
longitude under `GEOLAT` (the comment fields were correct, and the barycentric correction was
always computed with the right assignment, so no measurement was affected). Current versions
write both correctly. The reference outputs bundled in `TESTING_FOLDER/TEST_OUTPUT/` were
generated before the change and still carry the old convention — regenerate them if you rely on
those keywords.

---

# C. Sanity checks to run on your own data

| Stage | Check |
|---|---|
| 1 | Number of OBs matches the observing log; the exposure count over all OBs equals the number of input files; `Shapes_of_Data` is constant within a block. |
| 2 | `SKY.png`: sky fibres consistent within each exposure. `Substract_SKY.png`: no emission residuals at sky-line wavelengths. Spot-check one output file: `MCOUNT × flux` should look like raw counts. Verify the `wave` column is in Å (i.e. `wave_to_angstrom` was right). |
| 3 | `CCFS.png`: single clean peaks; lines aligned with the template. The RV histogram of a cluster should be a narrow single peak — a bimodal one usually means a template/velocity-range problem. Compare `RVBARY` against `Fast_usage/barycentric_vels.py` for one epoch. |
| 4 | Compare interactive vs automatic normalisation for a few stars. Check that `nflux ≈ 1` in line-free regions and that broad line wings are not clipped. |
| all | Re-run one OB from scratch in a fresh output folder and confirm the products are byte-comparable in the parts that should be deterministic (the seeds are fixed). |

---

# D. Getting help

Errors, unexpected behaviour, or a data format you cannot make fit: contact the author,
Ian Baeza — <i.baezamardones@uandresbello.edu>. Include the stage, your `INPUT.py` settings,
the full traceback, and the output of `python INPUT.py`.
