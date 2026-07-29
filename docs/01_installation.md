# 01 — Installation and requirements

[← back to index](README.md)

The package is a collection of plain Python scripts. There is nothing to compile, nothing to
install, and no package metadata: you clone the repository and run the scripts from inside
`Package_Data_Reduction/`.

---

## 1. Python version

Developed and tested with **Python 3.11.7**.

The practical minimum is **Python 3.10**, because the code uses PEP-604 type unions in
signatures and annotations, e.g.

```python
wave_range : float|tuple[float,float] = (-1, np.inf)
min_wave   : float|None               = None
```

These are evaluated at definition time in several modules, so 3.9 and older will raise
`TypeError: unsupported operand type(s) for |`.

## 2. Required packages

| Package | Author's version | Used for |
|---|---|---|
| `numpy` | 1.26.4 | arrays, masks, percentiles, `.npy` persistence |
| `scipy` | 1.11.4 | `interp1d`, `UnivariateSpline`, `curve_fit` |
| `astropy` | 5.3.4 | `io.fits`, `time.Time`, `coordinates.SkyCoord`/`EarthLocation`/`Angle`, `units` |
| `matplotlib` | 3.9.0 | all figures; `widgets.Slider`/`Button` for the two GUIs; `gridspec` |
| `prettytable` | 3.10.0 | human-readable tables in `Observing_blocks.txt` and the isochrone output |

```bash
python -m pip install "numpy>=1.24" "scipy>=1.10" "astropy>=5.2" "matplotlib>=3.7" prettytable
```

`tkinter` (standard library) is used **optionally**, for the "choose a folder / choose a file"
dialogs. If it is unavailable the pipeline falls back to a numbered prompt on the terminal, so
a missing `tkinter` is not fatal for `MAIN.py`. It *is* required by the two isochrone tools,
which call `tkinter.filedialog` directly at import time.

> On Linux, `tkinter` often ships separately: `sudo apt install python3-tk`.
> On macOS with conda: `conda install tk`. On the system Python from python.org it is included.

### Optional dependencies (stand-alone tools only)

| Tool | Extra requirement |
|---|---|
| `Fast_usage/barycentric_vels.py` | `astroquery` (SIMBAD name resolution) |
| `plotingRVS.py` | Both optional: **`gls`** — the generalised Lomb–Scargle of Zechmeister & Kürster (2009), A&A 496, 577; <https://github.com/mzechmeister/GLS> (Python package: <https://github.com/mzechmeister/GLS/tree/master/python>). Either `pip install gls` or set `GLS_PATH` to the folder holding `gls.py`; if absent the script falls back to `astropy.timeseries.LombScargle`, and if that is missing it just skips the period search. **`scikit-learn`** — only for the fixed-bandwidth Epanechnikov kernel density; silently skipped if absent. |
| `Review_Multiple_Epochs.py` | none beyond the core list. Needs an interactive matplotlib backend (`TkAgg`) like the other GUIs. |
| `multiple_isochrone_fitting.py` | none beyond the core list, but it uses `multiprocessing` and benefits strongly from many cores |

## 3. Matplotlib backend

Both interactive components — the continuum GUI (`Continuum_fitting.py`) and the isochrone
GUIs — rely on live `Slider`/`Button` widgets and on blocking `plt.show()` calls.

`MAIN.py` contains this guard before opening the Tk folder picker:

```python
if matplotlib.get_backend().lower() != 'tkagg':
    raise Exception('There is a bug that whenever another backend is being initialized '
                    'this may crash so better to check it by terminal')
```

If the backend is not `TkAgg`, the exception is caught and the pipeline falls back to asking
you to type the OB number on the terminal — the pipeline still works, only the graphical
folder picker is skipped.

Recommendations:

* Run the interactive stages from a **normal terminal**, not from Spyder/PyCharm's embedded
  console, not from Jupyter, and not over an X-less SSH session.
* Force the backend if needed:
  ```bash
  MPLBACKEND=TkAgg python MAIN.py
  ```
* The non-interactive stages (1, 2, 3, and the *automatic* continuum) work fine with `Agg`.

## 4. Getting the code

```bash
git clone https://github.com/Inebama/FLAMES-MOS-pipeline-reduction.git
cd FLAMES-MOS-pipeline-reduction/Package_Data_Reduction
```

All commands in this documentation assume that `Package_Data_Reduction/` is the working
directory. `INPUT.py` locates the package automatically with

```python
path_package = __file__.replace('INPUT.py', '')
```

so relative paths inside `TESTING_FOLDER` resolve regardless of where you launch from — but
for your own data always give **absolute paths**.

## 5. The self-test: `python INPUT.py`

`INPUT.py` is a configuration file, but it also has an `if __name__ == '__main__'` block that
runs a set of checks. **Running it never processes data and never writes anything.**

```bash
python INPUT.py
```

What it does, depending on which stage flags you enabled:

| Condition | Check performed |
|---|---|
| always | Imports every package the pipeline needs. An `ImportError` here means a missing dependency. |
| `LIST` or `MERGE` | Globs `path_input + key`, calls your `read_table()` on the **first** file, and prints mean RA, mean DEC, MJD, number of sky fibres and the filename. |
| `MERGE` | Inspects the shape of `OUTPUT['FLUX']` and exits with a warning if `AXIS_TO_MERGE` looks wrong for that shape. |
| `RADIAL_VEL` | Either confirms that you passed `WAVE_TEMP`/`FLUX_TEMP` directly, or calls `TEMPLATE_READER(complete_path_template)` and verifies that wavelength and flux have equal length. |
| `CONTINUUM` | Only if you set a custom `read_fits`: opens `path_output` with `astropy.io.fits` and checks that your reader returns equal-length arrays. |

Expected output for the bundled test data with `LIST = True`:

```
 mean RA	 mean DEC	 DATE	 How Many Skyes?	 First_File_Found
6.1507...  -72.0798...  58411.1588...  8  .../DATA_GIRAFFE/ADP.2018-11-23T16_36_07.288.fits

If this match what you expect (check please) the LIST and MERGE method should work fine

DONE tests, ...  everything should work well
```

The `How Many Skyes?` column is computed as `sum(~OUTPUT['NOT_SKY'])`, i.e. the number of sky
fibres found in that exposure. For the bundled GIRAFFE data this is **8** (of 124 fibres, so
116 science targets).

> ⚠️ The `CONTINUUM` branch of the self-test opens `path_output`, which is a *directory*, not a
> file. It therefore only makes sense — and only runs — if you have supplied a custom
> `read_fits`; with the default `read_fits = None` the check is skipped. See
> [14 — Troubleshooting](14_troubleshooting.md).

## 6. Running the pipeline

```bash
python MAIN.py
```

`MAIN.py` does `from INPUT import *` and then executes whichever of the four stages you set to
`True`. There is a commented-out block at the top of `MAIN.py` that lets you drive the stage
flags from the command line instead (`-list`, `-merge`, `-continuum`, `-ravel`); uncomment it
if you prefer that style.

## 7. Disk-space and performance notes

* Stage 2 resamples onto a grid ~3× finer than the native pixel sampling and extended by
  ±50 km s⁻¹, then writes **one FITS file per science fibre per OB**. For the bundled GIRAFFE
  test data that is 116 files × 2 OBs, each holding 14 713 double-precision rows × 2 columns
  (≈ 240 kB each) → ≈ 55 MB per OB. Plan accordingly for large programmes.
* Stage 3 is the CPU-heavy stage: the default 3-pass scan evaluates the CCF
  100 + 100 + 500 = 700 times per spectrum, each evaluation being one `interp1d` call over the
  template plus one dot product.
* Stage 4 (interactive) deliberately keeps **only one spectrum in RAM at a time**; each
  spectrum is re-read from disk when displayed. The expensive continuum guess is cached as an
  extra FITS extension (`CONTINUUM GUESS`) inside the rest-frame file the first time it is
  computed, so revisiting a star is instantaneous.
* `multiple_isochrone_fitting.py` loads the **entire isochrone grid into RAM** and can use
  several GB with large PARSEC grids; inference takes of order seconds per star.

## 8. Platform notes

* Developed on Linux; the bundled paths in the `__main__` blocks of several stand-alone
  scripts still point at the author's machine (`/home/ian/...`, `/scratch/home/ibaeza/...`)
  and must be edited before use. They do not affect `MAIN.py`.
* On macOS, `TkAgg` figures sometimes open behind the terminal window — bring them to the
  front manually.
* Paths must **always end with `/`** wherever the documentation says so; the code concatenates
  strings rather than using `os.path.join`.
