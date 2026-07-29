# 02 — Quick start tutorial

[← back to index](README.md)

This tutorial walks through a complete reduction of the **bundled GIRAFFE test dataset**, from
raw ESO products to normalised rest-frame spectra. Every step tells you what to change, what
to run, and what you should see.

The golden rule of this package:

> **You edit `INPUT.py`. You run `MAIN.py`. You do not touch anything else.**
> Never rename a variable or a function in `INPUT.py` — the other modules import them by name.
> Only change their *contents*.

---

## 0. The test dataset

```
Package_Data_Reduction/TESTING_FOLDER/
├── DATA_GIRAFFE/                 6 exposures, GIRAFFE HR14b, 47 Tuc field
│   ├── ADP.2018-11-23T16_36_07.288.fits   MJD 58411.1588   ┐
│   ├── ADP.2018-11-23T16_36_07.289.fits   MJD 58411.1698   ├ night 2
│   ├── ADP.2018-11-23T16_36_07.290.fits   MJD 58411.1808   ┘
│   ├── ADP.2018-11-23T16_40_50.318.fits   MJD 58410.0128   ┐
│   ├── ADP.2018-11-23T16_40_50.319.fits   MJD 58410.0238   ├ night 1
│   └── ADP.2018-11-23T16_40_50.320.fits   MJD 58410.0348   ┘
├── DATA_UVES/                    2 EsoReflex UVES-fibre exposure folders
├── SYNTH_TEMPLATE/               4 ATLAS9/AUTOKUR synthetic templates
└── TEST_OUTPUT/
    └── Test_example_try_to_reproduce_it/    ← reference results
```

Each GIRAFFE file is a `MOSSPECTRA` product: a primary image of shape `(4861, 124)`
(4861 wavelength pixels × 124 fibres) plus a `FIBER_SETUP` binary table with one row per
fibre. Of the 124 fibres, **8 are sky** (`TYPE == 'S'`) and **116 are science targets**.

The reference output was produced with `Range_in_hours_to_be_consider_same_Observing_Block =
24`, which yields exactly two Observing Blocks (one per night).

> The reference numbers may differ slightly from what you obtain — they were generated with an
> earlier version of the code (differences of order 0.5 km s⁻¹ in RV, or a few points in
> spectrum length, are expected and documented by the author).

---

## 1. Sanity check before touching anything

```bash
cd Package_Data_Reduction
python INPUT.py
```

Out of the box `INPUT.py` ships with

```python
path_input  = path_package + 'TESTING_FOLDER/DATA_GIRAFFE/'
path_output = path_package + 'TESTING_FOLDER/TEST_OUTPUT/Test_your_attemp/'

LIST       = True
MERGE      = False
RADIAL_VEL = False
CONTINUUM  = False
```

so this first run already tests exactly what we need. You should see the mean RA/DEC of the
field, the MJD of the first file, and `8` sky fibres, followed by

```
DONE tests, ... everything should work well
```

If you get an `ImportError`, install the missing package ([01 — Installation](01_installation.md)).
If `read_table` raises, your data is not in the expected format ([10 — Extending](10_extending_user_functions.md)).

---

## 2. Stage 1 — define the Observing Blocks

Relevant `INPUT.py` settings:

```python
LIST      = True
FITS_FILE = True     # search '*.fit*' ; set False to search subfolders '*/' (UVES/EsoReflex)
Range_in_hours_to_be_consider_same_Observing_Block = 24
read_table = TABLES_GIRAFFE_MOS_REDUCED
SAVE_TXT   = True
```

Run:

```bash
python MAIN.py
```

The terminal prints one PrettyTable per OB, and the output folder is created with:

```
Test_your_attemp/
├── Observing_blocks.npy
└── Observing_blocks.txt
```

`Observing_blocks.txt` should look like:

```
# OBS0
#           File_Name              Date             Diff_in_Hour   Diff_in_Minutes   Shapes_of_Data
 ADP.2018-11-23T16_40_50.318.fits  58410.012881928  0.0            0.0               (4861,124)
 ADP.2018-11-23T16_40_50.319.fits  58410.023881644  0.26399…       15.8396…          (4861,124)
 ADP.2018-11-23T16_40_50.320.fits  58410.034862726  0.52754…       31.6523…          (4861,124)

# OBS1
 ADP.2018-11-23T16_36_07.288.fits  58411.158871859  0.0            0.0               (4861,124)
 ADP.2018-11-23T16_36_07.289.fits  58411.169848420  0.26344…       15.8062…          (4861,124)
 ADP.2018-11-23T16_36_07.290.fits  58411.180835398  0.52712…       31.6275…          (4861,124)
```

**Checkpoint.** Two OBs, three exposures each, identical data shapes. If your science case
needs finer time sampling, lower the threshold (e.g. `0.25` h groups individual exposures) and
re-run — this stage is cheap and non-destructive.

The `.npy` file is a plain dict and can be inspected at any time:

```python
import numpy as np
OB = np.load('…/Observing_blocks.npy', allow_pickle=True).item()
print(OB.keys(), OB['OB-0'])
```

Algorithm details: [06 — Stage 1](06_stage1_observing_blocks.md).

---

## 3. Stage 2 — sky subtraction and merging

```python
LIST  = True      # can be False once Observing_blocks.npy exists; harmless to leave True
MERGE = True

function_used_to_merge = np.nanmedian   # robust; np.nanmean converges faster
min_wave = None                          # x-range of the QA figure, in *native* units
max_wave = None                          #   (e.g. 655 / 657.5 nm for HR14b)
wave_to_angstrom = 10                    # GIRAFFE pipeline gives nm → ×10 = Å
AXIS_TO_MERGE = 0
```

> `wave_to_angstrom` is the single most important instrument-dependent number here. GIRAFFE
> MOSSPECTRA wavelengths are in **nanometres** (`CRVAL2 = 638.3`, `CDELT2 = 0.005`), so the
> factor is 10. UVES products are already in Ångström → factor 1.

```bash
python MAIN.py
```

For each OB you will see the per-exposure science/sky shapes printed, then:

```
Merging and substracting sky  OB-0
(116, 8)
(116, 4861)	(8, 4861)
…
Saving Figure
```

Result:

```
Test_your_attemp/OB-0/
├── MED_RGB_102=6.345874999999894=-71.91016666666545=2018-10-19T00:34:22.837.fits
├── … 116 files total, one per science fibre …
├── Oficial_list.list      116 lines — the authoritative list for stages 3 and 4
├── INFO_FILES.txt         file · RA · DEC · MAG · median counts
├── Substract_SKY.png      all 116 merged spectra, colour-coded by magnitude
├── SKY.png                individual sky fibres and the master sky, per exposure
└── skies/
    ├── SKY_58410.012881928.fits
    ├── SKY_58410.023881644.fits
    └── SKY_58410.034862726.fits
```

**Checkpoint — look at the two PNGs.**

* `SKY.png`: left column = the individual sky fibres of each exposure; right column = the
  master sky actually subtracted (black). They should look like the same sky spectrum;
  a wildly discrepant fibre means a contaminated sky fibre.
* `Substract_SKY.png`: all merged spectra overplotted, coloured by magnitude. Emission
  residuals at the sky-line positions, or spectra sitting at very different levels, point at a
  problem. Use `min_wave`/`max_wave` to zoom on a diagnostic region.

Filename anatomy: `<OBJECT>=<RA_deg>=<DEC_deg>=<ISO date of the mean MJD>.fits`. The `=`
separator is deliberate — it lets you recover the target name with `name.split('=')[0]`,
which several helper scripts do.

Algorithm details: [07 — Stage 2](07_stage2_sky_and_merging.md).

---

## 4. Stage 3 — radial velocities

```python
RADIAL_VEL = True
Across_OB  = True          # loop over every OB automatically

EXTRA_FOLDER_RV = ''       # look directly inside OB-n/  (i.e. RV before continuum)

complete_path_template = path_package + \
    'TESTING_FOLDER/SYNTH_TEMPLATE/GIRAFFE_TEMPLATE_HR14b_R36900.dat'
reads_template = TEMPLATE_READER

RADVEL_RANGE = (-250, 250)
MIN_WAVE_CCF = None        # ⚠ leave at least one of these as None — see docs/14
MAX_WAVE_CCF = None
OBJ_Template = 'Synth_Spectra'

PLOT_EACH_ONE     = False  # True → pops up the full diagnostic plot for every star
SAVE_SUMMARY_PLOT = True
```

The bundled template is an ATLAS9/AUTOKUR synthetic spectrum with
`TEFF = 4100 K, log g = 1.00, [M/H] = −0.60, [α/Fe] = +0.4`, degraded to `R = 36 900` — a
sensible representative of a 47 Tuc red giant observed with HR14b.

```bash
python MAIN.py
```

You get a live table:

```
OBJ               RV       	 RV_ERR	 BarycenRV	 RV_corr  	relative_to_What?
MED_RGB_107        81.424	 0.384	  -13.341	   68.083	relative_to_Synth_Spectra
MED_RGB_3283        0.092	 0.384	  -13.329	  -13.237	relative_to_Synth_Spectra
…
```

and, per OB:

```
OB-0/
├── OUTPUT_RV.txt        the RV catalogue (8 columns, see below)
├── CCFS.png             summary: all CCFs + all rest-framed spectra vs the template
└── RESTFRAME/
    └── … 116 rest-frame FITS files …
```

`OUTPUT_RV.txt` columns:

```
#Object  RA  DEC  RV  RV_ERR  Barycentric_RV  RV_corrected  template_used_in_CCF
```

**Which column is the physical radial velocity?** `RV_corrected` (column 7). Stage 2 already
moved the wavelengths into the barycentric frame, so the CCF result is *already* barycentric.
Column 4 (`RV`) is that value minus the barycentric correction, i.e. the topocentric-equivalent
number. The same two numbers are written into each FITS header as `RVTRUE` and `RV`
respectively. `plotingRVS.py` reads `RV_corrected`. Full sign discussion:
[08 — Stage 3](08_stage3_radial_velocities.md#6-sign-conventions).

**Checkpoint — open `CCFS.png`.** Every CCF should have a single, clean, well-separated peak,
and in the three spectrum panels the stellar lines of all stars should line up with the black
template. Stars whose lines are systematically offset are either genuine RV outliers (field
stars, binaries) or failed correlations. Set `PLOT_EACH_ONE = True` to inspect one at a time.

Stage 3 also **modifies the Stage-2 files in place**, adding `RV`, `RVERR` and `RVTRUE` to
their primary headers. This is intentional (it keeps the RV attached to the spectrum) but
means Stage 3 is not read-only.

---

## 5. Stage 4 — continuum normalisation

### 5a. Interactive (recommended)

```python
CONTINUUM         = True
AUTOMATED_CONT    = False
EXTRA_FOLDER_CONT = 'RESTFRAME/'    # normalise the rest-frame spectra
WAVE_RANGE        = 0.04            # trim 4 % of the range at each end
read_fits         = None            # use the built-in reader
```

```bash
python MAIN.py
```

Because `Across_OB` handling only applies to the RV stage, the interactive continuum stage
asks you **which OB to work on** — either through a Tk folder dialog or, if Tk is
unavailable/backend is wrong, through a numbered prompt:

```
Accordingly to the following list:
0 OB-0
1 OB-1

 Please select which Observing block do you want to work with (e.g 0):
```

A window opens with:

* the spectrum (blue), the statistical **initial guess** (green, always visible for reference),
  and the **current continuum** (red);
* a vertical slider on the left listing the spectra, plus *Previous*/*Next* buttons;
* four horizontal sliders — *Sigma below*, *Sigma above*, *Iterations*, *Smoother factor* —
  each with `+`/`−` fine-adjust buttons.

Iteration semantics (this is the key to using the tool efficiently):

| `Iterations` | Behaviour |
|---|---|
| `-1` | Accept the statistical initial guess as-is. All other sliders ignored. |
| `0` | Smooth the initial guess with a cubic spline. Only *Smoother factor* matters. |
| `1` | Fit a cubic spline to the points inside the initial acceptance band. Only *Smoother factor* matters. |
| `> 1` | Classical asymmetric sigma-clipping, starting from the initial band, `n` times. All sliders matter. |

Work through every star, then **close the figure window** — that is how you end the program.
Parameters are saved continuously to `RESTFRAME/Stored_params.npy`, so you can stop and resume
(the class is instantiated with `Continue=True`, which resumes at the first star that has no
stored parameters).

When you close the window, `MAIN.py` calls `CONT.store_cont_sub(Force=False)`, which refuses
to write unless **every** spectrum has stored parameters. Finish the OB, or edit the call to
`Force=True` if you deliberately want partial output.

Result:

```
OB-0/RESTFRAME/
├── Stored_params.npy
└── Continuum_sub/
    └── … 116 normalised FITS files (same names) …
```

Each output holds `wave`, `nflux` (normalised), `cont` (the continuum used), `flux`
(pre-normalisation), plus the initial guess and its mask, and the four fitting parameters in
the extension header — the fit is fully reproducible.

### 5b. Automatic

```python
CONTINUUM      = True
AUTOMATED_CONT = True
```

No GUI, no OB prompt: every OB is processed in a loop. The algorithm iterates the
shift-and-percentile estimate (raising the percentile from 84 to 92 across iterations) until
the number of accepted continuum points stabilises, then smooths with a cubic spline.

Result: `OB-n/RESTFRAME/CONT_SUB/CS-<original name>.fits` (note the `CS-` prefix).

> The automatic mode is good but not infallible: be careful near the spectrum edges and around
> very strong/broad lines, and always compare against the interactive result on a few stars.

Algorithm details, including the figures behind the method: [09 — Stage 4](09_stage4_continuum.md).

---

## 6. Choosing the order: RV first or continuum first?

Both orders are supported; you select them with the two "extra folder" variables.

| Order | `EXTRA_FOLDER_RV` | `EXTRA_FOLDER_CONT` | Resulting tree |
|---|---|---|---|
| **RV → continuum** (default) | `''` | `'RESTFRAME/'` | `OB-n/RESTFRAME/Continuum_sub/` |
| **continuum → RV** | `'Continuum_sub/'` | `''` | `OB-n/Continuum_sub/RESTFRAME/` |

Notes for the second order:

* The interactive continuum writes files with the **same names**, so `Oficial_list.list` still
  resolves — this is why `Continuum_sub/` works as `EXTRA_FOLDER_RV` while `CONT_SUB/` (which
  prefixes `CS-`) does not.
* `Getting_RV.py` reads the columns `wave` and `flux` from extension 1. In a continuum-
  subtracted file `flux` is the **pre-normalisation** flux, so the CCF would run on the
  un-normalised data. If you want the CCF to use `nflux`, run the RV stage manually with a
  custom reader (see [11 — Manual mode](11_manual_mode.md)).

---

## 7. Comparing with the reference output

```bash
diff <(cut -f1,4 TESTING_FOLDER/TEST_OUTPUT/Test_your_attemp/OB-0/OUTPUT_RV.txt) \
     <(cut -f1,4 TESTING_FOLDER/TEST_OUTPUT/Test_example_try_to_reproduce_it/OB-0/OUTPUT_RV.txt)
```

Small differences are expected (the reference was made with an older version). What *should*
match: the OB grouping, the number and names of the output files, the barycentric corrections
(to ~10⁻³ km s⁻¹), and the RVs to better than ~1 km s⁻¹.

---

## 8. Checklist for your own data

1. Choose or write a `read_table` function ([10 — Extending](10_extending_user_functions.md)).
2. Set `FITS_FILE` (`True` = files, `False` = folders per exposure, as in EsoReflex/UVES).
3. Set `wave_to_angstrom` for your instrument.
4. Get a template that **fully covers** your spectral range after shifting by
   `RADVEL_RANGE` (otherwise the CCF class will trim your observed spectrum, and tell you so).
5. Run stage by stage, look at every QA figure before enabling the next stage.
