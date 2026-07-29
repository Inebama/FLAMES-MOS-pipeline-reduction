# 04 — Data organisation: inputs, outputs and how the stages chain

[← back to index](README.md)

This document describes **what the pipeline expects to find on disk**, **what it produces**,
and **how each stage locates the products of the previous one**. For the byte-level structure
of each file (HDUs, columns, keywords) see [05 — File formats](05_file_formats.md).

---

## 1. The mental model

The pipeline converts a *detector-centric* organisation into a *target-centric* one:

```
INPUT                                     OUTPUT
one file per exposure                     one file per star per observing block
all fibres inside (science + sky)          sky already subtracted
instrument wavelength solution             common oversampled grid, barycentric, Ångström
no grouping                                grouped in OB-0, OB-1, …
```

Everything is driven by two paths — `path_input` and `path_output` — and by a **reader
function** that knows how to open your particular flavour of data.

---

## 2. Input layout A — GIRAFFE MOSSPECTRA (one FITS per exposure)

Set `FITS_FILE = True` (→ `key = '*.fit*'`) and `read_table = TABLES_GIRAFFE_MOS_REDUCED`.

```
path_input/
├── ADP.2018-11-23T16_36_07.288.fits
├── ADP.2018-11-23T16_36_07.289.fits
├── ADP.2018-11-23T16_36_07.290.fits
├── ADP.2018-11-23T16_40_50.318.fits
├── ADP.2018-11-23T16_40_50.319.fits
└── ADP.2018-11-23T16_40_50.320.fits
```

Each file (ESO Phase-3 `MOSSPECTRA` product):

| HDU | Type | Content |
|---|---|---|
| 0 | `PrimaryHDU` | Image `(NAXIS2, NAXIS1) = (4861, 124)` = (wavelength pixels, fibres), plus the full ESO header |
| 1 | `BinTableHDU` `FIBER_SETUP` | 124 rows × 18 columns, one row per fibre |

Header keywords the reader uses:

| Keyword | Used for |
|---|---|
| `CRVAL2`, `CDELT2`, `NAXIS2` | wavelength solution along axis 2: `λ = CRVAL2 + i·CDELT2` (**nanometres**) |
| `MJD-OBS` | epoch of the exposure |
| `HIERARCH ESO TEL GEOLAT` / `GEOLON` / `GEOELEV` | observatory location, for the barycentric correction |

`FIBER_SETUP` columns the reader uses:

| Column | Used for |
|---|---|
| `OBJECT` | target name → output filename and `OBJECT` keyword |
| `RA`, `DEC` | fibre sky position (deg) → OB consistency test, fibre matching, output header |
| `MAGNITUDE` | optional; colour-codes the QA figure and is stored as `MAG` |
| `TYPE` | fibre class; `'S'`/`'s'` ⇒ **sky**, anything else ⇒ science (`NOT_SKY`) |
| `RP` | retractor/fibre id; read into `FIBER` (informational) |

> Note the transposition. `F[0].data` has shape `(4861, 124)`, i.e. *pixels × fibres*, whereas
> the rest of the pipeline wants *fibres × pixels*. `MERGING_MULTIPLE_SPECTRA` detects this
> automatically: it tries `Flux[0][NOT_SKY[0]]`, and if that fails it transposes every flux
> array. You do not need to do anything, but it explains the `try/except` at the top of that
> function.

## 3. Input layout B — UVES fibre mode / EsoReflex (one folder per exposure)

Set `FITS_FILE = False` (→ `key = '*/'`) and `read_table = TABLES_UVES_MOS_REDUCED`.

```
path_input/
├── UVES.2018-10-19T00_18_33.188/
│   ├── Epoch1_2_FIB_SCI_INFO_TAB_REDU.fits     ← fibre table, upper arm
│   ├── Epoch1_2_FIB_SCI_INFO_TAB_REDL.fits     ← fibre table, lower arm
│   ├── Epoch1_2_MWXB_SCI_REDU.fits             ← fibre 0, upper arm   ┐
│   ├── Epoch1_2_MWXB_SCI_REDU_1.fits           ← fibre 1              │ one 1-D
│   ├── …                                                              │ spectrum
│   ├── Epoch1_2_MWXB_SCI_REDU_6.fits           ← fibre 6              ┘ per file
│   ├── Epoch1_2_ERR_MWXB_SCI_REDU*.fits        ← error arrays (ignored)
│   └── … the same set for REDL …
└── UVES.2018-10-19T00_34_19.963/
    └── …
```

`TABLES_UVES_MOS_REDUCED(epoch, ARM='U')`:

1. Globs `*FIB_SCI_INFO_TAB_RED{ARM}.fits` and reads the first match.
   * Rows with an empty `TYPE` are dropped (the pipeline emits blank rows).
   * `TYPE == 'S'`/`'s'` ⇒ sky fibre.
   * `RA`, `DEC`, `OBJECT` come from the same table; `MJD-OBS` and the `HIERARCH ESO TEL GEO*`
     keywords from its primary header.
2. Globs `*_MWXB_SCI_RED{ARM}*.fits`, **excluding** anything containing `ERR`, and sorts by
   name. Each file's primary HDU is a 1-D spectrum with its own `CRVAL1`/`CDELT1`/`NAXIS1`
   (**Ångström**).

> **Critical assumption:** the alphabetically sorted spectrum files must correspond, in order,
> to the valid rows of the fibre table. With the ESO naming this holds (`…REDU.fits` sorts
> before `…REDU_1.fits` because `.` < `_` in ASCII), but verify it for your own data.
>
> The reader also returns `WAVE` as the wavelength array of the **last** file read, on the
> assumption that all fibres of one exposure share the same solution — which is true within an
> exposure but not necessarily between exposures calibrated on different nights. That case is
> handled at merge time by `MISMATCH_SECOND_DIMENSION` (see [07](07_stage2_sky_and_merging.md)).

No `MAGNITUDE` column exists for UVES, so `MAGS` is not returned; Stage 2 falls back to using
the median flux per fibre for the colour scale of the QA figure.

## 4. What a reader function must return

Any reader you plug into `read_table` must return a `dict` with these keys
(the contract is documented in `User_Functions.KEYS_FOR_TABLES_READER_FUNCTIONS`):

| Key | Required | Shape / meaning |
|---|---|---|
| `MJD` | ✔ | scalar; modified Julian date of the exposure (days) |
| `RAS` | ✔ | array `(N_fibres,)`, degrees |
| `DECS` | ✔ | array `(N_fibres,)`, degrees |
| `GLAT` | ✔ | scalar; observatory latitude (deg) |
| `GLON` | ✔ | scalar; observatory longitude (deg) |
| `HTEL` | ✔ | scalar; observatory elevation (m) |
| `WAVE` | ✔ | array `(N_pixels,)` (shared by all fibres of the exposure) |
| `FLUX` | ✔ | array `(N_fibres, N_pixels)` — or `(N_pixels, N_fibres)`, auto-detected |
| `NOT_SKY` | optional | boolean array `(N_fibres,)`, `True` = science target. If absent, *everything* is treated as science and no sky is subtracted. |
| `OBJS` | optional | array `(N_fibres,)` of target names. If absent, targets are numbered. |
| `SHAPE` | optional | anything printable; shown in `Observing_blocks.txt` |
| `MAGS` | optional | array `(N_fibres,)` of magnitudes; used for the QA colour scale and stored in headers |
| `FIBER` | informational | fibre/retractor id (read by the GIRAFFE reader, not used downstream) |

`MAIN.merge_dicts()` then stacks the per-exposure dictionaries into arrays of objects, in the
fixed order `['MJD','RAS','DECS','NOT_SKY','GLAT','GLON','HTEL','WAVE','FLUX','OBJS','SHAPE','MAGS']`,
inserting `None` for keys no reader supplied.

Worked example of writing a new reader: [10 — Extending](10_extending_user_functions.md).

---

## 5. Output tree

```
path_output/
│
├── Observing_blocks.npy      ← Stage 1: dict {'OB-0': array_of_filenames, 'OB-1': …}
├── Observing_blocks.txt      ← Stage 1: same content, PrettyTable format
│
├── OB-0/
│   │
│   ├── <OBJ>=<RA>=<DEC>=<DATE>.fits    ← Stage 2: one per science fibre
│   ├── Oficial_list.list               ← Stage 2: the authoritative file list
│   ├── INFO_FILES.txt                  ← Stage 2: file · RA · DEC · MAG · counts
│   ├── Substract_SKY.png               ← Stage 2: all merged spectra
│   ├── SKY.png                         ← Stage 2: sky fibres vs master sky
│   ├── skies/
│   │   └── SKY_<MJD>.fits              ← Stage 2: one per exposure
│   │
│   ├── OUTPUT_RV.txt                   ← Stage 3: RV catalogue for this OB
│   ├── CCFS.png                        ← Stage 3: summary figure
│   └── RESTFRAME/
│       ├── <same names>.fits           ← Stage 3: rest-frame spectra
│       │                                  (+ cached CONTINUUM GUESS extension after Stage 4)
│       ├── Stored_params.npy           ← Stage 4 interactive: parameters per spectrum
│       ├── Continuum_sub/
│       │   └── <same names>.fits       ← Stage 4 interactive: normalised spectra
│       └── CONT_SUB/
│           └── CS-<same names>.fits    ← Stage 4 automatic: normalised spectra
│
└── OB-1/ …
```

### Naming conventions

| Product | Pattern | Notes |
|---|---|---|
| Observing block | `OB-<n>`, `n = 0,1,2,…` | Order of the **seed** exposure in time. |
| Merged spectrum | `<OBJECT>=<RA_deg>=<DEC_deg>=<ISO_datetime>.fits` | `=` is the separator on purpose: `name.split('=')[0]` recovers the target name (used by `cheating.py`). RA/DEC are the *unrounded* degree values of the fibre; the date is the ISO form of the **mean MJD of the OB**, so all files of one OB share the same date string. |
| Sky file | `SKY_<MJD>.fits` | One per **exposure** (not per OB), MJD in full precision. |
| Automatic continuum | `CS-<original name>` | The `CS-` prefix is why `CONT_SUB/` cannot be used as `EXTRA_FOLDER_RV`, while `Continuum_sub/` can. |

> Filenames can be long (target + two full-precision coordinates + timestamp) — around 80
> characters. This is fine on ext4/APFS/NTFS but can hit limits when nested deeply on some
> network shares. Keep `path_output` short.

---

## 6. How the stages chain

The mechanism is deliberately simple and file-based; there is no database and no state kept in
memory between stages.

| Stage | Reads | Writes | Chaining variable |
|---|---|---|---|
| 1 | `path_input + key` | `path_output/Observing_blocks.{npy,txt}` | — |
| 2 | `Observing_blocks.npy` + the raw files again | `path_output/OB-n/…` | — |
| 3 | `path_output/OB-n/Oficial_list.list`, spectra in `OB-n/ + EXTRA_FOLDER_RV` | `…/RESTFRAME/`, `OUTPUT_RV.txt`, and RV keywords **written back into the input spectra** | `EXTRA_FOLDER_RV` |
| 4 | `path_output/OB-n/Oficial_list.list`, spectra in `OB-n/ + EXTRA_FOLDER_CONT` | `…/Continuum_sub/` or `…/CONT_SUB/`, `Stored_params.npy`, and a cached `CONTINUUM GUESS` extension **written back into the input spectra** | `EXTRA_FOLDER_CONT` |

Two consequences worth internalising:

1. **`Oficial_list.list` is always read from the OB root**, even when the spectra themselves
   live in a sub-folder. That is why any sub-folder used as `EXTRA_FOLDER_*` must preserve the
   original filenames.
2. **Stages 3 and 4 modify their input files in place** (adding header keywords / an extra
   extension). Nothing is destroyed — the original data arrays are untouched — but if you want
   a pristine copy, snapshot the OB folder before running them.

### The two supported orders

```
RV first (default)                     Continuum first
OB-n/                                  OB-n/
└── RESTFRAME/                         └── Continuum_sub/
    └── Continuum_sub/                     └── RESTFRAME/

EXTRA_FOLDER_RV   = ''                 EXTRA_FOLDER_RV   = 'Continuum_sub/'
EXTRA_FOLDER_CONT = 'RESTFRAME/'       EXTRA_FOLDER_CONT = ''
```

---

## 7. Supporting data that is not produced by the pipeline

| Item | Where | Format |
|---|---|---|
| Synthetic templates | `TESTING_FOLDER/SYNTH_TEMPLATE/*.dat` | AUTOKUR/ATLAS9 ASCII: a `#`-commented header (37–41 lines) then 4 columns — wavelength (Å), normalised flux, line flux, continuum flux. `TEMPLATE_READER` returns columns 0 and 1. |
| PARSEC isochrones | user-supplied | CMD 3.x ASCII tables. `Header=13` for the old format, `Header=14` for PARSEC v2 / CMD 3.8. Read with `np.genfromtxt(..., names=True, skip_header=Header, comments='#')`, so the code needs the column-name line to be the first line after the skipped header. |
| Photometric catalogue | user-supplied | Any FITS/table you can load into `MAG`, `COLOR`, `OBJ`, `RA`, `DEC` arrays. |

> `TEMPLATE_READER`'s `skip_header` is a parameter you set to match your template format —
> AUTOKUR/ATLAS9 header lengths vary with version and options, and the bundled templates carry
> 37–41 comment lines. See [10 §4](10_extending_user_functions.md#4-custom-template-reader-for-stage-3-template_reader).
