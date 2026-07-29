# 05 — File formats: every product, HDU by HDU

[← back to index](README.md)

Complete reference of every file the pipeline writes. All FITS products are written with
`astropy.io.fits` and `overwrite=True`. All binary-table wavelength columns are in
**Ångström** (the conversion factor `wave_to_angstrom` is applied at write time).

Quick map:

| Product | Written by | Location |
|---|---|---|
| [`Observing_blocks.npy`](#1-observing_blocksnpy) | Stage 1 | `path_output/` |
| [`Observing_blocks.txt`](#2-observing_blockstxt) | Stage 1 | `path_output/` |
| [merged spectrum](#3-merged-sky-subtracted-spectrum) | Stage 2 | `OB-n/` |
| [`SKY_<MJD>.fits`](#4-sky_mjdfits) | Stage 2 | `OB-n/skies/` |
| [`Oficial_list.list`](#5-oficial_listlist) | Stage 2 | `OB-n/` |
| [`INFO_FILES.txt`](#6-info_filestxt) | Stage 2 | `OB-n/` |
| [`Substract_SKY.png`, `SKY.png`](#7-qa-figures) | Stage 2 | `OB-n/` |
| [rest-frame spectrum](#8-rest-frame-spectrum) | Stage 3 | `OB-n/RESTFRAME/` |
| [`OUTPUT_RV.txt`](#9-output_rvtxt) | Stage 3 | `OB-n/` (or `OB-n/EXTRA_FOLDER_RV/`) |
| [`CCFS.png`](#7-qa-figures) | Stage 3 | same as above |
| [`CONTINUUM GUESS` extension](#10-continuum-guess-cache) | Stage 4 | appended to the rest-frame files |
| [normalised spectrum](#11-normalised-continuum-subtracted-spectrum) | Stage 4 | `Continuum_sub/` or `CONT_SUB/` |
| [`Stored_params.npy`](#12-stored_paramsnpy) | Stage 4 (interactive) | next to the spectra |
| [`Isochrone_values.fits` / `.txt`](#13-isochrone-outputs) | isochrone tools | current working directory |

---

## 1. `Observing_blocks.npy`

A pickled Python `dict`, written with `np.save`. Read it back with:

```python
OB = np.load(path_output+'Observing_blocks.npy', allow_pickle=True).item()
# {'OB-0': array(['ADP.…318.fits', 'ADP.…319.fits', 'ADP.…320.fits'], dtype='<U32'),
#  'OB-1': array([...])}
```

* Keys: `'OB-0'`, `'OB-1'`, … in the order the blocks were seeded (chronological).
* Values: a 1-D array of **file/folder names relative to `path_input`** (not full paths).

## 2. `Observing_blocks.txt`

Human- and machine-readable transcription of the same information, one PrettyTable per block,
all lines prefixed with `#` so `np.genfromtxt` skips them by default.

```
# OBS0
#           File_Name              Date             Diff_in_Hour   Diff_in_Minutes   Shapes_of_Data
 ADP.2018-11-23T16_40_50.318.fits  58410.012881928  0.0            0.0               (4861,124)
 …
```

| Column | Meaning |
|---|---|
| `File_Name` | exposure file/folder |
| `Date` | MJD |
| `Diff_in_Hour`, `Diff_in_Minutes` | time offset from the **seed** exposure of the block |
| `Shapes_of_Data` | whatever the reader put in `SHAPE` (spaces stripped) — only present if the reader supplies it |

## 3. Merged, sky-subtracted spectrum

`OB-n/<OBJECT>=<RA>=<DEC>=<DATE>.fits` — written by `Merging_and_substract_sky.new_fits`.

```
No.  Name                 Type          Dimensions
  0  PRIMARY              PrimaryHDU    ()                 ← header only
  1  SPECTRUM (sky sub)   BinTableHDU   N rows × 2 cols
```

### Extension 1 — `SPECTRUM (sky sub)`

| Column | Format | Unit | Content |
|---|---|---|---|
| `wave` | `D` | Angstrom | Common wavelength grid, **barycentric-corrected**, ~3× oversampled relative to the native pixels and extended by ±50 km s⁻¹ (see [07](07_stage2_sky_and_merging.md)). |
| `flux` | `D` | relative flux | Sky-subtracted, per-spectrum median-scaled, co-added flux. Zeros where an exposure contributed no coverage. |

Extension header comments: `TCOMM1 = 'wavelenght in angstrom corrected by barycentric vel'`,
`TCOMM2 = 'relative Flux wiht sky substracted'`, plus `COMMENT Object: <OBJ>`.

### Primary header

| Keyword | Type | Meaning |
|---|---|---|
| `OBJECT` | str | Target name taken from the fibre table of the first exposure. |
| `RA` | float | Fibre right ascension, degrees. The comment field shows the sexagesimal form. |
| `DEC` | float | Fibre declination, degrees. |
| `MJD` | float | **Mean** MJD of the exposures in the block; the comment carries the ISO date. |
| `MAG` | float | Magnitude from the fibre table (only if the reader supplies `MAGS`). |
| `OFIL0…OFILn` | str | Name of each parent exposure file. |
| `DFIL0…DFILn` | float | MJD of each parent exposure. |
| `GEOLON` | float | Observatory longitude (deg). |
| `GEOLAT` | float | Observatory latitude (deg). |
| `GEOLEV` | float | Observatory elevation (m). |
| `MCOUNT` | float | Multiply `flux` by this to recover counts instead of relative flux. It is the combination (`function_used_to_merge`, across exposures) of the per-exposure median counts of this fibre. |
| `RVBARY` | float | Mean barycentric correction over the block, km s⁻¹. |
| `BARY0…BARYn` | float | Barycentric correction of each individual exposure, km s⁻¹. |
| `MJD0…MJDn` | float | MJD of each individual exposure (comment: ISO date). |
| `COMMENT` | str | `File done by I. Baeza`. |
| `RV`, `RVERR`, `RVTRUE` | float | **Added later by Stage 3** — see [§9](#9-output_rvtxt). |

> Header keyword length: `OFIL`/`DFIL`/`BARY`/`MJD` are numbered without padding, so blocks
> with more than 999 exposures would exceed the 8-character FITS keyword limit. Not a concern
> in practice.

## 4. `SKY_<MJD>.fits`

`OB-n/skies/SKY_<MJD>.fits` — one per **exposure**, written by
`Merging_and_substract_sky.save_sky`.

```
No.  Name          Type          Dimensions
  0  PRIMARY       PrimaryHDU    ()
  1  SKY-0         ImageHDU      (N_pixels,)     ← individual sky fibre 0
  2  SKY-1         ImageHDU      (N_pixels,)
  …                                              ← one per sky fibre
  k  Merged Sky    BinTableHDU   N rows × 2 cols ← the master sky
```

Primary header: `OBJECT = 'SKY'`, `MJD` (of that exposure), `OFIL` (parent file),
`GEOLON`/`GEOLAT`/`GEOLEV` (observatory longitude, latitude, elevation).

Each `SKY-i` image extension header carries `RA`, `DEC` of that sky fibre and an IRAF-style
linear solution `CRVAL1`, `CDELT1` in Ångström.

`Merged Sky` binary table: `wave` (Å) and `counts` (the median over the sky fibres — i.e. the
spectrum actually subtracted from the science fibres of that exposure).

> The sky spectra are **not** barycentric-corrected — they are stored in the observed frame,
> as they should be for sky diagnostics.

## 5. `Oficial_list.list`

Plain text, one filename per line, no header — the list of merged spectra of that OB, in the
order the fibres appear. Stages 3 and 4 read it with

```python
spectrums = np.genfromtxt(path_folder+'Oficial_list.list', dtype='U100')
```

so filenames must stay under 100 characters when you re-use this list.

## 6. `INFO_FILES.txt`

Tab-separated, one line per merged spectrum:

```
#File_name	RA	DEC	MAGS	N_counts
MED_RGB_107=…fits	5.994791666666566	-71.90288888888769	12.062000274658203	2831.107940673828
```

| Column | Meaning |
|---|---|
| `File_name` | as in `Oficial_list.list` |
| `RA`, `DEC` | degrees |
| `MAGS` | catalogue magnitude (or median flux if the reader gave no magnitudes) |
| `N_counts` | the `MCOUNT` scale factor |

Convenient for cross-matching targets across OBs (`cheating.py` recipe 2 does exactly that).

## 7. QA figures

| File | Stage | Contents |
|---|---|---|
| `OB-n/Substract_SKY.png` | 2 | All merged spectra of the OB overplotted, colour-coded by magnitude (or by median counts when magnitudes are unavailable), with a colour bar. X-range set by `min_wave`/`max_wave`; Y-range from the count scale. |
| `OB-n/SKY.png` | 2 | One row per exposure. Left: every individual sky fibre. Right: the master sky actually subtracted (black), labelled with the parent filename. |
| `OB-n/CCFS.png` | 3 | Left: all CCFs of the OB overplotted. Right, three stacked panels: all rest-framed spectra over the template (black) — full range, then two zoom windows. |
| `<name>_continuum.png`, `<name>_points_cont.png`, `<name>_converge.png` | 4 (automatic, `save_plots=True`) | Continuum guess per iteration; the accepted continuum points; and the convergence of the number of accepted points. Off by default. |

## 8. Rest-frame spectrum

`OB-n/RESTFRAME/<same name>.fits` — written by `Getting_RV.new_fits`.

```
No.  Name                    Type          Dimensions
  0  PRIMARY                 PrimaryHDU    ()             ← copy of the parent header, incl. RV*
  1  SPECTRUM (restframe)    BinTableHDU   N rows × 2 cols
  2  CONTINUUM GUESS         BinTableHDU   M rows × 3 cols ← appended later by Stage 4
```

### Extension 1 — `SPECTRUM (restframe)`

| Column | Format | Unit | Content |
|---|---|---|---|
| `wave` | `D` | Angstrom | `wave_parent / (1 + RV/c)` — the wavelength grid shifted to the stellar rest frame. |
| `flux` | `D` | adu | Unchanged flux from the parent file. |

The primary header is copied verbatim from the parent (so it already contains `RV`, `RVERR`,
`RVTRUE`, `MCOUNT`, `RVBARY`, `OFIL*`…), except that `CRVAL1`/`CDELT1` are removed if present:
after the Doppler shift the grid is no longer described by the parent's linear solution.

## 9. `OUTPUT_RV.txt`

Written once per OB by `Getting_RV.py`, in `OB-n/EXTRA_FOLDER_RV/`.

```
# Method used Synthetic
# Path to sythetic spectra: /…/GIRAFFE_TEMPLATE_HR14b_R36900.dat
# Radial velocities are in km/s
#Object	RA	DEC	RV	RV_ERR	Barycentric_RV	RV_corrected	template_used_in_CCF
MED_RGB_107	5.9947916…	-71.9028888…	81.4242903…	0.3840443…	-13.3408221…	68.0834682…	relative_to_Synth_Spectra
```

| # | Column | Content |
|---|---|---|
| 1 | `Object` | `OBJECT` keyword |
| 2 | `RA` | degrees |
| 3 | `DEC` | degrees |
| 4 | `RV` | CCF result **minus** the barycentric correction → topocentric-equivalent value. Mirrored in the FITS header as `RV`. |
| 5 | `RV_ERR` | total 1σ uncertainty (km s⁻¹), see [08](08_stage3_radial_velocities.md#5-the-error-budget) |
| 6 | `Barycentric_RV` | `RVBARY`, the mean barycentric correction of the block |
| 7 | `RV_corrected` | the CCF result itself. Because Stage 2 already put the wavelengths in the barycentric frame, **this is the physical radial velocity**. Mirrored in the header as `RVTRUE`. |
| 8 | `template_used_in_CCF` | `relative_to_<OBJ_Template>` |

With `method='Template'` two extra comment lines record the template star and its RV.

## 10. `CONTINUUM GUESS` cache

Appended by `Continuum_fitting.cont_fit.get_spectra` **to the input spectrum file** the first
time a spectrum is displayed, so the (expensive) statistical guess is computed only once.

| Column | Format | Unit | Content |
|---|---|---|---|
| `wave` | `D` | Angstrom | the wavelength grid **after** applying `WAVE_RANGE` |
| `init_mask` | `L` | BOOL | pixels selected as "continuum candidates" by the initial guess |
| `Cont_guess` | `D` | Relative Flux | the statistical continuum guess |

Extension header:

| Keyword | Meaning |
|---|---|
| `Nsamp` | number of velocity shifts sampled (2500) |
| `Seed` | numpy random seed (170499) — guarantees reproducibility |
| `Shift` | characteristic shift σ in km s⁻¹ (180 interactive / 200 automatic) |
| `wMin`, `wMax` | wavelength range actually used |
| `PerC` | percentile adopted as the continuum (86 interactive) |
| `PerL`, `PerU` | percentile band used to build `init_mask` (80, 95 interactive) |

## 11. Normalised (continuum-subtracted) spectrum

Both continuum modes write the **same 6-column structure**; they differ only in the folder,
the filename prefix and a few header keywords.

| Mode | Path | Written by |
|---|---|---|
| Interactive | `…/Continuum_sub/<same name>.fits` | `cont_fit.UPDATE_fits` |
| Automatic | `…/CONT_SUB/CS-<same name>.fits` | `automated_continuum.UPDATE_fits` |

```
No.  Name                   Type          Dimensions
  0  PRIMARY                PrimaryHDU    ()               ← parent primary header
  1  SPECTRUM (cont sub)    BinTableHDU   N rows × 6 cols
```

| Column | Format | Unit | Content |
|---|---|---|---|
| `wave` | `D` | Angstrom | wavelength (after `WAVE_RANGE` trimming) |
| `nflux` | `D` | Normalized flux | `flux / cont` — **the normalised spectrum** |
| `cont` | `D` | Relative flux | the continuum that was divided out |
| `flux` | `D` | Relative flux | input flux, rescaled by its median |
| `Cont_guess` | `D` | Relative Flux | the statistical initial guess |
| `init_mask` | `L` | BOOL | the continuum-point mask (final mask in automatic mode) |

Extension header — **interactive mode** records the exact GUI settings, so the fit is
reproducible:

| Keyword | Meaning |
|---|---|
| `SIGMAL`, `SIGMAU` | lower / upper sigma-clipping thresholds |
| `ITERS` | number of iterations (`-1`, `0`, `1` or `>1`; see [09](09_stage4_continuum.md)) |
| `SMOOTH` | spline smoothing factor `s` |
| `Nsamp`, `Seed`, `Shift`, `wMin`, `wMax`, `PerC`, `PerL`, `PerU` | as in the guess cache |

**Automatic mode** instead records:

| Keyword | Meaning |
|---|---|
| `medext` | the median flux that was divided out at the start (multiply back to recover the input scale) |
| `Nsamp`, `Seed`, `Shift`, `wMin`, `wMax` | as above |
| `PerC` | the **final** percentile of the ramp |
| `PerL`, `PerU` | `percentile_range` |

If the primary header contains `CRVAL1`/`CDELT1`, the interactive writer also stores `nflux`
as the primary image data, giving you an IRAF-readable 1-D file for free.

## 12. `Stored_params.npy`

Pickled `dict` mapping **filename → `[sigma_min, sigma_max, iterations, smooth]`**
(or `None` for spectra not yet fitted).

```python
p = np.load('…/RESTFRAME/Stored_params.npy', allow_pickle=True).item()
p['MED_RGB_107=…fits']      # → [1.25, 2.5, 0, 0.03]
```

* Saved on **every** exit path of `cont_fit` (including exceptions) — you never lose work.
* `cont_fit(..., Continue=True)` reloads it and reopens the GUI at the first spectrum whose
  entry is `None`. If all entries are filled it reopens at the last one.
* If the keys do not match the spectrum list exactly (different number or different names),
  the class prints a warning and exits — set `Continue=False` to start over.
* `store_cont_sub(Force=False)` refuses to write any output until every entry is non-`None`.
* This file is the unit of exchange used by `cheating.py` recipe 1 to copy one OB's
  normalisation parameters onto the other OBs.

## 13. Isochrone outputs

### `multiple_isochrone_fitting.py` → `Isochrone_values.fits`

```
No.  Name     Type          Content
  0  PRIMARY  PrimaryHDU    configuration of the fit
  1  Data     BinTableHDU   one row per star
```

Primary header: `Iso` (isochrone file), `MH_l`/`MH_u` and `logA_l`/`logA_u` (the selected
ranges), `AV`, `UM` (visual absorption and distance modulus adopted), `AbsMag`, `AbsCol`
(the extinction coefficients of the photometric system).

Table columns: `OBJ`, `RA`, `DEC`, `MAG`, `COLOR`, then, **for every column present in the
PARSEC isochrone** (`Mass`, `logL`, `logTe`, `logg`, `MH`, `logAge`, every magnitude, …):

| Column | Meaning |
|---|---|
| `<key>` | posterior weighted mean |
| `<key>_err_l` | lower asymmetric weighted dispersion — stored as a **negative** number |
| `<key>_err_u` | upper asymmetric weighted dispersion |

plus the quality diagnostics:

| Column | Meaning |
|---|---|
| `good` | `1` = inferred, `0` = star too far from the model grid (skipped), `-1` = default/failure |
| `d_min` | minimum Mahalanobis distance to the model grid |
| `NP_mod` | number of model points that contributed (≤ `N_bins`) |
| `Skepticism_score` | reliability metric; **lower is better** |
| `scale_err` | suggested multiplicative rescaling of the uncertainties |

Definitions and caveats: [12 — Isochrone fitting](12_isochrone_fitting.md).

### `single_isochrone_fitting.py` → `Isochrone_values.txt`

A `#`-commented PrettyTable. The header block records the isochrone path, `MH`, `logAGE`,
`AV`, `UM` and the two extinction coefficients; the table gives, per star, the input `MAG`
and `COLOR` (plus `Obj`/`RA`/`DEC` when supplied) followed by **every column of the nearest
isochrone point**.

Both files are written to the **current working directory**, not to `path_output`.
