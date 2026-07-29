# Documentation — FLAMES–MOS Pipeline Reduction

This folder contains the complete reference documentation for the package. It is
self-contained: every algorithm, file format, configuration variable and auxiliary tool of the
repository is described here.

## How to read this

**If you have never run the package** → read in this order:

1. [01 — Installation and requirements](01_installation.md)
2. [02 — Quick start tutorial](02_quickstart.md)
3. [04 — Data organisation](04_data_organization.md)
4. [03 — `INPUT.py` reference](03_input_reference.md)

**If you want to understand what the code actually computes** →
[06](06_stage1_observing_blocks.md), [07](07_stage2_sky_and_merging.md),
[08](08_stage3_radial_velocities.md), [09](09_stage4_continuum.md).

**If you want to use only one engine (CCF or continuum) on your own spectra** →
[11 — Manual mode](11_manual_mode.md).

**If you want to plug in a new instrument** → [10 — Extending](10_extending_user_functions.md).

**If you came here for the isochrone tools** → [12 — Isochrone fitting](12_isochrone_fitting.md).

**If something breaks** → [14 — Troubleshooting](14_troubleshooting.md).

## Index

| File | Contents |
|---|---|
| [01_installation.md](01_installation.md) | Python and package requirements, matplotlib backend, the `INPUT.py` self-test, platform notes |
| [02_quickstart.md](02_quickstart.md) | Stage-by-stage tutorial on the bundled GIRAFFE test dataset, with the expected output at every step |
| [03_input_reference.md](03_input_reference.md) | Complete reference of every variable in `INPUT.py`, grouped by stage, with types, defaults and effects |
| [04_data_organization.md](04_data_organization.md) | What the input data must look like (GIRAFFE MOSSPECTRA, UVES/EsoReflex), the output directory tree, file-naming conventions, and how the stages chain through sub-folders |
| [05_file_formats.md](05_file_formats.md) | HDU-by-HDU description of every FITS product, all binary-table columns, all header keywords, and the format of every text/npy product |
| [06_stage1_observing_blocks.md](06_stage1_observing_blocks.md) | The Observing-Block grouping algorithm, the positional consistency test, and how to choose the time threshold |
| [07_stage2_sky_and_merging.md](07_stage2_sky_and_merging.md) | Fibre classification, master-sky construction, flux scaling, barycentric correction, wavelength resampling, co-addition, QA figures |
| [08_stage3_radial_velocities.md](08_stage3_radial_velocities.md) | The `CCF` class: range checks, cosmic-ray clipping, 3-pass scan, peak finder, the four analytic fits, the error budget, sign conventions, output plots |
| [09_stage4_continuum.md](09_stage4_continuum.md) | The shift-and-percentile continuum guess, the iteration semantics (−1 / 0 / 1 / >1), the interactive GUI manual, the automatic convergence mode |
| [10_extending_user_functions.md](10_extending_user_functions.md) | The reader-function contract, the required dictionary keys, worked example for a new instrument, custom template and FITS readers |
| [11_manual_mode.md](11_manual_mode.md) | `Fast_RV.py`, `Fast_continuum.py`, `barycentric_vels.py`, `stellar_oscillation_params.py`, `convective_turnover.py` |
| [12_isochrone_fitting.md](12_isochrone_fitting.md) | `single_isochrone_fitting.py` and `multiple_isochrone_fitting.py`: the Bayesian formalism, covariance matrix, cone sampling, weighted statistics, quality flags, outputs, performance |
| [13_auxiliary_tools.md](13_auxiliary_tools.md) | `plotingRVS.py` (multi-epoch RV review + periodograms), `Review_Multiple_Epochs.py` (interactive multi-epoch viewer), `cheating.py` (bookkeeping recipes) |
| [14_troubleshooting.md](14_troubleshooting.md) | Error messages, known bugs and quirks, workarounds, sanity checks |
| [15_api_reference.md](15_api_reference.md) | Signatures and return values of every public function and class |

## Conventions used in these documents

* **Wavelengths.** The pipeline is unit-agnostic internally: it works in whatever units the
  reader function returns (nm for GIRAFFE, Å for UVES). Only when writing FITS binary tables
  is the factor `wave_to_angstrom` applied, so the `wave` column of every output product is in
  **Ångström**. Stage 3 and Stage 4 therefore always see Ångström.
* **Velocities.** Always km s⁻¹; `c = 299 792.458` km s⁻¹ hard-coded throughout.
* **Times.** Always MJD (days).
* **Coordinates.** Degrees, FK5 by default (`coord_sys='fk5'`, the ESO convention).
* **Flux.** After Stage 2 the flux is *relative*: each spectrum is divided by its own median.
  The header keyword `MCOUNT` stores the factor that converts back to the original counts.
* `<OBJ>=<RA>=<DEC>=<DATE>.fits` denotes the standard output filename pattern of Stage 2.
* Code paths are given relative to `Package_Data_Reduction/`.
