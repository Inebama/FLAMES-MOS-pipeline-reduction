# 12 — Isochrone fitting: stellar parameters from a CMD

[← back to index](README.md)

Two independent, interactive tools that infer stellar parameters (T_eff, log g, mass,
luminosity, age, [M/H], …) by placing **PARSEC isochrones** on an observed colour–magnitude
diagram. They are **not** part of the four-stage spectroscopic pipeline: they are not called by
`MAIN.py`, take no settings from `INPUT.py`, and can be used on any photometric catalogue.

| Module | Method | Use it when |
|---|---|---|
| `single_isochrone_fitting.py` | one isochrone; each star inherits the parameters of the **geometrically nearest** isochrone point | you have a well-determined single stellar population and want a quick, transparent answer |
| `multiple_isochrone_fitting.py` | a **grid** of isochrones; full Bayesian inference with a Mahalanobis-distance likelihood and error propagation to every isochrone quantity | literature values for age/metallicity disagree, or no single isochrone reproduces the whole CMD |

The motivation for the second tool is astrophysical. In 47 Tuc — the system this was written
for — standard PARSEC isochrones match the main sequence and subgiants at higher metallicity
while the RGB/AGB drift towards lower-metallicity tracks, so **no single model reproduces the
whole CMD**. Rather than forcing a global (age, [M/H]) fit, the method abandons the "one
population" goal and infers *individual* stellar parameters, marginalising over a physically
motivated set of models.

---

# Part A — `multiple_isochrone_fitting.py` (Bayesian)

## A.1 Formalism

For a model point `θ` (one row of one isochrone) and an observed star `D` (one point in the
CMD),

$$P(\theta\,|\,D)\;\propto\;\mathcal{L}(D\,|\,\theta)\;\cdot\;P(\theta)$$

### The prior — the user's visual fit

The user selects, interactively, a **range of metallicities and ages** and places the grid on
the CMD by adjusting the distance modulus and the visual absorption until the isochrones match
the features stellar evolution says they must reproduce: the turn-off, the RGB bump and the
horizontal branch. Formally,

$$P(\theta)=\begin{cases}1 & \text{if the model was selected}\\ 0 & \text{otherwise}\end{cases}$$

This restricts the model space to a finite, physically motivated subset and simultaneously
fixes the distance modulus and reddening. It is a deliberate design choice: the astronomer's
knowledge of the system enters as a hard prior instead of as extra free parameters.

### The likelihood — generalised Mahalanobis distance

With `Δ = (Δmag, Δcol)ᵀ` the offset between the observed point and a model point,

$$D_M^2=\Delta^{T}\,\Sigma^{-1}\,\Delta \qquad\qquad \log\mathcal{L}\;\propto\;-\tfrac{N}{2}\,D_M^2$$

`N` (`degree` in the code, default `2.0`) tunes how fast the likelihood falls with distance:
`N = 2` is the Gaussian case and is the recommended value; larger `N` privileges the nearest
model points over densely sampled regions.

### The covariance matrix

For magnitude `m₁` and colour `c = m₂ − m₃`, with `A_x/A_V` the extinction coefficients:

$$\sigma^2_m=\sigma^2_{m_1,\rm phot}+\sigma^2_\mu+\left(\frac{A_{m_1}}{A_V}\delta A_V\right)^2$$

$$\sigma^2_c=\sigma^2_{m_2,\rm phot}+\sigma^2_{m_3,\rm phot}+\left(\frac{A_{m_2}-A_{m_3}}{A_V}\delta A_V\right)^2$$

$$\mathrm{Cov}(m,c)=f\,\sigma^2_{m_1,\rm phot}+\frac{A_{m_1}}{A_V}\cdot\frac{A_{m_2}-A_{m_3}}{A_V}\,\delta A_V^2$$

The covariance term encodes two real correlations: **(i)** a filter shared between the
magnitude and the colour, and **(ii)** the coherent shift that reddening and distance impose on
both axes simultaneously. The error ellipse is therefore tilted, which matters a lot along the
RGB where the isochrones run nearly parallel to that tilt.

The shared-filter factor `f` is set from the extinction coefficients, which identify which
filters are the same:

```python
if   A_mag1 == A_mag2: factor =  1     # magnitude also appears as the BLUE side of the colour
elif A_mag1 == A_mag3: factor = -1     # magnitude also appears as the RED  side of the colour
else:                  factor =  0     # no filter in common
```

For the built-in `'ubvir'` system (`Mag1 = I`, colour `V − I`, so `m₁ ≡ m₃`) this gives
`f = −1`: an error in `I` moves the star *up* in magnitude and *left* in colour, an
anti-correlation that tilts the error ellipse along the RGB. For `'gaia'` (`G` vs `BP − RP`) no
filter is shared and `f = 0`, leaving only the reddening term.

Defaults (representative of ~5 % flux errors, 1 % distance error, δE(B−V) ≈ 0.015):

| Argument | Default | Meaning |
|---|---|---|
| `err_mag1` | 0.05 | photometric error on the CMD magnitude [mag] |
| `err_mag2`, `err_mag3` | 0.05 | photometric errors on the two colour filters |
| `err_dist` | 0.03 | distance-modulus error [mag] (`δ(m−M) ≈ 2.17·δd/d`) |
| `err_redd` | 0.05 | error in visual absorption `δA_V` [mag] |

Each may be a scalar (same for all stars) or an array (one per star); in the latter case the
code builds and inverts one covariance matrix per star automatically.

### Cone sampling — removing the model-density bias

A naive sum over all nearby model points would weight regions where the isochrones happen to
be densely sampled. Instead, the CMD neighbourhood around each star is divided into
`N_bins = 360` angular sectors (1° cones) and **only the closest model point in each cone**
contributes:

```python
theta = np.arctan2(delta_mag, delta_col)
# one point per angular bin: the one with the smallest radius
```

Each direction on the CMD therefore contributes exactly one measurement, and at most 360 model
points enter the posterior. Reducing `N_bins` coarsens the sampling and speeds things up.

### Propagation to any parameter

The normalised weights `ω_i = exp(logL_i − max logL) / Σ` are attached to the model points, so
**any** column of the isochrone can be marginalised, preserving the non-linear mapping between
the CMD and, say, the Kiel diagram:

$$\bar X=\frac{\sum_i \omega_i X_i}{\sum_i \omega_i}\qquad
\sigma_{X\pm}=\sqrt{\frac{\sum_i \omega_i (X_i-\bar X)^2}{\sum_i \omega_i}}$$

with the sums restricted to `X_i ≥ X̄` for `σ₊` and `X_i ≤ X̄` for `σ₋` — i.e. **asymmetric**
weighted dispersions, which is what you need when the posterior is multi-modal (a star's
CMD position can be compatible with RGB, HB and AGB solutions simultaneously). In the output
table `<key>_err_l` is stored as a **negative** number and `<key>_err_u` as positive.

### What `d_min` is, and why it is written that way

`log_likelihood_cmd` returns, as its second value,

```python
MAX = logL.max()          # attained at the SMALLEST D^2, since logL = -0.5*N*D^2
...
return weights, (MAX/(-0.5*N))**0.5, N_w
```

This is the exact **inverse of the likelihood transform**: it converts the maximum
log-likelihood back into the Mahalanobis distance of the closest model point,

`MAX = −(N/2)·D²_min` → `MAX/(−N/2) = D²_min` → `(…)**0.5 = D_M,min`.

Returning `MAX` itself would not do, for two reasons:

| `degree` N | `MAX = logL.max()` | `(MAX/(-0.5*N))**0.5` |
|---|---|---|
| 1 | −0.180 | 0.600000 |
| 2 | −0.360 | 0.600000 |
| 5 | −0.900 | 0.600000 |
| 10 | −1.800 | 0.600000 |

1. `MAX` is negative and scales with the tunable hyperparameter `degree`, so it is not
   comparable between runs; dividing by `−N/2` cancels `N` **exactly** (same geometry in every
   row above). The *weights* still depend on `N`, which is the intended effect of `degree`.
2. Eq. B.9 compares a **distance** against an allowed radius, so the square root is what makes
   `d_min` the right kind of object.

`d_min` is therefore a genuine, `degree`-independent Mahalanobis distance, and it is the true
minimum over the box-selected candidates: within a 1° cone the direction is essentially fixed
and `D² = r²·ûᵀΣ⁻¹û` is monotonic in `r`, so the per-cone Euclidean-nearest point is also the
Mahalanobis-nearest one.

`d_min` is consumed in exactly two places — the `d_min` output column and the
`Skepticism_score` line — and never re-enters the inference.

### Quality diagnostics

Two columns accompany every star. Neither feeds back into the inference — they exist so that
you can filter or re-weight the results afterwards.

**Skepticism score** — lower is better:

$$S=\frac{D_{M,\min}}{R}\;\Big(1-\frac{N_{\rm used}}{N_{\rm cones}}\Big)^{1/2},
\qquad R=\sqrt{\texttt{min\_mag\_dist}^2+\texttt{min\_color\_dist}^2}$$

```python
Skepticism_score = (d_min/((min_mag_dist**2+min_color_dist**2)**0.5)) * ((1-(NP_mod/N_bins))**0.5)
```

The first factor measures how close the star lies to the model grid, normalised by the radius
of the acceptance region (`R = 0.0828` mag with the defaults, so `1/R = 12.08`); values near
zero indicate a well-matched, densely sampled region, larger values poor coverage or
extrapolation. The second factor measures how completely the neighbourhood was sampled — see
[§ The coverage factor](#the-coverage-factor) below. Typical behaviour:

| case | D_M,min | N_used | S |
|---|---|---|---|
| well embedded | 0.05 | 355 | 0.071 |
| typical RGB star | 0.50 | 300 | 2.466 |
| near edge of grid | 1.50 | 355 | 2.136 |
| sparse coverage, close | 0.20 | 120 | 1.973 |
| outside the grid | 2.00 | 50 | 22.424 |

Because `S` is a product, either factor going to zero zeroes the score: a star surrounded by
models in all `N_bins` directions scores 0 whatever its distance.

**Scale error** — `scale_err = 1 + √S`, a heuristic factor for rescaling the asymmetric
uncertainties of stars that cannot be discarded but whose formal errors are clearly
under-estimated, typically because only a handful of model points contributed and the resulting
dispersion is artificially small. It has no formal statistical justification and is intended
for individual problem cases (a star lying outside the grid), not for blanket application to a
whole catalogue — at `S = 0.07` it already inflates errors by 26 %.

The author's usage, left as a note at the end of the module, is in TOPCAT column syntax:

```
exp10(logTe + logTe_err_U * pow(scale_err,1))     # rescaled upper Teff bound
logg + logg_err_L * pow(scale_err,1)              # rescaled lower log g bound
```

### The coverage factor

`N_used` (`NP_mod` in the output) is the number of **populated** cones: `line_of_sight_of_point`
does `continue` on any angular bin that contains no model point, so `idx` — and hence
`N_w = len(weights)` — has exactly one entry per direction that actually had a model point.
Empty cones therefore raise the score, which is the intended meaning: the factor measures how
*completely* the star is surrounded by models.

| situation | `N_used` | `(1 − N_used/N_cones)^½` |
|---|---|---|
| dense grid, star fully surrounded | 360/360 | 0.000 |
| star on the edge of the grid (models on one side only) | 180/360 | 0.707 |
| sparse grid, ~40 model points nearby | 39/360 | 0.944 |

Two consequences worth knowing, and one harmless edge case:

* **`N_used` depends on the pre-selection box.** Candidates are first cut to the rectangle
  `±max_mag_dist` × `±max_color_dist` (default `±1` mag × `±0.3` mag), which is strongly
  anisotropic: a cone along the colour axis is truncated at 0.3 mag while one along the
  magnitude axis reaches 1.0 mag. A cone can therefore be empty because the *box* is narrow in
  that direction, not because the grid is. `Skepticism_score` is only comparable between runs
  that share the same `max_*_dist` and the same `N_bins`.
* **The product form lets either factor zero the score.** A star sitting in a hole in the grid
  but surrounded by models in all 360 directions gets `N_used = N_bins` → `S = 0`, whatever
  `d_min` is. Rare, but the distance term cannot flag it.
* **Edge case:** `np.digitize` assigns `theta` exactly `= +π` to bin `N_bins+1`, which the loop
  `range(1, nbins+1)` does not cover, so such a point is dropped. It requires `Δmag == 0` and
  `Δcol < 0` exactly — no practical impact.

> The disabled weight filter in `log_likelihood_cmd`
> (`# mask = (weights <= 0)`; `N_w = len(weights)`) should **stay** commented out: enabling it
> would redefine `N_used` as "cones whose weight did not underflow to zero" instead of
> "populated cones", silently changing what the coverage factor measures.

**`good` flag**

| Value | Meaning |
|---|---|
| `1` | parameters were inferred |
| `0` | no model point inside the search box, or the star is farther than the minimum-distance criterion → all parameters `nan` |
| `-1` | initial/failure value |

The rejection criterion is

```python
(abs(d_mag/min_mag_dist) + abs(d_col/min_color_dist)).min() > 1.3
```

an L1 approximation of an ellipse (cheaper than the exact test) around the star, with
`min_mag_dist = 0.075`, `min_color_dist = 0.035` by default.

## A.2 Using it

```python
from multiple_isochrone_fitting import ISO_FIT

ISO_FIT(MAG   = I,              # CMD magnitude (m1)
        COLOR = V - I,          # CMD colour   (m2 - m3)
        OBJ   = ids,            # required in practice — see the caveat below
        RA    = ra,
        DEC   = dec,
        iso_path   = '/path/to/PARSEC_grid.txt',   # None → a file dialog opens
        photometry = 'ubvir',   # 'gaia' | 'gaia_old' | 'ubvir' | 'ubvir_old'
        Header     = 14,        # 13 for old PARSEC, 14 for PARSEC v2 / CMD 3.8
        degree     = 2.0,       # N in the likelihood
        N_bins     = 360,       # angular cones
        err_mag1=0.05, err_mag2=0.05, err_mag3=0.05, err_dist=0.03, err_redd=0.05,
        Stages     = [1,2,3,4,7,8])
```

### The GUI

```
┌──────────────────────────────────┬──────────────────────────┐
│                                  │ [ Open other Isochrone ] │
│         Colour–Magnitude         │                          │
│            Diagram               │   MH        logAge       │
│                                  │  ▓ range   ▓ range       │
│  · grey  = your stars            │  ▓ slider  ▓ slider      │
│  + blue  = isochrone points OUT  │                          │
│  + red   = isochrone points IN   │  MH range   Age range    │
│            the selected ranges   │  -0.78 dex  10.00 Gyr    │
│                                  │  -0.41 dex  13.00 Gyr    │
├──────────────────────────────────┤                          │
│  Visual Absorption [──●──] - +   │   ┌──────────────────┐   │
│  Distance Modulus  [───●─] - +   │   │     Obtain       │   │
│                                  │   │ stellar parameters│  │
│                                  │   └──────────────────┘   │
└──────────────────────────────────┴──────────────────────────┘
```

Workflow:

1. Open a PARSEC grid file (one file containing many ages/metallicities).
2. Use the two **range sliders** to restrict [M/H] and log(Age) to the literature-motivated
   interval — this *is* the prior. Selected points turn red.
3. Slide **Visual Absorption** and **Distance Modulus** (with `+`/`−` fine steps) until the
   selected isochrones match the turn-off, the RGB bump and the horizontal branch.
4. Press **Obtain stellar parameters**. The window closes and the computation starts, printing
   a percentage counter. All CPU cores are used (`multiprocessing.Pool` with
   `imap_unordered`).
5. `Isochrone_values.fits` is written in the **current working directory**.

### Output

See [05 §13](05_file_formats.md#13-isochrone-outputs). Every isochrone column appears three
times (`<key>`, `<key>_err_l`, `<key>_err_u`) plus `good`, `d_min`, `NP_mod`,
`Skepticism_score`, `scale_err`. The primary header records the isochrone file, the selected
[M/H] and log Age ranges, `AV`, `UM` and the extinction coefficients — everything needed to
reproduce the run.

`save_plain_txt()` writes the same table as a PrettyTable text file if you prefer ASCII (call
it instead of `save_fits` in `match_nearest`).

## A.3 Practical notes

1. **Subsetting hooks.** The first two lines of `__init__` let you restrict the catalogue
   without touching your input arrays:
   ```python
   mask = ()      # a boolean array here selects a subset; () means "everything"
   nnn  = None    # an integer here keeps only the first nnn entries
   self.MAG = MAG[mask][:nnn]
   ```
   Handy for testing the fit on a hundred stars before launching the full run.
2. **`OBJ`, `RA` and `DEC` are effectively required.** They are sliced in the constructor, so
   passing `None` raises `TypeError`. Supply dummy arrays if you have no identifiers.
3. **The initial `AV` and `UM` are hard-coded** to `0.01` and `13.37` (47 Tuc values) in
   `isochrone_read`, overriding the automatically estimated values computed just above. You
   adjust them with the sliders anyway, but do not read them as defaults appropriate to your
   system — change those two literals if you work on another cluster regularly.
4. **Memory.** The whole isochrone grid is loaded with `np.genfromtxt` and kept in RAM, and the
   arrays are re-sent to every worker process. Large PARSEC grids plus a large photometric
   catalogue can reach several GB; the GUI may lag on the first draw.
5. **Speed.** Of order seconds per star. With `imap_unordered` over all cores, a few thousand
   stars is minutes, not hours; but reduce `N_bins`, tighten `max_mag_dist`/`max_color_dist`,
   or restrict the CMD region if you need it faster.
6. **PARSEC-only.** The column names (`MH`, `logAge`, `label`, `logTe`, `logg`, the magnitude
   names) are hard-coded to the CMD web-interface output. Other isochrone families require
   real code changes.

---

# Part B — `single_isochrone_fitting.py` (nearest point)

## B.1 What it does

1. Reads one PARSEC file, keeps only the evolutionary `Stages` you asked for, and lists the
   available `MH` and `logAge` values.
2. Opens a GUI with **discrete** sliders for [M/H] and log Age (they step through the values
   actually present in the file) plus continuous sliders for `A_V` and the distance modulus.
   The initial `A_V`/`UM` are estimated from the mean offset between your data and the
   isochrone — a sane starting point.
3. On **Obtain stellar parameters**:
   * `increment_points()` densifies the current isochrone to ~10 000 points by linear
     interpolation of *every* column against point index, **splitting the isochrone wherever
     consecutive points jump by more than 4σ** of the typical spacing (so the phase transitions
     of the track are not bridged by spurious interpolated points), and distributing the new
     points proportionally to the arc length of each segment;
   * each star is assigned the parameters of the isochrone point minimising the plain Euclidean
     distance `√(Δmag² + Δcol²)` — no errors, no weighting, no marginalisation;
   * `Isochrone_values.txt` is written (PrettyTable with a `#` header recording the isochrone,
     `MH`, `logAGE`, `AV`, `UM` and the two extinction coefficients);
   * a **Kiel diagram** (T_eff vs log g) of the result pops up, with both axes inverted.

## B.2 Using it

```python
from single_isochrone_fitting import ISO_FIT

ISO_FIT(MAG   = I, COLOR = V - I,
        OBJ = ids, RA = ra, DEC = dec,      # all optional here
        iso_path   = '/path/to/one_isochrone.txt',
        photometry = 'ubvir',
        Header     = 14,
        Stages     = [1,2,3,4,7,8])
```

**Open other Isochrone** restarts the whole application with a new file, keeping your data
loaded.

> Running `python single_isochrone_fitting.py` directly does **not** launch the fit: its
> `__main__` block is the author's 47 Tuc analysis (differential-reddening map from Pancino
> 2024, extinction corrections, diagnostic plots) and it calls `sys.exit()` before reaching the
> `ISO_FIT(...)` line. Import the class instead, or delete that block in your own copy.
> `multiple_isochrone_fitting.py` has a similar personal `__main__` block, but there the
> `ISO_FIT(...)` call *is* reached.

---

# Part C — shared reference

## C.1 Photometric systems

Set with `photometry=` (case-insensitive). Each entry defines which isochrone columns are
`Mag1` (the CMD magnitude), `Mag2` and `Mag3` (the colour), together with the extinction
coefficients `A_λ/A_V`:

| Key | Mag1 | Mag2 − Mag3 | Notes |
|---|---|---|---|
| `'GAIA'` | `G_fSBmag` | `G_BP_fSBmag − G_RP_fSBmag` | PARSEC v2 column names (DR3/EDR3) |
| `'GAIA_OLD'` | `Gmag` | `G_BPmag − G_RPmag` | PARSEC v1.2 column names |
| `'UBVIR'` | `I_fSBmag` | `V_fSBmag − I_fSBmag` | PARSEC v2; used with Stetson photometry (Stetson et al. 2019) |
| `'UBVIR_OLD'` | `Imag` | `Vmag − Imag` | PARSEC v1.2 |
| `'system_random'` | `magnitude1` | `magnitude2 − magnitude3` | template to copy for your own system |

Coefficients in use: `A_U/A_V = 1.55814`, `A_B = 1.32616`, `A_V = 1.00096`, `A_R = 0.80815`,
`A_I = 0.59893`, `A_J = 0.28688`, `A_H = 0.18103`, `A_K = 0.11265`; `A_G = 0.83627`,
`A_BP = 1.08337`, `A_RP = 0.63439` (PARSEC v2 CMD 3.8 output page, R_V = 3.1, Cardelli et al.
1989 / O'Donnell 1994 for a G2V star).

**To add a system**, copy the `'GAIA'` block in `photometric_system` and edit the column names
and the coefficients. The keys of `absorption_coeff` must match the names used in `system`.
The class then derives

```python
abs_mag = A[Mag1]                  # applied to the magnitude axis
abs_col = A[Mag2] - A[Mag3]        # applied to the colour axis
```

and the isochrone is placed with

```python
mag_fit   = mag_iso   + UM + AV*abs_mag
color_fit = color_iso +      AV*abs_col
```

## C.2 Evolutionary stages (PARSEC `label`)

| `label` | Stage | Default |
|---|---|---|
| 0 | pre-main sequence | excluded |
| 1 | main sequence | included |
| 2 | subgiant branch | included |
| 3 | red giant branch | included |
| 4 | core He burning (CHEB) | included |
| 5 | red CHEB | excluded |
| 6 | blue CHEB | excluded |
| 7 | early AGB | included |
| 8 | TP-AGB | included |
| 9 | post-AGB | excluded |

Change with the `Stages` argument.

## C.3 Preparing an isochrone file

Download from the [PARSEC CMD web interface](http://stev.oapd.inaf.it/cmd) as a single ASCII
file:

* choose the photometric system matching your `photometry` key;
* for `multiple_isochrone_fitting.py`, request a **grid**: a range of ages and a range of
  metallicities in one file — bracket the literature values for your system;
* for `single_isochrone_fitting.py` a single (or small) set is enough;
* keep the `label` column (evolutionary stage) and the physical columns you care about
  (`logTe`, `logg`, `logL`, `Mass`, …);
* set `Header=14` for PARSEC v2 / CMD 3.8 output, `Header=13` for the older format. The
  parameter is `skip_header` for `np.genfromtxt(..., names=True)`, so the **column-name line
  must be the first line after the skipped block**. If you get `KeyError: 'MH'` or
  `KeyError: 'label'`, your `Header` is off by one — count the leading `#` lines in your file.

## C.4 Choosing between the two tools

| | `single_isochrone_fitting` | `multiple_isochrone_fitting` |
|---|---|---|
| Models used per star | 1 isochrone, 1 point | up to `N_bins` points from a whole grid |
| Uncertainties | none | full covariance → asymmetric weighted dispersions |
| Handles multimodality (RGB/HB/AGB degeneracy) | no | yes — that is the point |
| Speed | instantaneous | seconds per star, parallelised |
| Output | text (PrettyTable) | FITS binary table + quality flags |
| Best for | a quick look, a clean single population | publication-grade individual parameters |
