# 10 — Extending the pipeline: `User_Functions.py`

[← back to index](README.md)

Everything instrument-specific lives in **one file**: `User_Functions.py`. To support a new
instrument or a new data layout you write one function there and point `read_table` at it in
`INPUT.py`. Nothing else in the pipeline needs to change.

The file currently provides four functions:

| Function | Role |
|---|---|
| `TABLES_GIRAFFE_MOS_REDUCED(file)` | reader for GIRAFFE MOSSPECTRA files (Stages 1–2) |
| `TABLES_UVES_MOS_REDUCED(epoch, ARM='U')` | reader for UVES fibre-mode EsoReflex folders (Stages 1–2) |
| `READ_FITS_CONTINUUM(hdul)` | example spectrum reader for Stage 4 |
| `TEMPLATE_READER(path)` | example template reader for Stage 3 |

---

## 1. The metadata-reader contract (Stages 1 and 2)

```python
def MY_INSTRUMENT_READER(file_or_folder) -> dict:
    ...
    return DATA_FROM_FILE
```

It receives **the full path of one exposure** — a file if `FITS_FILE = True`, a folder (with a
trailing `/`) if `FITS_FILE = False` — and returns a dictionary. The authoritative key list is
`User_Functions.KEYS_FOR_TABLES_READER_FUNCTIONS`:

| Key | Required | Type / shape | Rules |
|---|---|---|---|
| `MJD` | ✔ | scalar float | Modified Julian Date **in days**. Convert with `astropy.time.Time` if your header stores something else. |
| `RAS` | ✔ | array `(N_fib,)` | Right ascension per fibre, **degrees**. |
| `DECS` | ✔ | array `(N_fib,)` | Declination per fibre, **degrees**. `len(RAS) == len(DECS)`. |
| `GLAT` | ✔ | scalar | Observatory latitude, degrees. |
| `GLON` | ✔ | scalar | Observatory longitude, degrees. |
| `HTEL` | ✔ | scalar | Observatory elevation, metres. |
| `WAVE` | ✔ | array `(N_pix,)` | Wavelength grid shared by all fibres of this exposure, in the instrument's native units. |
| `FLUX` | ✔ | array `(N_fib, N_pix)` | Must be a **numpy array**. `(N_pix, N_fib)` is also accepted — Stage 2 transposes automatically. |
| `NOT_SKY` | optional | bool array `(N_fib,)` | `True` = science target, `False` = sky fibre. **If omitted, no sky is subtracted.** |
| `OBJS` | optional | array `(N_fib,)` of str | Target names. If omitted, targets are numbered `0…N−1` and the output filenames use those numbers. |
| `SHAPE` | optional | anything printable | Shown in `Observing_blocks.txt` — handy for spotting inconsistent exposures. |
| `MAGS` | optional | array `(N_fib,)` | Magnitudes; colour-code the QA figure and populate the `MAG` header keyword. |
| `FIBER` | informational | array `(N_fib,)` | Fibre/retractor id. Read by the GIRAFFE reader; not consumed downstream. |

Only the first eight are strictly needed; the optional ones make the products far more useful,
and supplying `NOT_SKY` and `OBJS` is strongly recommended.

### Golden rules

1. **Return numpy arrays**, not lists, for `FLUX`, `RAS`, `DECS`, `NOT_SKY`, `MAGS` — they are
   used with boolean masks and fancy indexing.
2. **Keep the fibre ordering consistent** between `RAS`, `DECS`, `OBJS`, `MAGS`, `NOT_SKY` and
   the rows of `FLUX`. Everything downstream relies on positional correspondence.
3. **Do not pre-subtract the sky and do not shift wavelengths** — the pipeline does that.
4. Set `wave_to_angstrom` in `INPUT.py` to match the units you return in `WAVE`.

## 2. Worked example: adding a new instrument

Suppose a hypothetical multi-fibre spectrograph writes one FITS per exposure with

* HDU 0: header only,
* HDU 1 `SPECTRA`: image `(N_fib, N_pix)`, wavelengths in Å via `CRVAL1`/`CDELT1`,
* HDU 2 `FIBRES`: binary table with `NAME`, `ALPHA`, `DELTA`, `KIND` (`'SKY'` or `'OBJ'`), `VMAG`.

```python
# ── add this to User_Functions.py ───────────────────────────────────────────
def TABLES_MYSPEC(file):
    """Reader for MYSPEC multi-fibre exposures. Returns the standard dictionary."""
    D = dict()
    with fits.open(file) as F:
        # --- flux and wavelength -------------------------------------------
        flux  = np.array(F[1].data)                      # (N_fib, N_pix)
        CRVAL = F[1].header['CRVAL1']
        CDELT = F[1].header['CDELT1']
        NPIX  = F[1].header['NAXIS1']

        D['FLUX'] = flux
        D['WAVE'] = CRVAL + np.arange(NPIX)*CDELT        # Angstrom → wave_to_angstrom = 1

        # --- per-fibre metadata --------------------------------------------
        tab        = F[2].data
        D['OBJS']  = tab['NAME']
        D['RAS']   = tab['ALPHA']                        # degrees
        D['DECS']  = tab['DELTA']                        # degrees
        D['MAGS']  = tab['VMAG']
        D['NOT_SKY'] = (tab['KIND'] != 'SKY')            # True = science

        # --- epoch and site -------------------------------------------------
        D['MJD']  = Time(F[0].header['DATE-OBS'], format='isot').mjd
        D['GLAT'] = F[0].header['SITELAT']
        D['GLON'] = F[0].header['SITELONG']
        D['HTEL'] = F[0].header['SITEALT']

        D['SHAPE'] = flux.shape                          # optional, printed in the OB table
    return D
```

Then in `INPUT.py`:

```python
from User_Functions import TABLES_MYSPEC       # add to the imports at the top
...
FITS_FILE  = True
read_table = TABLES_MYSPEC
wave_to_angstrom = 1
```

and validate with `python INPUT.py`, which will call your reader on the first file and print
the mean RA/DEC, the MJD and the number of sky fibres. If those numbers look right, Stages 1
and 2 will work.

### If your exposures are folders

Set `FITS_FILE = False` (so `key = '*/'`) and write the reader to accept a folder path,
globbing inside it — that is exactly what `TABLES_UVES_MOS_REDUCED` does:

```python
table = glob.glob(epoch + f'*FIB_SCI_INFO_TAB_RED{ARM}.fits')[0]
files = sorted(f for f in glob.glob(epoch+'*.fits')
               if f'_MWXB_SCI_RED{ARM}' in f and 'ERR' not in f)
```

Note the two implicit assumptions in that pattern, which you should reproduce carefully for
your own data: **the sorted spectrum files must line up, one-to-one and in order, with the
valid rows of the fibre table**, and all fibres of an exposure must share one wavelength
solution.

### Passing extra arguments to a reader

`MAIN.py` calls `read_table(path)` with a single argument. If your reader needs options (like
the UVES `ARM`), bind them in `INPUT.py`:

```python
read_table = lambda p: TABLES_UVES_MOS_REDUCED(p, ARM='L')
# or
from functools import partial
read_table = partial(TABLES_UVES_MOS_REDUCED, ARM='L')
```

The *name* `read_table` is what the pipeline imports; what it points to is up to you.

## 3. Custom spectrum reader for Stage 4 (`read_fits`)

Signature: **receives an open `HDUList`, returns `(wave, flux)`**.

```python
def READ_FITS_CONTINUUM(hdul):
    try:                                    # binary-table layout
        wave = hdul[1].data['wave']
        flux = hdul[1].data['flux']
    except:
        try:
            wave = hdul[1].data['wavelength']
            flux = hdul[1].data['Flux']
        except:                             # 1-D IRAF layout
            flux  = hdul[0].data
            CRVAL = hdul[0].header['CRVAL1']
            CDELT = hdul[0].header['CDELT1']
            wave  = CRVAL + np.arange(len(flux))*CDELT
    return wave, flux
```

Set `read_fits = READ_FITS_CONTINUUM` (or your own) in `INPUT.py`; `None` uses the equivalent
built-in `cont_fit.default_read`.

Useful variant — normalising already-normalised or multi-extension files:

```python
def read_fits(hdul):                       # as used in Fast_usage/Fast_continuum.py
    try:
        return hdul[-1].data.wave, hdul[-1].data.flux
    except:
        return hdul[-2].data.wave, hdul[-2].data.flux
```

(The `-1`/`-2` fallback exists because Stage 4 appends its own `CONTINUUM GUESS` extension to
the file the first time it opens it, which shifts what "the last extension" means.)

> The **automatic** continuum mode uses a different contract: `READ_FUNCTION(full_path)` →
> `(wave, flux, header)`. See `automated_cont_fit.default_read`.

## 4. Custom template reader for Stage 3 (`TEMPLATE_READER`)

Signature: **receives a path, returns `(wave, flux)`**.

```python
def TEMPLATE_READER(path):
    syn = np.genfromtxt(path, skip_header=47)
    return syn[:,0], syn[:,1]        # wavelength (Å), normalised flux
```

The bundled version reads AUTOKUR/ATLAS9 output, whose columns are
`wavelength | normalised flux | line flux | continuum flux`.

Two things to set for your own templates:

1. **`skip_header` must match your file.** It is a parameter, not a fixed property: AUTOKUR and
   ATLAS9 write headers of different lengths depending on version and options, and the bundled
   templates carry 37–41 comment lines. Note that `np.genfromtxt` also strips `#` comments
   independently, so a `skip_header` larger than the real header removes that many *data* rows
   from the blue end (with the default `47` and the 41-line HR14b template, the first 6 rows).
   That is inconsequential when the template comfortably exceeds the observed range — the `CCF`
   class checks the coverage and tells you if it does not — but set it to your own header
   length, or drop the argument and rely on `comments='#'`, if your template is only marginally
   wider than your spectra.
2. ⚠️ `MAIN.py` calls `TEMPLATE_READER` **by name**, not through the `reads_template` variable.
   If you write a different reader, either name it `TEMPLATE_READER` in `User_Functions.py`
   and import it in `INPUT.py`, or bypass the reader entirely by setting `WAVE_TEMP` and
   `FLUX_TEMP` in `INPUT.py`:

```python
w, f = my_own_reader('/path/to/template.fits')
WAVE_TEMP, FLUX_TEMP = w, f
```

## 5. Extending further

| You want to… | Where to look |
|---|---|
| change the OB positional tolerance (0.5″) or the equal-count requirement | `List_files_related.OBSERVING_BLOCK` |
| implement distance-weighted sky subtraction | `Merging_and_substract_sky.MERGING_MULTIPLE_SPECTRA` — the normalised inverse-distance array `dist` is already computed there but unused |
| change the output grid oversampling (currently ×3) or the ±50 km s⁻¹ padding | `MERGING_MULTIPLE_SPECTRA` (`3 * N_pixels`) and `extend_linear_range(..., rv_kms=50)` |
| add error propagation | `MERGING_MULTIPLE_SPECTRA` (UVES `ERR_*` files are currently skipped) |
| change the CCF definition (e.g. normalise by the template norm) | `CCF.cross_correlation` — a commented "template-match" variant is already there |
| add more analytic fits to the CCF peak | `CCF.find_rv` (the four `try/except` blocks) and `CCF.get_errs_info` |
| change the continuum guess (shift width, number of samples, percentile) | `cont_fit.__init__` (`self.shift`, `self.Nsamples`, `self.percentile_of_cont`, `self.percentile_range`) or the corresponding arguments of `automated_continuum` |
| add a photometric system to the isochrone tools | the `photometric_system` method of either isochrone class — copy the `'GAIA'` block and edit the filter names and A_λ/A_V coefficients |

Because `INPUT.py` is imported with `from INPUT import *`, **any** name you define there is
visible to `MAIN.py` and to `Getting_RV.py` (which is `exec`'d in the same namespace). That is
the intended extension point for new options: define them in `INPUT.py`, consume them in your
modified module.
