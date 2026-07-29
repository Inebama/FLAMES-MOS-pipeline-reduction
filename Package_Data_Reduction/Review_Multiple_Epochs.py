# This file was written by AI (Claude Opus, Anthropic), inspired by the core of
# I. Baeza's `Joining_Multiple_epochs.py` visualiser. The layout, the reference/residual idea,
# the sigma bands and the line-bisector diagnostic are the author's; the AI rebuilt them as a
# clean, dataset-independent tool. The original script is untouched.
"""
Review_Multiple_Epochs.py -- interactive review of multi-epoch spectra
=====================================================================

An eyeball tool for the output of the FLAMES/MOS pipeline (or any set of per-epoch 1-D
spectra): step through your stars and see, for each one, every epoch overplotted together
with a residual against a reference of your choice.

    +--------------------------------------------------+---------------------+
    |  spectrum: all epochs + reference (black)        |  parameter panel    |
    |                                                  |  e.g. Teff vs logL  |
    +--------------------------------------------------+  (or RA vs DEC)     |
    |  residual: flux / reference, with sigma bands    |  click to jump      |
    +--------------------------------------------------+---------------------+
    |  [<] [ star slider ] [>]     (arrow keys work too)                     |
    +------------------------------------------------------------------------+

Everything is configured in the CONFIGURATION block. Reading your own auxiliary data
(temperatures, luminosities, ...) is your job: build a plain dict {star_name: value} and
assign it to PARAM_X / PARAM_Y. If you do not, the parameter panel shows RA/DEC instead.

Designed not to break in the awkward cases:
  * a single epoch                -> residual is flat, sigma bands are skipped
  * no template / no reference    -> reference falls back to the mean of the epochs
  * a star missing from some OBs  -> that epoch is simply not drawn
  * spectra on different grids    -> resampled onto a common grid per star

Controls
--------
    slider / [<] [>] / left-right arrow keys   change star
    click a point in the parameter panel       jump to that star
    'f'                                        fit the line profile in FIT_WINDOW
    'b'                                        toggle the line bisector
    's'                                        save the current figure to PNG
    'h'                                        print this help on the terminal
"""

# ---------------------------------------------------------------------------------------------
#  Packages
# ---------------------------------------------------------------------------------------------
import os
import sys
import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as U

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


# =============================================================================================
#  CONFIGURATION
# =============================================================================================

# --- Where the data is ------------------------------------------------------------------------
PATH = './TESTING_FOLDER/TEST_OUTPUT/Test_example_try_to_reproduce_it/'   # must end with '/'

# Observing blocks (epochs) to show.
#   None -> read them from PATH+'Observing_blocks.npy'
#   list -> exactly these, e.g. ['OB-0','OB-1','OB-4']
LIST = None

# Sub-folder inside each OB holding the spectra ('' = the OB folder itself).
SUBFOLDER = '/RESTFRAME/Continuum_sub/'

# Filename prefix added by some stages (the automatic continuum writes 'CS-<name>').
FILE_PREFIX = ''

# Wavelength range to display. None on either side = no limit.
W_MIN, W_MAX = None, None

# Resample every epoch onto a common grid of this many points (per star).
#   None -> use the first available epoch's own grid
N_SAMPLE = None

# --- What to divide by ------------------------------------------------------------------------
# 'mean'   : mean of the epochs shown for that star  (default; always available)
# 'median' : median of the epochs
# '<OB>'   : a specific epoch used as reference, e.g. 'OB-0'
# 'template': the spectrum in TEMPLATE_PATH (read with READ_TEMPLATE below)
REFERENCE = 'mean'
TEMPLATE_PATH = ''      # only used when REFERENCE == 'template'

# Re-normalise each spectrum before comparing. A callable f(wave, flux) -> float, or None.
#   None                                        -> leave the flux as read
#   lambda w, f: np.median(f)                   -> median
#   lambda w, f: np.percentile(f, 90)           -> high percentile (good for crowded spectra)
NORMALISE = lambda w, f: np.median(f)

# --- Which stars ------------------------------------------------------------------------------
# Match the same star across epochs by 'name' (the OBJECT part of the filename) or 'coords'.
MATCH_MODE          = 'name'
MATCH_RADIUS_ARCSEC = 1.75

# Drop these stars entirely (by name).
DISCARD_STARS = [
    # 'MED_RGB_370',
]

# Drop single epochs: {'OB name': [star, star, ...]}, or '*' for the whole epoch.
DISCARD_OBS = {
    # 'OB-2': ['MED_RGB_370', 'MED_RGB_402'],   # bad observations of these two stars
    # 'OB-5': '*',                              # skip this epoch entirely
}

# Order in which you step through the stars:
#   'name' | 'ra' | 'param_x' | 'param_y' | 'n_epochs' | None (file order)
SORT_BY = 'name'

# --- Parameter panel (top right) ---------------------------------------------------------------
# Read your own catalogue however you like and build dicts {star_name: value}.
# Leave both as None to show RA/DEC instead.
#
#   e.g.  import numpy as np
#         t = np.genfromtxt('params.txt', names=True, dtype=None, encoding=None)
#         PARAM_X = dict(zip(t['OBJ'], t['Teff']))
#         PARAM_Y = dict(zip(t['OBJ'], t['logL']))
PARAM_X, PARAM_Y = None, None
PARAM_X_LABEL, PARAM_Y_LABEL = r'$T_{\rm eff}$ [K]', r'$\log(L/L_\odot)$'
PARAM_X_INVERT, PARAM_Y_INVERT = True, False     # True for Teff (hot on the left), etc.
PARAM_BACKGROUND = None                          # optional (x, y) arrays drawn in grey behind

# --- Diagnostics you can switch on -------------------------------------------------------------
SHOW_SIGMA_BANDS = True          # +-1/3/5 sigma envelope on the residual panel
SMOOTH_N         = 0             # boxcar width in pixels for display only (0 = off)

# Vertical reference lines: [(wavelength, label), ...]
MARK_LINES = [
    (6562.819, r'H$\alpha$'),
]

# Shaded regions of interest / known trouble: [(w1, w2, label, colour), ...]
HIGHLIGHT = [
    # (6563.6, 6565.3, 'sky', 'red'),
    # (6552.0, 6554.0, 'sky', 'red'),
]

# Line-profile fit (key 'f'). Window around FIT_CENTER used for the fit.
FIT_CENTER = 6562.819
FIT_WINDOW = 5.0                 # +- angstrom
FIT_ON_START = False             # fit automatically for every star

# Line bisector (key 'b'); same window as the fit.
BISECTOR_ON_START = False
BISECTOR_LEVELS   = 25           # number of flux levels sampled across the profile

# --- Cosmetics ---------------------------------------------------------------------------------
CMAP        = 'viridis'          # colour per epoch
FIGSIZE     = (15, 8)
SAVE_PREFIX = './review_'        # key 's' writes SAVE_PREFIX + <star>.png
VERBOSE     = True


# =============================================================================================
#  Readers -- adapt only if your FITS layout differs
# =============================================================================================

def READ_SPECTRUM(hdul):
    """
    Return (wave, flux) from an open HDUList.

    Tries, in order: the normalised column written by the continuum stage, the plain flux
    column, then a 1-D IRAF-style image with CRVAL1/CDELT1.
    """
    for ext in (-1, 1):
        try:
            data = hdul[ext].data
            wave = data['wave']
            for col in ('nflux', 'flux', 'Flux'):
                if col in data.names:
                    return np.asarray(wave, float), np.asarray(data[col], float)
        except Exception:
            pass
    try:                                        # older layout
        data = hdul[1].data
        return np.asarray(data['wavelength'], float), np.asarray(data['Flux'], float)
    except Exception:
        pass
    flux = np.asarray(hdul[0].data, float)      # 1-D image
    h = hdul[0].header
    wave = h['CRVAL1'] + np.arange(len(flux)) * h['CDELT1']
    return wave, flux


def READ_TEMPLATE(path):
    """Return (wave, flux) for the reference template. Default: AUTOKUR/ATLAS9 ASCII."""
    if path.lower().endswith(('.fits', '.fit')):
        with fits.open(path) as hdul:
            return READ_SPECTRUM(hdul)
    data = np.genfromtxt(path, comments='#')
    return data[:, 0], data[:, 1]


# =============================================================================================
#  Small utilities
# =============================================================================================

def boxcar(y, n):
    """Running mean over n pixels, edges preserved. Display only."""
    if n is None or n < 2:
        return y
    k = np.ones(int(n)) / float(int(n))
    pad = int(n) // 2
    padded = np.concatenate([np.repeat(y[0], pad), y, np.repeat(y[-1], pad)])
    return np.convolve(padded, k, mode='same')[pad:pad + len(y)]


def gauss_profile(x, center, sigma, depth, offset):
    """Absorption line: offset - depth * exp(-(x-center)^2 / 2 sigma^2)."""
    return offset - depth * np.exp(-(x - center) ** 2 / (2.0 * sigma ** 2))


def fit_line_profile(wave, flux, center, halfwidth):
    """
    Fit a Gaussian absorption profile in [center-halfwidth, center+halfwidth].

    Returns a dict with center, sigma, fwhm, depth, offset, equivalent width and the model
    curve, or None if the fit fails or the window is empty.
    """
    m = (wave > center - halfwidth) & (wave < center + halfwidth)
    if m.sum() < 6:
        return None
    x, y = wave[m], flux[m]
    p0 = [center, halfwidth / 4.0, max(np.nanmax(y) - np.nanmin(y), 1e-3), np.nanmax(y)]
    try:
        popt, pcov = curve_fit(gauss_profile, x, y, p0=p0, maxfev=20000)
    except Exception:
        return None
    c, s, d, o = popt
    err = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else np.full(4, np.nan)
    fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * abs(s)
    return {'center': c, 'center_err': err[0], 'sigma': abs(s), 'fwhm': fwhm,
            'depth': d, 'offset': o,
            'ew': d * abs(s) * np.sqrt(2 * np.pi) / o if o else np.nan,
            'rv_shift': (c - center) / center * 299792.458,
            'x': x, 'model': gauss_profile(x, *popt)}


def line_bisector(wave, flux, center, halfwidth, n_levels=25):
    """
    Bisector of an absorption line: the midpoint of the two wings at each flux level.

    A vertical bisector means a symmetric line; curvature indicates asymmetry, which for a
    resolved line traces velocity fields (granulation, expansion, blends).

    Returns (bisector_wavelengths, levels) or (None, None).
    """
    m = (wave > center - halfwidth) & (wave < center + halfwidth)
    if m.sum() < 8:
        return None, None
    x, y = wave[m], flux[m]
    core = np.nanargmin(y)
    if core < 2 or core > len(x) - 3:
        return None, None

    blue_x, blue_y = x[:core + 1], y[:core + 1]
    red_x,  red_y  = x[core:],     y[core:]
    order_b = np.argsort(blue_y)
    order_r = np.argsort(red_y)

    lo, hi = np.nanmin(y), min(np.nanmax(blue_y), np.nanmax(red_y))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return None, None
    levels = np.linspace(lo + 0.05 * (hi - lo), lo + 0.95 * (hi - lo), int(n_levels))

    try:
        wb = interp1d(blue_y[order_b], blue_x[order_b], bounds_error=False)(levels)
        wr = interp1d(red_y[order_r],  red_x[order_r],  bounds_error=False)(levels)
    except Exception:
        return None, None
    return (wb + wr) / 2.0, levels


# =============================================================================================
#  Building the star index
# =============================================================================================

def discover_obs(path, requested):
    if requested:
        return [o for o in requested if DISCARD_OBS.get(o) != '*']
    npy = os.path.join(path, 'Observing_blocks.npy')
    if os.path.exists(npy):
        keys = list(np.load(npy, allow_pickle=True).item().keys())
    else:
        keys = sorted(d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)))
    return [o for o in keys if DISCARD_OBS.get(o) != '*']


def star_name(filename):
    """The pipeline writes '<OBJECT>=<RA>=<DEC>=<DATE>.fits'."""
    base = os.path.basename(filename)
    if FILE_PREFIX and base.startswith(FILE_PREFIX):
        base = base[len(FILE_PREFIX):]
    return base.split('=')[0]


def coords_from_name(filename):
    """RA/DEC encoded in the filename, or (nan, nan)."""
    parts = os.path.basename(filename).split('=')
    try:
        return float(parts[1]), float(parts[2])
    except Exception:
        return np.nan, np.nan


def build_index():
    """
    Return (obs_ids, stars) where stars is a list of
        {'name', 'ra', 'dec', 'files': {obs: full_path}}
    """
    obs_ids = discover_obs(PATH, LIST)
    if not obs_ids:
        sys.exit(f'No observing blocks found under {PATH}.')

    per_obs = {}
    for obs in obs_ids:
        folder = PATH + obs + SUBFOLDER
        files = sorted(glob.glob(folder + FILE_PREFIX + '*.fits'))
        if not files:
            if VERBOSE:
                print(f'  ! {obs}: no FITS in {folder}')
            continue
        per_obs[obs] = files
    if not per_obs:
        sys.exit(f'No spectra found. Check PATH, SUBFOLDER ({SUBFOLDER}) and FILE_PREFIX.')

    obs_ids = list(per_obs)
    ref_obs = obs_ids[0]

    stars = []
    if MATCH_MODE == 'name':
        index = {}
        for obs, files in per_obs.items():
            for f in files:
                index.setdefault(star_name(f), {})[obs] = f
        for name, files in index.items():
            ra, dec = coords_from_name(files[next(iter(files))])
            stars.append({'name': name, 'ra': ra, 'dec': dec, 'files': files})

    elif MATCH_MODE == 'coords':
        base = per_obs[ref_obs]
        coords0 = np.array([coords_from_name(f) for f in base])
        first = SkyCoord(coords0[:, 0] * U.degree, coords0[:, 1] * U.degree)
        for f, (ra, dec) in zip(base, coords0):
            stars.append({'name': star_name(f), 'ra': ra, 'dec': dec,
                          'files': {ref_obs: f}})
        for obs in obs_ids[1:]:
            files = per_obs[obs]
            cc = np.array([coords_from_name(f) for f in files])
            other = SkyCoord(cc[:, 0] * U.degree, cc[:, 1] * U.degree)
            j, d2d, _ = first.match_to_catalog_sky(other)
            for i, (jj, dd) in enumerate(zip(j, d2d)):
                if dd < MATCH_RADIUS_ARCSEC * U.arcsec:
                    stars[i]['files'][obs] = files[jj]
    else:
        sys.exit("MATCH_MODE must be 'name' or 'coords'")

    # per-epoch and per-star rejections
    for st in stars:
        for obs, drop in DISCARD_OBS.items():
            if drop and drop != '*' and st['name'] in np.asarray(drop, dtype=str):
                st['files'].pop(obs, None)
    stars = [s for s in stars if s['name'] not in DISCARD_STARS and s['files']]

    if not stars:
        sys.exit('Every star was discarded -- check DISCARD_STARS / DISCARD_OBS.')

    # sorting
    keys = {'name':     lambda s: s['name'],
            'ra':       lambda s: (np.nan_to_num(s['ra'], nan=1e9)),
            'param_x':  lambda s: (PARAM_X or {}).get(s['name'], np.inf),
            'param_y':  lambda s: (PARAM_Y or {}).get(s['name'], np.inf),
            'n_epochs': lambda s: -len(s['files'])}
    if SORT_BY in keys:
        stars.sort(key=keys[SORT_BY])

    if VERBOSE:
        print(f'\n{len(stars)} stars over {len(obs_ids)} epochs: {", ".join(obs_ids)}')
        n = np.array([len(s['files']) for s in stars])
        print(f'epochs per star: min {n.min()}, median {int(np.median(n))}, max {n.max()}')
    return obs_ids, stars


# =============================================================================================
#  The viewer
# =============================================================================================

class Viewer:

    def __init__(self, obs_ids, stars):
        self.obs_ids = obs_ids
        self.stars = stars
        self.cache = {}                 # star index -> (wave, {obs: flux})
        self.current = 0
        self.show_bisector = BISECTOR_ON_START
        self.show_fit = FIT_ON_START
        self._warned_reference = False
        self._warned_trivial = False

        cmap = plt.get_cmap(CMAP)
        n = max(len(obs_ids) - 1, 1)
        self.colors = {o: cmap(i / n) for i, o in enumerate(obs_ids)}

        self.template = None
        if REFERENCE == 'template':
            if not TEMPLATE_PATH or not os.path.exists(TEMPLATE_PATH):
                print(f'  ! TEMPLATE_PATH not usable ({TEMPLATE_PATH!r}); '
                      f'falling back to the mean of the epochs.')
            else:
                self.template = READ_TEMPLATE(TEMPLATE_PATH)

        self._build_figure()
        self.update(0)

    # -----------------------------------------------------------------------------------------
    #  Data
    # -----------------------------------------------------------------------------------------
    def load(self, i):
        """Load and resample one star's epochs. Cached, so revisiting is instant."""
        if i in self.cache:
            return self.cache[i]

        star = self.stars[i]
        raw = {}
        for obs, path in star['files'].items():
            try:
                with fits.open(path) as hdul:
                    w, f = READ_SPECTRUM(hdul)
            except Exception as err:
                print(f'  ! could not read {os.path.basename(path)}: {err}')
                continue
            good = np.isfinite(w) & np.isfinite(f)
            w, f = w[good], f[good]
            if W_MIN is not None:
                m = w >= W_MIN; w, f = w[m], f[m]
            if W_MAX is not None:
                m = w <= W_MAX; w, f = w[m], f[m]
            if w.size < 10:
                continue
            if NORMALISE is not None:
                scale = NORMALISE(w, f)
                if scale:
                    f = f / scale
            raw[obs] = (w, f)

        if not raw:
            self.cache[i] = (None, {})
            return self.cache[i]

        # common grid: the overlap of all epochs
        lo = max(w.min() for w, _ in raw.values())
        hi = min(w.max() for w, _ in raw.values())
        first = raw[next(iter(raw))][0]
        if N_SAMPLE:
            grid = np.linspace(lo, hi, int(N_SAMPLE))
        else:
            grid = first[(first >= lo) & (first <= hi)]
        if grid.size < 10:                       # no overlap: keep the epochs on their own grid
            grid = first

        fluxes = {}
        for obs in self.obs_ids:
            if obs not in raw:
                continue
            w, f = raw[obs]
            fluxes[obs] = interp1d(w, f, bounds_error=False, fill_value=np.nan)(grid)

        self.cache[i] = (grid, fluxes)
        return self.cache[i]

    def reference_of(self, grid, fluxes):
        """Return the reference spectrum on `grid`, and a label for it."""
        if not fluxes:
            return None, ''
        stack = np.array(list(fluxes.values()))
        if REFERENCE == 'template' and self.template is not None:
            tw, tf = self.template
            ref = interp1d(tw, tf, bounds_error=False, fill_value=np.nan)(grid)
            if NORMALISE is not None:
                ok = np.isfinite(ref)
                scale = NORMALISE(grid[ok], ref[ok]) if ok.any() else None
                if scale:
                    ref = ref / scale
            return ref, 'template'
        if REFERENCE in fluxes:
            return fluxes[REFERENCE].copy(), f'reference {REFERENCE}'
        if REFERENCE == 'median':
            return np.nanmedian(stack, axis=0), 'median of epochs'
        if REFERENCE not in ('mean', 'median', 'template') and not self._warned_reference:
            print(f"  ! REFERENCE {REFERENCE!r} is not among the epochs loaded; "
                  f"using the mean instead.")
            self._warned_reference = True
        return np.nanmean(stack, axis=0), 'mean of epochs'

    # -----------------------------------------------------------------------------------------
    #  Figure
    # -----------------------------------------------------------------------------------------
    def _build_figure(self):
        self.fig = plt.figure(figsize=FIGSIZE)
        gs = gridspec.GridSpec(3, 4, height_ratios=[1.0, 0.45, 0.12],
                               width_ratios=[1, 0.06, 0.06, 0.55],
                               left=0.06, right=0.93, bottom=0.08, top=0.92,
                               wspace=0.05, hspace=0.05)

        self.ax_spec = self.fig.add_subplot(gs[0, 0])
        self.ax_res  = self.fig.add_subplot(gs[1, 0], sharex=self.ax_spec)
        self.ax_par  = self.fig.add_subplot(gs[0:2, 3])
        ax_slider    = self.fig.add_subplot(gs[2, 0])
        ax_prev      = self.fig.add_subplot(gs[2, 1])
        ax_next      = self.fig.add_subplot(gs[2, 2])

        self.ax_spec.set_ylabel('Flux')
        self.ax_res.set_ylabel('Flux / ref')
        self.ax_res.set_xlabel(r'Wavelength [$\AA$]')
        plt.setp(self.ax_spec.get_xticklabels(), visible=False)
        self.ax_res.axhline(1.0, color='k', lw=0.8, zorder=-5)

        # marks and highlighted regions on both panels
        for ax in (self.ax_spec, self.ax_res):
            for w, _lab in MARK_LINES:
                ax.axvline(w, color='k', lw=0.8, ls=':', zorder=-6)
            for band in HIGHLIGHT:
                w1, w2 = band[0], band[1]
                col = band[3] if len(band) > 3 else 'orange'
                ax.axvspan(w1, w2, color=col, alpha=0.15, zorder=-7)
        for w, lab in MARK_LINES:
            self.ax_spec.annotate(lab, (w, 1.005), xycoords=('data', 'axes fraction'),
                                  ha='center', fontsize=9)

        # one line per epoch, created empty and updated in place
        self.lines_spec = {o: self.ax_spec.plot([], [], lw=1.0, color=self.colors[o],
                                                label=o)[0] for o in self.obs_ids}
        self.lines_res  = {o: self.ax_res.plot([], [], lw=1.0,
                                               color=self.colors[o])[0] for o in self.obs_ids}
        self.line_ref,  = self.ax_spec.plot([], [], color='k', lw=1.6, zorder=5, label='ref')
        self.line_fit,  = self.ax_spec.plot([], [], color='fuchsia', lw=2.2, zorder=6)
        self.line_bis,  = self.ax_spec.plot([], [], color='cyan', lw=2.2, zorder=6)
        self.bands = []

        # Shown when the residual carries no information (one epoch compared with itself),
        # so a flat line at exactly 1 is never mistaken for "the epochs agree perfectly".
        self.trivial_note = self.ax_res.text(
            0.5, 0.5, '', transform=self.ax_res.transAxes, ha='center', va='center',
            fontsize=11, color='0.35',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.7', alpha=0.85))

        self.ax_spec.legend(loc='upper left', bbox_to_anchor=(1.005, 1.0),
                            fontsize=8, ncols=1 + len(self.obs_ids) // 18)

        # ---- parameter panel -----------------------------------------------------------------
        self._build_param_panel()

        # ---- widgets --------------------------------------------------------------------------
        self.slider = Slider(ax_slider, 'Star', 0, max(len(self.stars) - 1, 1e-9),
                             valinit=0, valstep=1)
        self.slider.on_changed(lambda v: self.update(int(v)))

        self.btn_prev = Button(ax_prev, '<')
        self.btn_next = Button(ax_next, '>')
        self.btn_prev.on_clicked(lambda e: self.step(-1))
        self.btn_next.on_clicked(lambda e: self.step(+1))

        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)

    def _build_param_panel(self):
        ax = self.ax_par
        names = [s['name'] for s in self.stars]

        if PARAM_X and PARAM_Y:
            x = np.array([PARAM_X.get(n, np.nan) for n in names], dtype=float)
            y = np.array([PARAM_Y.get(n, np.nan) for n in names], dtype=float)
            xlab, ylab = PARAM_X_LABEL, PARAM_Y_LABEL
            inv_x, inv_y = PARAM_X_INVERT, PARAM_Y_INVERT
            missing = np.isnan(x) | np.isnan(y)
            if missing.any() and VERBOSE:
                print(f'  ! {missing.sum()} star(s) have no PARAM_X/PARAM_Y value; '
                      f'they are not shown in the parameter panel.')
        else:
            x = np.array([s['ra'] for s in self.stars], dtype=float)
            y = np.array([s['dec'] for s in self.stars], dtype=float)
            xlab, ylab = 'RA [deg]', 'DEC [deg]'
            inv_x, inv_y = False, False
            if VERBOSE:
                print('  i PARAM_X / PARAM_Y not set -> the panel shows RA/DEC.')

        self.px, self.py = x, y

        if PARAM_BACKGROUND is not None:
            ax.scatter(PARAM_BACKGROUND[0], PARAM_BACKGROUND[1],
                       s=4, color='lightgray', zorder=-5)

        self.sc = ax.scatter(x, y, s=18, c='tab:blue', picker=True, pickradius=6)
        self.vline = ax.axvline(np.nan, color='red', lw=1)
        self.hline = ax.axhline(np.nan, color='red', lw=1)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        if inv_x and np.isfinite(x).any():
            ax.invert_xaxis()
        if inv_y and np.isfinite(y).any():
            ax.invert_yaxis()
        ax.yaxis.set_label_position('right')
        ax.yaxis.tick_right()
        ax.grid(ls=':', color='0.85')

    # -----------------------------------------------------------------------------------------
    #  Interaction
    # -----------------------------------------------------------------------------------------
    def step(self, delta):
        self.slider.set_val(int(np.clip(self.current + delta, 0, len(self.stars) - 1)))

    def on_pick(self, event):
        if event.artist is self.sc and len(event.ind):
            self.slider.set_val(int(event.ind[0]))

    def on_key(self, event):
        if event.key in ('right', 'up'):
            self.step(+1)
        elif event.key in ('left', 'down'):
            self.step(-1)
        elif event.key == 'f':
            self.show_fit = not self.show_fit
            self.update(self.current)
        elif event.key == 'b':
            self.show_bisector = not self.show_bisector
            self.update(self.current)
        elif event.key == 's':
            name = f'{SAVE_PREFIX}{self.stars[self.current]["name"]}.png'
            self.fig.savefig(name, dpi=150)
            print(f'  saved {name}')
        elif event.key == 'h':
            print(__doc__)

    # -----------------------------------------------------------------------------------------
    #  Draw
    # -----------------------------------------------------------------------------------------
    def update(self, i):
        self.current = int(np.clip(i, 0, len(self.stars) - 1))
        star = self.stars[self.current]
        grid, fluxes = self.load(self.current)

        for b in self.bands:
            b.remove()
        self.bands = []

        if grid is None:
            for o in self.obs_ids:
                self.lines_spec[o].set_visible(False)
                self.lines_res[o].set_visible(False)
            self.line_ref.set_visible(False)
            self.fig.suptitle(f'{star["name"]}  --  no readable spectrum', fontsize=13)
            self.fig.canvas.draw_idle()
            return

        ref, ref_label = self.reference_of(grid, fluxes)

        self.line_ref.set_data(grid, boxcar(ref, SMOOTH_N))
        self.line_ref.set_visible(True)
        self.line_ref.set_label(ref_label)

        residuals = []
        for o in self.obs_ids:
            ls, lr = self.lines_spec[o], self.lines_res[o]
            if o not in fluxes:
                ls.set_visible(False)
                lr.set_visible(False)
                continue
            f = fluxes[o]
            ls.set_data(grid, boxcar(f, SMOOTH_N))
            ls.set_visible(True)
            with np.errstate(divide='ignore', invalid='ignore'):
                r = f / ref
            residuals.append(r)
            lr.set_data(grid, boxcar(r, SMOOTH_N))
            lr.set_visible(True)

        # ---- is the residual informative at all? ----------------------------------------------
        # With a single epoch and REFERENCE='mean'/'median', the reference IS that spectrum, so
        # the residual is identically 1 and says nothing. Label it instead of letting a flat
        # line be misread as perfect agreement.
        trivial = (len(fluxes) == 1) and (REFERENCE in ('mean', 'median')
                                          or REFERENCE in fluxes)
        if trivial:
            self.trivial_note.set_text(
                'single epoch compared with itself: residual is 1 by construction\n'
                "set REFERENCE = 'template' (or another epoch) for a meaningful residual")
            if not self._warned_trivial:
                print("  i Only one epoch for this star and REFERENCE is that same spectrum, so "
                      "the residual\n    is identically 1. Set REFERENCE = 'template' (with "
                      "TEMPLATE_PATH) to compare against\n    a model, or pick another epoch.")
                self._warned_trivial = True
        else:
            self.trivial_note.set_text('')

        # ---- sigma bands ---------------------------------------------------------------------
        # Scatter of the residuals, propagated as 1/sqrt(reference): where the reference flux is
        # low (a line core) the same photon noise produces a larger relative deviation.
        if SHOW_SIGMA_BANDS and len(residuals) > 1:
            res = np.array(residuals)
            with np.errstate(invalid='ignore'):
                scaled = np.sqrt(np.abs(ref)) * (res - 1.0)
            sig = np.nanstd(scaled)
            if np.isfinite(sig) and sig > 0:
                denom = np.sqrt(np.abs(ref))
                denom[denom == 0] = np.nan
                for k, alpha in ((5, 0.12), (3, 0.16), (1, 0.22)):
                    self.bands.append(self.ax_res.fill_between(
                        grid, 1 - k * sig / denom, 1 + k * sig / denom,
                        color='gray', alpha=alpha, zorder=-10, linewidth=0))

        # ---- optional line diagnostics --------------------------------------------------------
        self.line_fit.set_visible(False)
        self.line_bis.set_visible(False)
        info = ''
        if self.show_fit:
            fit = fit_line_profile(grid, ref, FIT_CENTER, FIT_WINDOW)
            if fit:
                self.line_fit.set_data(fit['x'], fit['model'])
                self.line_fit.set_visible(True)
                info += (f"   fit: centre {fit['center']:.3f} A "
                         f"({fit['rv_shift']:+.1f} km/s), FWHM {fit['fwhm']:.3f} A, "
                         f"EW {fit['ew']:.3f} A")
        if self.show_bisector:
            bx, by = line_bisector(grid, ref, FIT_CENTER, FIT_WINDOW, BISECTOR_LEVELS)
            if bx is not None:
                self.line_bis.set_data(bx, by)
                self.line_bis.set_visible(True)

        # ---- crosshair on the parameter panel -------------------------------------------------
        if np.isfinite(self.px[self.current]):
            self.vline.set_xdata([self.px[self.current]])
        if np.isfinite(self.py[self.current]):
            self.hline.set_ydata([self.py[self.current]])

        title = f'{star["name"]}   [{self.current + 1}/{len(self.stars)}]   ' \
                f'{len(fluxes)} epoch(s): {", ".join(fluxes)}'
        self.fig.suptitle(title + info, fontsize=11)
        if VERBOSE:
            print(f'{self.current:4d}  {star["name"]:<24} {len(fluxes)} epochs{info}')

        self.ax_spec.legend(loc='upper left', bbox_to_anchor=(1.005, 1.0),
                            fontsize=8, ncols=1 + len(self.obs_ids) // 18)
        for ax in (self.ax_spec, self.ax_res):
            ax.relim()
            ax.autoscale_view()
        self.fig.canvas.draw_idle()


# =============================================================================================
#  Main
# =============================================================================================

def main():
    obs_ids, stars = build_index()
    viewer = Viewer(obs_ids, stars)
    print('\nControls: slider / [<] [>] / arrow keys = change star,  click the parameter panel '
          'to jump,\n          f = fit profile,  b = bisector,  s = save PNG,  h = help\n')
    plt.show()
    return viewer


if __name__ == '__main__':
    main()
