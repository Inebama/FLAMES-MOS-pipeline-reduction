# This file was rewritten by AI (Claude Opus, Anthropic) from I. Baeza's original script,
# to generalise it beyond the 47 Tuc / GIRAFFE science case it was written for.
# The intent, the statistics and the figures are the author's; the AI reorganised them into
# configurable, dataset-independent code. Please cross-check results before publishing.
"""
plotingRVS.py -- multi-epoch radial-velocity review
===================================================

Collects the `OUTPUT_RV.txt` catalogues written by Stage 3 for several Observing Blocks,
matches the stars across epochs, and produces:

  1. a mean radial velocity (and combined uncertainty) per star,
  2. a list of radial-velocity variable candidates,
  3. per-candidate time series + periodogram + phase-folded curve,
  4. the radial-velocity distribution of the whole sample (histogram + kernel densities).

Everything is configured in the CONFIGURATION block below. Nothing else needs editing.

Optional dependencies
---------------------
* GLS -- generalised Lomb-Scargle periodogram of Zechmeister & Kuerster (2009), A&A 496, 577
         (2009A&A...496..577Z).  Reference implementation:
             https://github.com/mzechmeister/GLS
         Python package directly:
             https://github.com/mzechmeister/GLS/tree/master/python
         Put the folder containing `gls.py` in GLS_PATH below, or `pip install gls`.
         If it is not found, the script falls back to `astropy.timeseries.LombScargle`,
         and if that is unavailable the period analysis is skipped (everything else still runs).
* scikit-learn -- only for the fixed-bandwidth Epanechnikov kernel density. Skipped if absent.
"""

# ---------------------------------------------------------------------------------------------
#  Packages
# ---------------------------------------------------------------------------------------------
import os
import sys
import glob
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as U

from scipy.optimize import curve_fit
from scipy.stats import norm

np.random.seed(170499)   # reproducible jitter in the rug plot


# =============================================================================================
#  CONFIGURATION -- this is the only block you need to edit
# =============================================================================================

# --- Where the pipeline output lives --------------------------------------------------------
PATH = './TESTING_FOLDER/TEST_OUTPUT/Test_example_try_to_reproduce_it/'   # must end with '/'

# Observing blocks to combine.
#   None  -> read them from PATH+'Observing_blocks.npy' (all of them, in order)
#   list  -> use exactly these, e.g. ['OB-0','OB-1','OB-4']
LIST = None

# Sub-folder of each OB that holds the spectra. Used ONLY to recover the epoch time (MJD)
# from a FITS header. Use '' to look directly in the OB folder. If no spectrum can be read,
# the script falls back to the OB order as the time axis and skips the period analysis.
SUBFOLDER = '/RESTFRAME/Continuum_sub/'

# Name of the RV catalogue inside each OB folder, and which of its columns to use.
# 'RV_corrected' is the barycentric (physical) radial velocity -- normally what you want.
RV_FILE   = 'OUTPUT_RV.txt'
RV_COLUMN = 'RV_corrected'      # or 'RV' (before removing the barycentric motion)

# --- How to match the same star across epochs -----------------------------------------------
#   'name'   -> by the OBJECT name in the RV catalogue (robust for this pipeline's own output)
#   'coords' -> by sky position, with the tolerance below (use for heterogeneous datasets)
MATCH_MODE          = 'name'
MATCH_RADIUS_ARCSEC = 1.75

# --- Rejecting bad data ---------------------------------------------------------------------
# Stars to drop entirely (by OBJECT name).
DISCARD_STARS = [
    # 'MED_RGB_370',
]

# Measurements to drop for one epoch only: {'OB name': [star, star, ...]}.
# Use the string '*' instead of a list to drop the whole epoch for every star.
# This is the place to record "that night was cloudy" / "this fibre was contaminated".
DISCARD_OBS = {
    # 'OB-2': ['MED_RGB_370', 'MED_RGB_402'],   # two stars had bad observations in this OB
    # 'OB-5': '*',                              # drop this epoch completely
}

# Also drop individual measurements whose formal RV error exceeds this value [km/s].
# None disables the cut.
MAX_RV_ERR = None

# --- Variability selection -------------------------------------------------------------------
# Threshold on the combined uncertainty (measurement + epoch-to-epoch scatter) above which a
# star is flagged as a radial-velocity variable candidate, in km/s.
#   None -> automatic: median(sigma) + N_SIGMA_AUTO * robust scatter of sigma
VAR_THRESHOLD = None
N_SIGMA_AUTO  = 5.0
MIN_EPOCHS    = 2       # a star needs at least this many surviving epochs (2 is the minimum
                        # for any scatter to exist at all; raise it for a stricter sample)

MAX_CANDIDATES_PLOTTED = 12   # keeps the time-series figure readable

# --- Period search ---------------------------------------------------------------------------
GLS_PATH     = ''            # folder containing gls.py; '' if installed or unavailable
PERIOD_RANGE = (2.0, 1000.0) # (Pbeg, Pend) in days
OVERSAMPLING = 2000          # 'ofac' of the periodogram; lower it if it feels slow

# --- Output ----------------------------------------------------------------------------------
SAVE_SUMMARY = 'RV_SUMMARY.txt'   # None to skip
SAVE_FIGURES = None               # e.g. './rv_figs/' to save instead of showing; None -> show
VERBOSE      = True


# =============================================================================================
#  Helper functions
# =============================================================================================

def gaussian(x, x0, sigma, a=1.0):
    """Simple Gaussian, used to fit the RV histogram."""
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def Folding(t, width, ref=None):
    """
    Fold a time array on a period.

    t     : times (days)
    width : folding period (days)
    ref   : reference epoch; the minimum of t if not given

    Returns phases in [0, 1).
    """
    t = np.asarray(t, dtype=float)
    if ref is None:
        ref = np.nanmin(t)
    phase = ((t - ref) % width) / width
    phase[phase < 0] += 1
    return phase


def reorderLegend(ax=None, order=None):
    """Reorder legend entries to follow `order` (entries not listed go last)."""
    if ax is None:
        ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None, None
    if order is not None:
        rank = {lab: i for i, lab in enumerate(order)}
        labels, handles = zip(*sorted(zip(labels, handles),
                                      key=lambda t: rank.get(t[0], np.inf)))
    ax.legend(handles, labels)
    return handles, labels


def variable_bandwidth_kde(X, bandwidths, query_points, kernel='gaussian'):
    """
    Kernel density estimate in which **each point carries its own bandwidth**.

    Using the individual RV uncertainties as bandwidths means precise measurements
    contribute narrow kernels and poor ones broad kernels -- an honest way of showing a
    velocity distribution built from heterogeneous errors.

    X, bandwidths : same length (data and their 1-sigma errors)
    query_points  : where to evaluate the density
    kernel        : 'gaussian' or 'tophat'
    """
    X = np.asarray(X).ravel()
    bandwidths = np.asarray(bandwidths).ravel()
    query_points = np.asarray(query_points).ravel()

    good = np.isfinite(X) & np.isfinite(bandwidths) & (bandwidths > 0)
    X, bandwidths = X[good], bandwidths[good]
    if X.size == 0:
        return np.zeros_like(query_points)

    densities = np.zeros(query_points.shape)
    for i, qp in enumerate(query_points):
        if kernel == 'gaussian':
            k = norm.pdf((qp - X) / bandwidths) / bandwidths
        elif kernel == 'tophat':
            k = np.where(np.abs(qp - X) <= bandwidths, 1.0 / (2 * bandwidths), 0.0)
        else:
            raise ValueError("Unsupported kernel type. Use 'gaussian' or 'tophat'.")
        densities[i] = np.sum(k) / X.size
    return densities


def robust_scatter(x):
    """Median absolute deviation scaled to a Gaussian sigma. Immune to a few outliers."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    return 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))


# =============================================================================================
#  Reading the pipeline output
# =============================================================================================

def discover_obs(path, requested):
    """Return the list of OB folder names to use."""
    if requested:
        return list(requested)
    npy = os.path.join(path, 'Observing_blocks.npy')
    if os.path.exists(npy):
        blocks = np.load(npy, allow_pickle=True).item()
        return list(blocks.keys())
    # last resort: any subfolder that looks like an OB and has an RV catalogue
    found = sorted(d for d in os.listdir(path)
                   if os.path.isdir(os.path.join(path, d))
                   and os.path.exists(os.path.join(path, d, RV_FILE)))
    if not found:
        sys.exit(f'No observing blocks found under {path}. Set LIST explicitly.')
    return found


def read_rv_table(filename):
    """
    Read one OUTPUT_RV.txt.

    The column names are taken from the last commented line of the header, so the reader keeps
    working if the column order ever changes. Falls back to the documented order if no usable
    header is found.

    Returns a dict of arrays keyed by column name.
    """
    default = ['Object', 'RA', 'DEC', 'RV', 'RV_ERR',
               'Barycentric_RV', 'RV_corrected', 'template_used_in_CCF']

    names = None
    with open(filename) as f:
        for line in f:
            if not line.startswith('#'):
                break
            fields = line.lstrip('#').strip().split('\t')
            if len(fields) > 3:                 # a real column header, not a comment sentence
                names = [x.strip() for x in fields]
    if names is None or len(names) < 5:
        names = default

    raw = np.genfromtxt(filename, dtype=str, comments='#', delimiter='\t')
    if raw.ndim == 1:                            # a single row
        raw = raw.reshape(1, -1)

    table = {}
    for i, name in enumerate(names):
        if i >= raw.shape[1]:
            break
        col = raw[:, i]
        try:
            table[name] = col.astype(float)
        except ValueError:
            table[name] = col
    return table


def epoch_time(path, obs, subfolder):
    """
    Recover the MJD of an epoch from the header of any spectrum in it.
    Returns np.nan if nothing can be read.
    """
    for pattern in (path + obs + subfolder + '*.fits', path + obs + '/*.fits'):
        for f in sorted(glob.glob(pattern)):
            try:
                with fits.open(f) as hdul:
                    return float(hdul[0].header['MJD'])
            except Exception:
                continue
    return np.nan


def build_catalogue():
    """
    Read every OB, match the stars, and return

      names   : (N_stars,)            object names
      coords  : SkyCoord              positions from the first epoch a star appears in
      mjd     : (N_epochs,)           epoch times (or indices if no MJD was readable)
      obs_ids : (N_epochs,)           OB names
      rv      : (N_epochs, N_stars)   radial velocities, NaN where missing/discarded
      rve     : (N_epochs, N_stars)   their uncertainties
    """
    obs_ids = discover_obs(PATH, LIST)

    tables, kept_obs = [], []
    for obs in obs_ids:
        fname = PATH + obs + '/' + RV_FILE
        if not os.path.exists(fname):
            if VERBOSE:
                print(f'  ! {obs}: no {RV_FILE}, skipped')
            continue
        if DISCARD_OBS.get(obs) == '*':
            if VERBOSE:
                print(f'  - {obs}: discarded by configuration')
            continue
        tables.append(read_rv_table(fname))
        kept_obs.append(obs)

    if not tables:
        sys.exit('No RV catalogues could be read. Check PATH, LIST and RV_FILE.')
    if RV_COLUMN not in tables[0]:
        sys.exit(f'Column "{RV_COLUMN}" not in {RV_FILE}. Available: {list(tables[0])}')

    # --- the first epoch defines the star list ------------------------------------------------
    names = np.array([str(n) for n in tables[0]['Object']])
    ra, dec = tables[0]['RA'], tables[0]['DEC']

    keep = np.array([n not in DISCARD_STARS for n in names])
    names, ra, dec = names[keep], ra[keep], dec[keep]
    first = SkyCoord(ra * U.degree, dec * U.degree)

    n_star = len(names)
    rv  = np.full((len(tables), n_star), np.nan)
    rve = np.full((len(tables), n_star), np.nan)

    for e, (obs, tab) in enumerate(zip(kept_obs, tables)):
        this_names = np.array([str(n) for n in tab['Object']])
        this_rv    = np.asarray(tab[RV_COLUMN], dtype=float)
        this_err   = np.asarray(tab.get('RV_ERR', np.full(len(this_names), np.nan)), dtype=float)

        if MATCH_MODE == 'name':
            lookup = {n: i for i, n in enumerate(this_names)}
            idx = np.array([lookup.get(n, -1) for n in names])
        elif MATCH_MODE == 'coords':
            other = SkyCoord(np.asarray(tab['RA']) * U.degree,
                             np.asarray(tab['DEC']) * U.degree)
            j, d2d, _ = first.match_to_catalog_sky(other)
            idx = np.where(d2d < MATCH_RADIUS_ARCSEC * U.arcsec, j, -1)
        else:
            sys.exit("MATCH_MODE must be 'name' or 'coords'")

        found = idx >= 0
        rv[e, found]  = this_rv[idx[found]]
        rve[e, found] = this_err[idx[found]]

        # per-epoch rejections requested by the user
        drop = DISCARD_OBS.get(obs, [])
        if drop and drop != '*':
            bad = np.isin(names, np.asarray(drop, dtype=str))
            rv[e, bad], rve[e, bad] = np.nan, np.nan
            if VERBOSE and bad.any():
                print(f'  - {obs}: dropped {bad.sum()} star(s) by configuration')

        if MAX_RV_ERR is not None:
            bad = rve[e] > MAX_RV_ERR
            rv[e, bad], rve[e, bad] = np.nan, np.nan

        if VERBOSE:
            print(f'  + {obs}: {np.isfinite(rv[e]).sum():4d}/{n_star} stars matched')

    mjd = np.array([epoch_time(PATH, obs, SUBFOLDER) for obs in kept_obs])
    if not np.isfinite(mjd).all():
        if VERBOSE:
            print('\n  ! Could not read an MJD for every epoch; using the OB order as the '
                  'time axis and skipping the period search.')
        mjd = np.arange(len(kept_obs), dtype=float)

    return names, first, mjd, np.array(kept_obs), rv, rve


# =============================================================================================
#  Period analysis (optional)
# =============================================================================================

def get_periodogram_backend():
    """
    Return ('gls', Gls) | ('astropy', LombScargle) | (None, None).

    GLS -- Zechmeister & Kuerster (2009), A&A 496, 577; https://github.com/mzechmeister/GLS
    """
    if GLS_PATH:
        sys.path.append(GLS_PATH)
    try:
        from gls import Gls
        return 'gls', Gls
    except Exception:
        pass
    try:
        from astropy.timeseries import LombScargle
        if VERBOSE:
            print('  ! gls not found -- falling back to astropy.timeseries.LombScargle.')
            print('    For the original generalised Lomb-Scargle see '
                  'https://github.com/mzechmeister/GLS/tree/master/python')
        return 'astropy', LombScargle
    except Exception:
        if VERBOSE:
            print('  ! No periodogram backend available; period analysis skipped.')
        return None, None


def best_period(t, y, dy, backend, cls):
    """
    Fit a sinusoid by periodogram.

    Returns (period, T0, amplitude) or (nan, nan, nan) if it cannot be done.
    T0 is defined so that the sinusoid is sin(2*pi*(t-T0)/P).
    """
    good = np.isfinite(t) & np.isfinite(y)
    if backend is None or good.sum() < 4 or np.ptp(t[good]) <= 0:
        return np.nan, np.nan, np.nan

    t, y = t[good], y[good]
    dy = dy[good]
    dy = np.where(np.isfinite(dy) & (dy > 0), dy, np.nanmedian(dy[np.isfinite(dy)]) or 1.0)

    Pbeg, Pend = PERIOD_RANGE
    Pend = min(Pend, 3 * np.ptp(t)) if np.ptp(t) > 0 else Pend
    if Pend <= Pbeg:
        return np.nan, np.nan, np.nan

    try:
        if backend == 'gls':
            g = cls((t, y, dy), Pbeg=Pbeg, Pend=Pend, ofac=OVERSAMPLING)
            if VERBOSE:
                g.info()
            return 1.0 / g.hpstat['f'], g.hpstat['T0'], g.hpstat['amp']
        else:
            ls = cls(t, y, dy)
            freq, power = ls.autopower(minimum_frequency=1.0 / Pend,
                                       maximum_frequency=1.0 / Pbeg,
                                       samples_per_peak=10)
            f = freq[np.argmax(power)]
            # amplitude and phase from the fitted sinusoid
            theta = ls.model_parameters(f)
            amp = np.hypot(theta[1], theta[2]) if len(theta) >= 3 else np.nanstd(y) * np.sqrt(2)
            phase = np.arctan2(theta[2], theta[1]) if len(theta) >= 3 else 0.0
            T0 = -phase / (2 * np.pi * f)
            return 1.0 / f, T0, amp
    except Exception as err:
        if VERBOSE:
            print(f'    periodogram failed: {err}')
        return np.nan, np.nan, np.nan


# =============================================================================================
#  Figures
# =============================================================================================

def show_or_save(fig, name):
    if SAVE_FIGURES:
        os.makedirs(SAVE_FIGURES, exist_ok=True)
        fig.savefig(os.path.join(SAVE_FIGURES, name), dpi=150)
        plt.close(fig)
    # the actual plt.show() happens once, at the end


def figure_variables(names, mjd, rv, rve, cand, backend, cls):
    """Per-candidate RV range, time series and phase-folded curve."""
    if len(cand) == 0:
        print('\nNo radial-velocity variable candidates above the threshold.')
        return

    fig = plt.figure(figsize=(15, 1.2 * len(cand) + 2))
    gs = gridspec.GridSpec(len(cand), 3,
                           height_ratios=[1] * len(cand), width_ratios=[1, 0.7, 0.7],
                           left=0.14, right=0.88, bottom=0.10, top=0.90,
                           wspace=0.05, hspace=0.02)

    ax_range = fig.add_subplot(gs[:, 0])
    ax_range.set_title('RV range of candidates')

    ax_ts = ax_ph = None
    for row, i in enumerate(cand):
        y  = rv[:, i]
        dy = rve[:, i]
        centre = np.nanmedian(y)

        ax_ts = fig.add_subplot(gs[row, 1], sharex=None if row == 0 else ax_ts)
        ax_ph = fig.add_subplot(gs[row, 2], sharex=None if row == 0 else ax_ph)
        if row == 0:
            ax_ts.set_title('Time series')
            ax_ph.set_title('Phase fold')

        ax_range.errorbar(y, np.full(len(y), len(cand) - 1 - row), xerr=dy,
                          elinewidth=2, linewidth=0, marker='o', markersize=6,
                          ecolor='gray', label=names[i])

        if VERBOSE:
            print(f'\n--- {names[i]} ---')
        P, T0, amp = best_period(mjd, y - centre, dy, backend, cls)

        ax_ts.errorbar(mjd, y - centre, yerr=dy, elinewidth=1, linewidth=0,
                       marker='o', color='k', markersize=5)
        if np.isfinite(P):
            x = np.linspace(np.nanmin(mjd), np.nanmax(mjd), 1000)
            ax_ts.plot(x, np.sin(2 * np.pi * (x - T0) / P) * amp)

            ph = Folding(mjd, P, T0)
            xx = np.linspace(0, 1, 200)
            ax_ph.plot(xx, np.sin(2 * np.pi * xx) * amp)
            ax_ph.errorbar(ph, y - centre, yerr=dy, elinewidth=1, linewidth=0,
                           marker='o', color='k', markersize=5)
            ax_ph.legend(title=f'Period\n{P:.2f} d', loc='upper left', bbox_to_anchor=(1, 1))
        else:
            ax_ph.text(0.5, 0.5, 'no period', ha='center', va='center',
                       transform=ax_ph.transAxes, color='gray')

        for ax in (ax_ts, ax_ph):
            ax.yaxis.set_visible(False)
            ax.xaxis.set_visible(row == len(cand) - 1)

    ax_range.set_yticks(np.arange(len(cand))[::-1], labels=names[cand], fontsize=11)
    ax_range.set_xlabel('RV [km/s]')
    ax_range.grid(ls=':', color='gray')
    if ax_ts is not None:
        ax_ts.set_xlabel('MJD')
        ax_ph.set_xlabel('Phase')

    show_or_save(fig, 'RV_variables.png')


def figure_distribution(rv_mean, rv_err):
    """Histogram of the mean RVs plus three density estimates and a rug plot."""
    good = np.isfinite(rv_mean)
    x, dx = rv_mean[good], rv_err[good]
    if x.size < 3:
        print('\nToo few stars with a valid mean RV to build the distribution figure.')
        return

    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.2], left=0.09, right=0.97,
                           bottom=0.10, top=0.92, hspace=0.0)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax2.yaxis.set_visible(False)

    nbins = max(10, min(40, int(np.sqrt(x.size) * 2)))
    hist = ax1.hist(x, bins=nbins, color='gray', edgecolor='None', density=True)
    width = np.mean(np.diff(hist[1]))
    hist[-1].set_label(f'Histogram RV (bin size = {width:.2f})')

    centres = (hist[1][1:] + hist[1][:-1]) / 2
    span = np.nanmax(x) - np.nanmin(x)
    grid = np.linspace(np.nanmin(x) - 0.1 * span, np.nanmax(x) + 0.1 * span, 400)

    try:
        popt, _ = curve_fit(gaussian, centres, hist[0],
                            p0=[np.nanmedian(x), max(np.nanstd(x), 1e-3), np.nanmax(hist[0])])
        ax1.plot(grid, gaussian(grid, *popt), color='k',
                 label='Gaussian Fit to Histogram Values')
        print(f'\nGaussian fit to the RV histogram: '
              f'centre = {popt[0]:.3f}, sigma = {abs(popt[1]):.3f} km/s')
    except Exception as err:
        print(f'\nGaussian fit to the histogram failed: {err}')

    ax1.plot(grid, variable_bandwidth_kde(x, dx, grid, kernel='gaussian'),
             label='Gaussian Kernel Density (bandwidth = individual errors)')

    try:
        from sklearn.neighbors import KernelDensity
        bw = 5 * np.nanmean(dx) if np.isfinite(np.nanmean(dx)) else np.nanstd(x) / 3
        kde = KernelDensity(kernel='epanechnikov', bandwidth=max(bw, 1e-3))
        kde.fit(x.reshape(-1, 1))
        ax1.plot(grid, np.exp(kde.score_samples(grid.reshape(-1, 1))),
                 label='Epanechnikov Kernel Density (bandwidth = 5 x mean error)')
    except Exception:
        pass    # scikit-learn is optional

    ax1.axvline(np.nanmedian(x), color='g', label='Median RV')
    ax1.axvline(np.nanmean(x),   color='r', label='Mean RV')
    ax1.errorbar([], [], xerr=[], elinewidth=2, linewidth=0, marker='|', color='k',
                 markersize=15, label='Data Measured')

    ax1.set_ylabel('Density', fontsize=12)
    ax2.set_xlabel('RV [km/s]', fontsize=12)
    ax1.legend(loc='upper right', fontsize=8)
    reorderLegend(ax1, [f'Histogram RV (bin size = {width:.2f})',
                        'Gaussian Fit to Histogram Values',
                        'Gaussian Kernel Density (bandwidth = individual errors)',
                        'Epanechnikov Kernel Density (bandwidth = 5 x mean error)',
                        'Median RV', 'Mean RV', 'Data Measured'])

    ax2.errorbar(x, np.random.normal(0, 0.3, x.size), xerr=dx, elinewidth=2, linewidth=0,
                 marker='|', color='k', markersize=15)
    ax2.set_ylim(-1, 1)

    show_or_save(fig, 'RV_distribution.png')


# =============================================================================================
#  Main
# =============================================================================================

def main():
    print('\nReading RV catalogues ...')
    names, coords, mjd, obs_ids, rv, rve = build_catalogue()

    n_epoch = np.isfinite(rv).sum(axis=0)

    # --- combine epochs -----------------------------------------------------------------------
    # np.nanmean/nanstd legitimately warn on all-NaN columns (a star missing from every epoch);
    # those columns are removed by the MIN_EPOCHS cut two lines below.
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        rv_mean = np.nanmean(rv, axis=0)
        # measurement error and epoch-to-epoch scatter added in quadrature:
        # a genuinely variable star therefore ends up with a large sigma, which is what the
        # variability selection below keys on.
        rv_err = np.sqrt(np.nanmean(rve, axis=0) ** 2 + np.nanstd(rv, axis=0) ** 2)

    enough = n_epoch >= MIN_EPOCHS
    rv_mean[~enough] = np.nan
    rv_err[~enough] = np.nan

    print(f'\n{len(names)} stars, {len(obs_ids)} epochs '
          f'({enough.sum()} stars with >= {MIN_EPOCHS} epochs)')
    if not enough.any():
        sys.exit(f'\nNo star has {MIN_EPOCHS} or more surviving epochs -- nothing to analyse.\n'
                 f'Lower MIN_EPOCHS, add more OBs to LIST, or relax DISCARD_OBS / MAX_RV_ERR.')

    print('Stats over stars with enough epochs')
    print('   mean       median     std        median(sigma)')
    print(f'   {np.nanmean(rv_mean):<10.3f} {np.nanmedian(rv_mean):<10.3f} '
          f'{np.nanstd(rv_mean):<10.3f} {np.nanmedian(rv_err):<10.3f}')

    # --- variability selection ----------------------------------------------------------------
    threshold = VAR_THRESHOLD
    if threshold is None:
        threshold = np.nanmedian(rv_err) + N_SIGMA_AUTO * robust_scatter(rv_err)
        print(f'\nAutomatic variability threshold: {threshold:.3f} km/s '
              f'(median sigma + {N_SIGMA_AUTO:g} x robust scatter)')
    else:
        print(f'\nVariability threshold: {threshold:.3f} km/s')

    cand = np.where(enough & (rv_err >= threshold))[0]
    cand = cand[np.argsort(rv_err[cand])[::-1]]          # worst offenders first
    if len(cand) > MAX_CANDIDATES_PLOTTED:
        print(f'  {len(cand)} candidates found; plotting the {MAX_CANDIDATES_PLOTTED} '
              f'with the largest scatter.')
        cand = cand[:MAX_CANDIDATES_PLOTTED]

    if len(cand):
        print(f'\n{"idx":>5} {"Object":<20} {"RA":>11} {"DEC":>11} {"sigma":>8} {"N_ep":>5}')
        for i in cand:
            print(f'{i:>5} {names[i]:<20} {coords.ra.value[i]:>11.6f} '
                  f'{coords.dec.value[i]:>11.6f} {rv_err[i]:>8.3f} {n_epoch[i]:>5}')

    # --- figures -------------------------------------------------------------------------------
    backend, cls = get_periodogram_backend()
    figure_variables(names, mjd, rv, rve, cand, backend, cls)
    figure_distribution(rv_mean, rv_err)

    # --- summary table -------------------------------------------------------------------------
    if SAVE_SUMMARY:
        out = os.path.join(PATH, SAVE_SUMMARY)
        with open(out, 'w') as f:
            f.write(f'# RV summary from {RV_COLUMN} over {len(obs_ids)} epochs: '
                    f'{", ".join(obs_ids)}\n')
            f.write(f'# variability threshold = {threshold:.4f} km/s\n')
            f.write('#OBJ\tRA\tDEC\tRV\tRV_ERR\tN_EPOCHS\tVARIABLE\n')
            for i in range(len(names)):
                f.write(f'{names[i]}\t{coords.ra.value[i]}\t{coords.dec.value[i]}\t'
                        f'{rv_mean[i]}\t{rv_err[i]}\t{n_epoch[i]}\t{int(i in cand)}\n')
        print(f'\nSummary written to {out}')

    if not SAVE_FIGURES:
        plt.show()


if __name__ == '__main__':
    main()
