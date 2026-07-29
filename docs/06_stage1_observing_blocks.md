# 06 — Stage 1: Observing-Block definition

[← back to index](README.md)

**Module:** `List_files_related.py` — function `OBSERVING_BLOCK`
**Driver:** the `if LIST:` branch of `MAIN.py`
**Products:** `Observing_blocks.npy`, `Observing_blocks.txt`

---

## 1. What an Observing Block is here

An **Observing Block (OB)** in this package is *a set of exposures that you are willing to
co-add into a single epoch*. It is a science decision, not an instrument one: an OB should
group exposures over a timescale during which you expect **no significant intrinsic variation**
of the spectral features you are studying.

Examples:

| Science case | Sensible threshold |
|---|---|
| Long-period variables, mass-loss episodes in red giants | 24–48 h (whole nights) |
| Close binaries, short-period pulsators | 0.25–1 h |
| Testing exposure-to-exposure repeatability | < 0.1 h (isolates individual exposures) |

The trade-off is the usual one: longer blocks give higher S/N and fewer epochs, shorter blocks
give finer time sampling and noisier spectra. Because Stage 1 is cheap and non-destructive,
the recommended practice is to try several thresholds and compare (this is exactly what the
author did — see the commented `SUPER_LIST` dictionaries in `cheating.py`, which map a fine OB
grouping onto coarser ones).

## 2. The grouping algorithm

Driven by `Range_in_hours_to_be_consider_same_Observing_Block` (call it `T`, in hours).

```
1.  Sort all exposures by MJD.
2.  seed ← the earliest exposure not yet assigned to any OB
3.  candidates ← every unassigned exposure with |MJD − MJD(seed)| < T/24   (days)
4.  Apply the positional consistency test (§3) to the candidates, in time order.
5.  The surviving candidates become OB-j; mark them assigned.
6.  If any exposure remains unassigned → j ← j+1, go to 2.
```

Two properties worth being explicit about:

* **It is a window around the seed, not a chained/gap clustering.** Exposures 25 h apart never
  share an OB when `T = 24`, even if a third exposure sits between them — that third exposure
  joins the first block, and the last one seeds a new block. This makes the result
  deterministic and easy to reason about, at the cost of some sensitivity to where the seed
  falls.
* **Assignment is greedy and chronological.** Blocks are numbered `OB-0`, `OB-1`, … in the
  order of their seeds, so `OB-0` is always the earliest epoch.

### Implementation note

The loop is written as `for j,i in enumerate(AUX)` where `AUX` starts as `[time[0]]` and grows
inside the loop (`AUX.append(time[temp][0])` whenever unassigned exposures remain). Iterating
over a list that is being appended to is what makes the "repeat until nothing is left" loop
work without an explicit `while`.

## 3. The positional consistency test

Grouping by time alone is not enough: two exposures taken the same night may point at a
*different plate configuration*. Stage 2 will co-add fibre-by-fibre, so it needs to be certain
that fibre *k* observes the same star in every exposure of the block.

For each consecutive pair of candidate exposures `(k, k+1)`:

```python
current = SkyCoord(RA[k]  [not_sky[k]],   DEC[k]  [not_sky[k]])
nxt     = SkyCoord(RA[k+1][not_sky[k+1]], DEC[k+1][not_sky[k+1]])
idx, d2d, _ = current.match_to_catalog_sky(nxt)
sep = d2d > 0.5*u.arcsec
keep = (sum(sep) == 0) and (len(current) == len(nxt))
```

That is: **every science fibre of exposure *k* must have a counterpart in exposure *k+1*
within 0.5″, and both exposures must contain the same number of science fibres.** Sky fibres
are excluded from the test — their number and positions are allowed to differ.

If the test fails, that exposure (and, through the cumulative `mask[current_index[k+1]] =
(ras & decs) & mask[current_index[k]]`, everything after it in the block) is removed from the
current OB and stays unassigned, so it will seed or join a later block.

**Consequences**

* A field observed with two different plate configurations on the same night produces two
  separate OBs, which is what you want.
* A single exposure with one dead/missing fibre will break the length test and be split off.
  If that is not what you want, either drop the exposure or relax the test in
  `List_files_related.py` (the 0.5″ threshold and the length check are both hard-coded there).

> The same test is repeated, independently and against the **first** exposure of the block,
> inside `MERGING_MULTIPLE_SPECTRA`; if it fails there, Stage 2 raises
> `Exception('The Values of RA DEC does not match between files …')`. See
> [14 — Troubleshooting](14_troubleshooting.md).

## 4. Function signature

```python
OBSERVING_BLOCK(path,                              # output folder (must end with '/')
                time,                              # iterable of MJD, one per exposure
                epochs,                            # iterable of file/folder names
                RA,                                # list of RA arrays, one per exposure
                DEC,                               # list of DEC arrays, one per exposure
                not_sky        = None,             # list of boolean arrays; None ⇒ no sky fibres
                OBS_BLOCK_TIME = 12,               # threshold in hours
                shapes         = None,             # optional, printed in the table
                Fname          = 'Observing_blocks.npy',
                SAVE_TXT       = True) -> None
```

It returns nothing; it writes the two files and prints one table per block.

`MAIN.py` calls it as:

```python
OBSERVING_BLOCK(OBS_BLOCK_TIME = Range_in_hours_to_be_consider_same_Observing_Block,
                path    = path_output,
                epochs  = Files,               # names relative to path_input
                time    = DATA['MJD'],
                RA      = DATA['RAS'],
                DEC     = DATA['DECS'],
                not_sky = DATA.get('NOT_SKY'),
                shapes  = DATA.get('SHAPE'),
                SAVE_TXT= SAVE_TXT)
```

so it is perfectly usable stand-alone on any dataset for which you can produce those five
arrays.

## 5. Reading the output

```python
import numpy as np
OB = np.load('/path/to/out/Observing_blocks.npy', allow_pickle=True).item()

for name, files in OB.items():
    print(name, len(files), files[0])
```

`Observing_blocks.txt` contains the same information plus the time offsets from the seed and
the data shapes — useful for spotting an exposure that sneaked in at the edge of the window.

## 6. Checklist before moving to Stage 2

* [ ] The number of OBs matches your expectation from the observing log.
* [ ] Within each block, `Diff_in_Hour` is compatible with the variability timescale you care about.
* [ ] `Shapes_of_Data` is identical for every exposure of a block.
* [ ] No block contains a single exposure by accident (unless that is intended) — a lone
      exposure usually signals a failed positional test.
* [ ] Total number of exposures across all blocks equals the number of input files.
