# SpexRock

**Asteroid reflectance spectra from NASA IRTF SpeX, end to end.**

SpexRock reduces raw SpeX frames — prism, SXD, or LXD mode, pre- or
post-upgrade instrument — to a normalized asteroid reflectance spectrum in
one command. It is a thin, batch-friendly pipeline built on
[pyspextool](https://github.com/pyspextool/pyspextool), the official Python
implementation of the Spextool reduction protocols (Cushing, Vacca & Rayner
2004), so every instrumental step follows the community-standard procedures
maintained by the instrument team. SpexRock adds the asteroid-specific
layer that Spextool leaves to the user:

* **Solar-analog telluric correction** — the object/analog ratio
  (pyspextool's `reflectance` mode, the Python descendant of
  `xtellcor_basic`) that removes the atmosphere and the solar spectrum in
  one division, with pyspextool's sub-pixel wavelength-shift optimization;
* **Normalized reflectance products** — unity at 1.20 µm (Bus-DeMeo
  convention) or any wavelength you choose, written as FITS + ASCII + CSV
  with full provenance (frame numbers, package versions, the exact config);
* **NEATM thermal-excess removal** (optional) — for LXD data beyond
  ~2.5 µm, the Harris (1998) model evaluated from your observing
  circumstances, with the beaming parameter either fixed or fit to the
  spectrum's own thermal tail;
* **One JSON config per reduction** — a `ReduceConfig` dataclass replaces
  the notebook top cell; the config used is snapshotted next to the
  products so any spectrum can be re-derived exactly.

## How it works

The pipeline runs the canonical Spextool sequence via pyspextool, then the
asteroid-specific steps:

1. **Calibrations** — master flat with order-edge tracing
   (`make_flat`), wavelength solution from argon arcs — supplemented by
   sky emission lines in LXD mode, where the argon lamp runs out of lines
   beyond ~4 µm (`make_wavecal`, `sky_files`).
2. **Extraction** — nonlinearity and bias corrections, A−B pair
   subtraction, flat fielding, spatial-profile-based aperture location and
   tracing, optimal extraction with bad-pixel rejection — identical
   settings for asteroid and analog, as required for the ratio to cancel
   instrument signatures.
3. **Combination** — robust weighted mean of the extracted spectra.
4. **Telluric correction / reflectance** — asteroid ÷ solar analog with
   `correction_type='reflectance'`; the airmass difference is checked and
   loudly warned about beyond a configurable threshold.
5. **Order merging** — for cross-dispersed modes, orders are scaled and
   merged on the telluric-corrected spectrum (prism skips this).
6. **Normalization** — to unity at 1.20 µm (prism/SXD default) or 3.55 µm
   (LXD default), configurable.
7. **NEATM thermal-excess removal** (optional, off by default) — the
   modeled excess `E(λ) = F_thermal / F_reflected` is subtracted in
   reflectance units; requires `p_v`, `r`, `α`, and either a fixed or
   fitted beaming parameter η.

Stages that already produced their outputs are skipped on rerun (`--fresh`
overrides), so a tweak to a late stage does not re-extract everything.
Multi-block extractions also resume at block granularity: frame-number
segments whose per-frame outputs already exist are skipped.

### Procedural details per stage

**Beam-position probe (before first extraction of a new night).** Nod
positions along the slit vary night to night, and automatic aperture
finding fails on faint targets and on junk orders. The documented recipe:
extract one A−B pair of the *standard* with `load_image` + `make_profiles`,
read the ± peak positions from the order profiles, and set
`aperture_find_method="fixed"`, `aperture_find_parameter=[pos_A, pos_B]`,
`aperture_signs=["+","-"]` in the config. Signs are pinned because noisy
profiles flip the automatic sign inference on faint targets.

**Bright vs. faint targets (important robustness rule).** For targets
invisible in a single A−B pair, set `stack_object_images=True`: all
pair-subtracted frames are median-stacked into one image which is
extracted once (`reduction_mode='A'` internally) — the approach also used
by Rivkin et al. (2022) and Arredondo et al. (2024). For targets with a
visible per-pair trace, leave stacking OFF: a robust median across many
bright, seeing-variable traces clips the profile core and suppresses flux,
worst in the thermal-background orders (we measured artificial band
deepening up to a factor ~2 on bright asteroids). When in doubt, reduce
both ways and compare.

**LXD saturation knobs.** LXD sky frames commonly saturate the
longest-wavelength order; set `ignore_wavecal_saturation=True` (the
wavelength solution proceeds) and `exclude_orders="<N>"` to drop that
order from extraction. LXD wavelength calibration always needs sky frames
(the argon lamp has no lines beyond ~4 µm); by default the object frames
themselves are used (`sky_files=None`).

**File-name quirks.** `object_prefix` / `analog_prefix` override
`source_prefix` when the observer named object and standard files
differently; prefixes must include any separator (e.g. `"spc-"`,
`"objectname."`). Pre-upgrade (pre-2014) data use 4-digit frame numbers
and `instrument="spex"`.

**Residual water corrections (`spexrock/water.py`).** Beyond standard-star
division, two residual corrections are available: (1) *airmass
regression* — Beer–Lambert fit of ln(flux) vs airmass across the night's
repeated standard visits yields the night's measured optical-depth
spectrum τ(λ), applied as exp(τ·ΔAM) for the object−standard airmass
offset (model-free; preferred — in our tests it removes ~18% of
water-window residuals in a worst-case ΔAM=0.46 pair, vs ~3% for ATRAN
scaling); (2) *ATRAN power scaling* — division by a model transmission
T(λ)^x with x fitted on water-dominated windows. Fit residual corrections
on high-S/N spectra only (the standards), never on a noisy target, where
exponent fits chase noise.

**Band parameters (`spexrock/bands.py`).** Reported quantities follow the
3-µm literature: linear continuum from the clean 2.0–2.45 µm windows;
R(2.90) and R(3.05) continuum-normalized reflectances (Takir & Emery
2012 convention); band depth over 2.95–3.10 µm and at the smoothed
minimum; band center from a parabola fit around the minimum — flagged as
an *upper limit* when the minimum sits at the blue edge of the 2.86 µm
telluric cutoff (the ground-based signature of sharp-type/phyllosilicate
bands whose true center lies inside the atmospheric gap); FWHM and
integrated band area over the accessible window; bootstrap uncertainties
throughout. Ground-based caveat: the 2.50–2.86 µm region is always masked.

**Compositional fitting (`spexrock/labfit.py`).** Interpretation follows
standard curve-matching practice (Vernazza et al.; Brunetto et al.;
Reddy & Sanchez; Rivkin et al.; Fornasier et al.; Takir & Emery 2012,
Appendix C): select laboratory species plausibly present at the target's
surface temperature (phyllosilicates, carbonates, sulfates, opaques,
anhydrous silicates, carbonaceous-chondrite analogs, ammoniated species,
ices), and fit non-negative linear (areal) combinations of their
reflectance spectra to the observed spectrum, exhaustively ranking 1–3
component mixtures by reduced χ². Lab spectra are resampled and smoothed
to the data's resolution. Stated caveats on every result: linear mixing
gives areal (not mass) fractions — intimate-mixing/Hapke modeling is out
of scope; library spectra are room-temperature, while volatile band
shapes at asteroid temperatures (150–220 K) differ — cryogenic ice
spectra (SSHADE; Mastrapa et al. optical constants) are the known gap;
grain size modulates band contrast. The library-building script and index
format are documented in `data/lab_library/README.md`.

**Archival data.** The IRTF Legacy Archive serves raw FITS at
`irtfdata.ifa.hawaii.edu/idals/<sem>/data/scrs1/bigdog/<prog>/<yymmdd>/`,
searchable by object/date at `/search/`; semester schedules
(`/observing/schedule/cy<sem>_sched.txt`) map any date to program number
and PI — consult them (and the data's PI) before publishing archival
reductions.

## Install

```bash
cd SpexRock
python3.11 -m venv .venv
.venv/bin/pip install -e .
```

pyspextool downloads two large instrument calibration files from its public
S3 bucket on first use (cached in `~/Library/Caches/pyspextool`).

## Usage

```bash
spexrock-run --init my_asteroid.json   # write a template config
# edit the JSON: paths, frame numbers, analog star, mode...
spexrock-run my_asteroid.json          # reduce
spexrock-run my_asteroid.json --fresh  # full rerun, ignore intermediates
```

Or from the source tree without installing:

```bash
.venv/bin/python scripts/spexrock_run.py my_config.json
```

As a library:

```python
from spexrock import ReduceConfig
from spexrock.pipeline import run

config = ReduceConfig(
    object_name="1 Ceres", analog_name="HD 28099",
    instrument="uspex", mode="prism",
    raw_dir="raw/20260815", output_dir="spexrock_out/ceres",
    source_prefix="spc-", flat_prefix="flat-", arc_prefix="arc-",
    object_files="10-17", analog_files="20-23",
    flat_files="30-34", arc_files="35")
run(config)
```

The key config fields (see `spexrock/config.py` for all of them, with
comments):

| Field | Meaning |
| --- | --- |
| `instrument`, `mode` | `uspex`/`spex`; `prism`/`SXD`/`LXD` |
| `*_files` | pyspextool frame-number strings, e.g. `"23-24"`, `"1,3,5-9"` |
| `sky_files` | LXD only: sky frames for the long orders (default: the object frames themselves, standard IRTF practice) |
| `analog_sptype/bmag/vmag` | optional; set all three to skip the SIMBAD query |
| `normalization_wavelength_um` | default 1.20 (prism/SXD) or 3.55 (LXD) |
| `neatm_*` | thermal-excess removal; needs `p_v`, `r_au`, `alpha_deg` |

## Worked example

The repository reduces the pyspextool sample data end to end (a star pair,
so the "reflectance" is a stellar flux ratio — flat-ish, which is exactly
what makes it a good pipeline check):

```bash
git clone https://github.com/pyspextool/test_data data/sample/test_data
.venv/bin/python -m pytest tests/test_spexrock.py -m slow -v
```

![LXD sample reduction](docs/images/lxd_sample_summary.png)

## Package layout

```
spexrock/
  config.py     ReduceConfig — the one dataclass that defines a reduction
  engine.py     pyspextool call sequence (setup, cals, extract, combine)
  telluric.py   solar-analog division (reflectance mode) + airmass check
  merging.py    order merging for cross-dispersed modes
  products.py   normalization, thermal correction, FITS/ASCII/CSV writers
  thermal.py    NEATM: subsolar temperature, disk integral, excess removal
  water.py      residual telluric-water corrections (airmass regression, ATRAN)
  bands.py      formal 3-um band parameters (continuum, center, depth, width)
  labfit.py     lab-spectra mixture fitting (NNLS curve matching, RELAB library)
  qa.py         final-spectrum summary plot (pyspextool writes per-stage QA)
  pipeline.py   orchestrator with resume-from-stage
  cli/run.py    spexrock-run
scripts/        source-tree shims
tests/          test_spexrock.py (fast unit tests + slow end-to-end)
```

## Tests

```bash
.venv/bin/python -m pytest tests/test_spexrock.py -m "not slow"   # fast, no data
.venv/bin/python -m pytest tests/test_spexrock.py -m slow         # sample reduction
```

## Scope and limitations

* The analog is assumed to be a good solar match; no color correction is
  applied for non-G2V analogs (the standard caveat for this method).
* NEATM excess removal assumes the geometric albedo at the normalization
  wavelength approximates `p_v`; the default solar spectrum is a 5772 K
  Planck curve scaled to the solar constant (good to a few percent in the
  2–5 µm continuum; supply `solar_spectrum_file` for better).
* Renaming the package: the name lives only in the `spexrock/` directory
  name, `pyproject.toml`, and `PACKAGE_NAME` in `spexrock/__init__.py`.

## References

* Cushing, M. C., Vacca, W. D., & Rayner, J. T. 2004, PASP, 116, 362 — Spextool
* Vacca, W. D., Cushing, M. C., & Rayner, J. T. 2003, PASP, 115, 389 — telluric method
* Rayner, J. T., et al. 2003, PASP, 115, 362 — the SpeX instrument
* Harris, A. W. 1998, Icarus, 131, 291 — NEATM
* DeMeo, F. E., et al. 2009, Icarus, 202, 160 — 1.2 µm normalization convention
* [pyspextool](https://github.com/pyspextool/pyspextool) — the reduction engine

*Additional method references:* Vernazza, P., et al. 2015, ApJ 806, 204 & 2017, AJ 153, 72 — compositional curve matching · Reddy, V., & Sanchez, J. A., et al., Asteroids IV spectral-analysis chapter · Fornasier, S., et al. 2014, Icarus 233, 163 · Brunetto, R., et al. — laboratory analog spectra · Smette, A., et al. 2015, A&A 576, A77 (Molecfit) · Ulmer-Moll, S., et al. 2019, A&A 621, A79 (telluric-method comparison)
