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
