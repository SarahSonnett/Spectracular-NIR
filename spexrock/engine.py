"""Thin wrappers translating a :class:`~spexrock.config.ReduceConfig` into
pyspextool calls.

Every pyspextool import in SpexRock happens here (and in the sibling
:mod:`spexrock.telluric` / :mod:`spexrock.merging` modules), so an API change
upstream is contained to these files.  The call sequence mirrors the
official uspex-prism / uspex-lxd tutorial notebooks, which are the canonical
usage of pyspextool 2.x:

``pyspextool_setup`` -> ``extract.make_flat`` -> ``extract.make_wavecal`` ->
``extract.extract`` (analog, then object) -> ``combine.combine``.

pyspextool keeps global state set by ``pyspextool_setup`` (instrument,
raw/cal/proc/qa paths, QA behavior), so :func:`setup` must run before any
other function here, and configs with different paths must not be
interleaved without re-running it.

References
----------
Cushing, M. C., Vacca, W. D., & Rayner, J. T. 2004, PASP, 116, 362
pyspextool: https://github.com/pyspextool/pyspextool
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from spexrock.config import ReduceConfig

# Headless-safe plotting unless the user already chose a backend.
os.environ.setdefault("MPLBACKEND", "Agg")

import pyspextool as ps  # noqa: E402

# pyspextool 2.0.1 inconsistency: extract() validates `aperture_signs` as a
# list but forwards it to override_aperturesigns(), which only accepts a
# comma-joined string ('+,-') — no value passes both checks.  Shim the
# override inside extract's namespace to accept either form.
import importlib  # noqa: E402

# NB: pyspextool.extract (the package) re-exports extract (the function), so
# the module itself must be fetched through importlib.
_extract_module = importlib.import_module("pyspextool.extract.extract")
from pyspextool.extract.override_aperturesigns import (  # noqa: E402
    override_aperturesigns as _override_aperturesigns)


def _override_signs_flexible(aperture_signs, verbose=None):
    if isinstance(aperture_signs, (list, tuple)):
        aperture_signs = ",".join(aperture_signs)
    return _override_aperturesigns(aperture_signs, verbose=verbose)


_extract_module.override_aperturesigns = _override_signs_flexible


def flat_name(config: ReduceConfig) -> str:
    """Name of the master flat produced by :func:`make_calibrations`."""
    return f"flat{config.flat_files}.fits"


def wavecal_name(config: ReduceConfig) -> str:
    """Name of the wavelength-calibration file produced by :func:`make_calibrations`."""
    return f"wavecal{config.arc_files}.fits"


def output_prefix(config: ReduceConfig) -> str:
    """Extraction output prefix encoding the applied corrections.

    Follows the pyspextool notebook convention so SpexRock products are
    named identically to a by-hand tutorial reduction (e.g.
    ``A-B_bc_lc_flat_wc_fxbdpx-``).
    """
    prefix = config.reduction_mode
    if config.correct_bias:
        prefix += "_bc"
    if config.correct_linearity:
        prefix += "_lc"
    if config.flat_field:
        prefix += "_flat"
    prefix += "_wc"
    if config.fix_badpixels:
        prefix += "_fxbdpx"
    return prefix + "-"


def setup(config: ReduceConfig) -> None:
    """Initialize pyspextool's global state (instrument, paths, QA policy)."""
    for directory in (config.cal_dir, config.proc_dir, config.qa_dir):
        directory.mkdir(parents=True, exist_ok=True)
    if not config.raw_dir.is_dir():
        raise FileNotFoundError(f"raw_dir does not exist: {config.raw_dir}")

    ps.pyspextool_setup(config.instrument,
                        raw_path=str(config.raw_dir) + os.sep,
                        cal_path=str(config.cal_dir) + os.sep,
                        proc_path=str(config.proc_dir) + os.sep,
                        qa_path=str(config.qa_dir) + os.sep,
                        verbose=config.verbose,
                        qa_write=config.qa_plots,
                        qa_show=False,
                        qa_showblock=False)


def make_calibrations(config: ReduceConfig) -> tuple[str, str]:
    """Build the master flat and the wavelength solution.

    For LXD the argon lamp runs out of lines beyond ~4 um, so pyspextool
    *requires* sky frames to calibrate the long orders.  The standard IRTF
    practice is to use the science frames themselves as the sky (they are
    sky-dominated at these wavelengths); that is the default here whenever
    ``config.sky_files`` is not set.
    """
    flat = flat_name(config)
    wavecal = wavecal_name(config)

    sky = None
    if config.mode == "LXD":
        sky = [config.object_file_prefix,
               config.sky_files or config.object_files]

    ps.extract.make_flat([config.flat_prefix, config.flat_files],
                         flat.removesuffix(".fits"),
                         verbose=config.verbose)

    ps.extract.make_wavecal([config.arc_prefix, config.arc_files],
                            flat,
                            wavecal.removesuffix(".fits"),
                            sky_files=sky,
                            use_stored_solution=config.use_stored_wavecal_solution,
                            ignore_saturation_error=config.ignore_wavecal_saturation,
                            verbose=config.verbose)
    return flat, wavecal


def stack_name(config: ReduceConfig) -> str:
    """File root of the combined pair-subtracted object image (faint mode)."""
    return f"{config.object_stem}_stack"


def stackspec_name(config: ReduceConfig) -> str:
    """File root of the spectra extracted from the stacked image."""
    return f"{config.object_stem}_stackspec"


def combine_object_images(config: ReduceConfig) -> str:
    """Faint-object mode step 1: robustly stack all pair-subtracted object
    frames into one image (written to the proc directory)."""
    ps.extract.combine_images([config.object_file_prefix, config.object_files],
                              stack_name(config),
                              beam_mode="A-B",
                              correct_nonlinearity=config.correct_linearity,
                              verbose=config.verbose)
    return stack_name(config)


def trace_snr_from_image(image: np.ndarray, nstrip: int = 32,
                         highpass: int = 31) -> float:
    """Strip-wise significance of narrow beam structure in a pair-subtracted
    image.

    The image is split into ``nstrip`` dispersion strips (short enough that
    order curvature is negligible); each strip is median-collapsed into a
    spatial profile and HIGH-PASSED with a running median of ``highpass``
    pixels, which removes slit-wide sky-residual order stripes while
    keeping the few-pixel-wide +/- beam peaks.  The returned value is the
    median peak significance over the stronger half of the strips.

    Calibration on real SpeX LXD nights (both detector eras): unambiguously
    bright targets score ~100-1400; faint targets and SOME bright nights
    both score ~15-50.  The metric is therefore ONE-SIDED: a high value
    proves a bright, per-pair-extractable trace, but a low value proves
    nothing (bright targets can score low through guiding smear and
    era-dependent frame structure).
    """
    sigs = []
    width = image.shape[1]
    edges = np.linspace(0, width, nstrip + 1, dtype=int)
    for a, b in zip(edges[:-1], edges[1:]):
        profile = np.nanmedian(image[:, a:b], axis=1)
        good = np.isfinite(profile)
        if good.sum() < 3 * highpass:
            continue
        profile = profile[good]
        padded = np.pad(profile, highpass // 2, mode="reflect")
        running = np.array([np.median(padded[i:i + highpass])
                            for i in range(profile.size)])
        resid = profile - running
        mad = np.median(np.abs(resid - np.median(resid))) * 1.4826
        if mad > 0:
            sigs.append(float(np.max(np.abs(resid)) / mad))
    if not sigs:
        return 0.0
    ordered = np.sort(sigs)
    return float(np.median(ordered[len(ordered) // 2:]))


def _raw_frame_path(config: ReduceConfig, prefix: str, n: int) -> Path | None:
    for pattern in (f"{prefix}{n:05d}*.fits*", f"{prefix}{n:04d}*.fits*"):
        hits = sorted(config.raw_dir.glob(pattern))
        if hits:
            return hits[0]
    return None


def object_trace_snr(config: ReduceConfig) -> float | None:
    """Trace significance of the first object A-B pair (None if unreadable)."""
    from astropy.io import fits

    numbers = expand_numbers(config.object_files)[:2]
    if len(numbers) < 2:
        return None
    paths = [_raw_frame_path(config, config.object_file_prefix, n)
             for n in numbers]
    if any(p is None for p in paths):
        return None
    frames = []
    for p in paths:
        with fits.open(p) as hdul:
            frames.append(np.asarray(hdul[0].data, dtype=float))
    if frames[0].shape != frames[1].shape:
        return None
    return trace_snr_from_image(frames[0] - frames[1])


# Above this, the per-pair trace is unambiguously detectable: every faint
# night measured scores below ~30, every score above ~100 was a bright
# target.  One-sided by design -- see trace_snr_from_image.
TRACE_SNR_BRIGHT_THRESHOLD = 60.0


def check_brightness_mode(config: ReduceConfig) -> float | None:
    """Warn when stacking is enabled on an unambiguously bright target.

    Median image stacking clips the trace cores of bright targets and
    corrupts band depths (measured: up to 2x band deepening/suppression),
    and it does so SILENTLY -- hence this gate.  The reverse mistake
    (per-pair on a faint target) fails loudly at extraction, and the metric
    cannot prove faintness anyway, so no warning is issued in that
    direction.  Advisory only: the config stays authoritative.
    """
    snr = object_trace_snr(config)
    if snr is None:
        return None
    print(f"Object trace significance (first A-B pair): {snr:.0f}")
    if snr > TRACE_SNR_BRIGHT_THRESHOLD and config.stack_object_images:
        print("WARNING: stack_object_images=True on an unambiguously bright "
              "target -- the stack's sigma-clipped combine rejects jittered "
              "trace cores and corrupts band depths; use per-pair "
              "extraction (stack_object_images=false).")
    return snr


def expand_numbers(numbers: str) -> list[int]:
    """Expand a pyspextool index string ('20-29,60,62-64') to frame numbers."""
    out: list[int] = []
    for part in numbers.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            out.extend(range(int(lo), int(hi) + 1))
        else:
            out.append(int(part))
    return out


def pending_numbers(config: ReduceConfig, numbers: str) -> str:
    """Drop comma-separated segments whose extracted spectra already exist.

    A long multi-block extraction can be interrupted; on resume, segments of
    the index string (e.g. each '20-29' block) whose per-frame output files
    are all present in the proc directory are skipped rather than redone.
    Returns the remaining index string ('' when nothing is left to do).
    """
    prefix = output_prefix(config)
    remaining = []
    for part in numbers.split(","):
        # pyspextool keeps the raw file's digit count: uspex writes 5-digit
        # frame numbers, classic spex 4-digit -- accept either
        done = all(any((config.proc_dir / f"{prefix}{n:0{d}d}.fits").exists()
                       for d in (5, 4))
                   for n in expand_numbers(part))
        if not done:
            remaining.append(part)
    return ",".join(remaining)


def extract_set(config: ReduceConfig, which: str) -> str:
    """Extract the object's or the analog's A-B pairs; returns the file numbers used.

    ``which`` is ``'object'`` or ``'analog'``.  Extraction settings
    (aperture finding, radii, background annulus, corrections) come straight
    from the config and are identical for both targets, as they must be for
    a ratio-based reflectance to cancel instrument signatures.  Frame-number
    segments already fully extracted are skipped (see :func:`pending_numbers`).
    """
    numbers = {"object": config.object_files, "analog": config.analog_files}[which]
    if not (which == "object" and config.stack_object_images):
        todo = pending_numbers(config, numbers)
        if todo != numbers:
            print(f"Resuming extraction: remaining frames {todo or '(none)'}")
        if not todo:
            return numbers
        numbers = todo

    # In 'auto' mode pyspextool's find parameter is the integer NUMBER of
    # apertures (np.empty chokes on a float there).
    find_parameter = config.aperture_find_parameter
    if config.aperture_find_method == "auto":
        find_parameter = int(find_parameter)

    stacked = which == "object" and config.stack_object_images
    if stacked:
        # The stack is a single, already pair-subtracted image: extract it in
        # 'A' mode (pyspextool's 'A-B' mode demands an even file count and
        # would try to pair-subtract again).  The +/- beam traces are handled
        # by the fixed aperture positions and pinned signs.
        reduction_mode = "A"
        files = stack_name(config) + ".fits"
        extra = {"load_directory": "proc",
                 "output_filenames": stackspec_name(config)}
    else:
        reduction_mode = config.reduction_mode
        prefix = {"object": config.object_file_prefix,
                  "analog": config.analog_file_prefix}[which]
        files = [prefix, numbers]
        extra = {}

    ps.extract.extract(reduction_mode,
                       files,
                       flat_name(config),
                       wavecal_name(config),
                       [config.aperture_find_method, find_parameter],
                       config.aperture_radius_arcsec,
                       flat_field=config.flat_field,
                       linearity_correction=config.correct_linearity,
                       output_prefix=output_prefix(config),
                       write_rectified_orders=False,
                       aperture_signs=config.aperture_signs,
                       include_orders=config.include_orders,
                       exclude_orders=config.exclude_orders,
                       bg_annulus=list(config.bg_annulus_arcsec),
                       fix_badpixels=config.fix_badpixels,
                       # object: psf_radius enables profile-weighted
                       # (optimal) extraction; analog: separately
                       # configured, default sum (bright-source profile
                       # systematics corrupt optimal weights)
                       psf_radius=(config.psf_radius_arcsec
                                   if which == "object"
                                   else config.analog_psf_radius_arcsec),
                       # bias-drift correction is a uSpeX (H2RG) detector
                       # option; pre-upgrade spex has no equivalent knob
                       detector_info=({"correct_bias": config.correct_bias}
                                      if config.instrument == "uspex" else None),
                       verbose=config.verbose,
                       **extra)
    return numbers


def combine_set(config: ReduceConfig, which: str) -> str:
    """Combine the extracted spectra for 'object' or 'analog'; returns the
    output name (without ``.fits``) written into the proc directory."""
    numbers = {"object": config.object_files, "analog": config.analog_files}[which]
    name = {"object": config.object_stem, "analog": config.analog_stem}[which]
    if which == "object" and config.stack_object_images:
        # Single stacked-extraction file: combine() averages its two (+/-)
        # apertures into the final one-aperture spectrum.
        files = stackspec_name(config) + ".fits"
    else:
        files = [output_prefix(config), numbers]
    ps.combine.combine(files, name, verbose=config.verbose)
    return name


def pyspextool_version() -> str:
    """Installed pyspextool version, for provenance headers."""
    try:
        from importlib.metadata import version
        return version("pyspextool")
    except Exception:
        return "unknown"
