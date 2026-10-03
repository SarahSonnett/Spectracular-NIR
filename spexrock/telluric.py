"""Telluric correction by solar-analog division.

For asteroid work the telluric standard is a solar analog, and the telluric
correction and the conversion to reflectance are the same operation: the
ratio asteroid/analog cancels both the atmospheric transmission and the
solar spectrum, leaving (up to normalization, and up to the analog's
deviation from the Sun's colors) the relative reflectance.  This is
pyspextool's ``correction_type='reflectance'`` — the Python descendant of
IDL Spextool's ``xtellcor_basic`` division, including pyspextool's
wavelength-shift optimization between object and standard that minimizes
telluric residuals at the band edges.

The airmass difference between asteroid and analog is checked here and
loudly warned about: at |dAM| beyond ~0.1-0.15 the residual telluric slope
starts to compete with real spectral slope differences between taxonomic
classes.

References
----------
Vacca, W. D., Cushing, M. C., & Rayner, J. T. 2003, PASP, 115, 389
Cushing, M. C., Vacca, W. D., & Rayner, J. T. 2004, PASP, 116, 362
"""

from __future__ import annotations

from pathlib import Path

from astropy.io import fits

from spexrock.config import ReduceConfig
from spexrock import engine  # ensures MPLBACKEND is set before pyspextool loads

import pyspextool as ps  # noqa: E402


def reflectance_name(config: ReduceConfig) -> str:
    """Output root of the telluric-corrected (reflectance) spectrum."""
    return f"{config.object_stem}_reflectance"


def _airmass(path: Path) -> float | None:
    """Average airmass recorded in a combined-spectrum FITS header, if any."""
    try:
        header = fits.getheader(path)
    except OSError:
        return None
    for key in ("AVE_AM", "AM", "TCS_AM", "AIRMASS"):
        if key in header:
            try:
                return float(header[key])
            except (TypeError, ValueError):
                continue
    return None


def correct(config: ReduceConfig) -> str:
    """Divide the combined object spectrum by the combined analog spectrum.

    Returns the reflectance file root written into the proc directory.
    The result is an unnormalized relative reflectance; normalization to
    unity at the config wavelength happens in :mod:`spexrock.products`.
    """
    object_file = config.object_stem + ".fits"
    analog_file = config.analog_stem + ".fits"

    am_object = _airmass(config.proc_dir / object_file)
    am_analog = _airmass(config.proc_dir / analog_file)
    if am_object is not None and am_analog is not None:
        delta = abs(am_object - am_analog)
        print(f"Airmass: object {am_object:.3f}, analog {am_analog:.3f} "
              f"(|dAM| = {delta:.3f})")
        if delta > config.airmass_warn_threshold:
            print(f"WARNING: |dAM| = {delta:.3f} exceeds "
                  f"{config.airmass_warn_threshold:.2f}; residual telluric "
                  "slopes may masquerade as spectral slope.")

    ps.telluric.telluric(object_file,
                         analog_file,
                         config.analog_info,
                         f"telluric_{config.object_files}",
                         reflectance_name(config),
                         correction_type=config.correction_type,
                         verbose=config.verbose)
    return reflectance_name(config)
