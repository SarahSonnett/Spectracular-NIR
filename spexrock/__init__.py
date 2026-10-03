"""SpexRock: asteroid reflectance spectra from NASA IRTF SpeX.

SpexRock is a thin, batch-friendly pipeline layer over :mod:`pyspextool`,
the official Python implementation of the Spextool reduction protocols
(Cushing, Vacca & Rayner 2004).  pyspextool does the heavy lifting --
nonlinearity and bias corrections, flat-field construction with order-edge
tracing, wavelength calibration from argon arcs (plus sky emission lines in
LXD mode), pair subtraction, spatial-profile-weighted optimal extraction,
spectrum combination, telluric correction, and order merging.  SpexRock
drives that machinery from a single JSON-serializable configuration
(:class:`~spexrock.config.ReduceConfig`, which replaces the "notebook top
cell") and adds the asteroid-specific steps that Spextool leaves to the
user: division by a solar-analog star (pyspextool's ``'reflectance'``
telluric mode) with normalization to unity at a chosen wavelength, tidy
reflectance products with full provenance, and optional NEATM thermal-excess
removal for wavelengths beyond ~2.5 microns (:mod:`spexrock.thermal`).

Supported observing modes are SpeX prism and cross-dispersed SXD/LXD for
both the pre-upgrade (Aladdin) and upgraded (H2RG) instrument eras --
whatever the installed pyspextool supports.

The pipeline stages, in order (:mod:`spexrock.pipeline`):

1. calibrations -- master flat and wavelength solution
   (:func:`spexrock.engine.make_calibrations`);
2. extraction -- object and solar analog, A-B pair mode
   (:func:`spexrock.engine.extract_set`);
3. combination of the extracted spectra (:func:`spexrock.engine.combine_set`);
4. telluric correction by solar-analog division
   (:func:`spexrock.telluric.correct`);
5. order merging for cross-dispersed modes (:func:`spexrock.merging.merge`);
6. optional NEATM thermal-excess removal (:func:`spexrock.thermal.remove_excess`);
7. final normalized-reflectance products (:mod:`spexrock.products`).

References
----------
Cushing, M. C., Vacca, W. D., & Rayner, J. T. 2004, PASP, 116, 362
Vacca, W. D., Cushing, M. C., & Rayner, J. T. 2003, PASP, 115, 389
Rayner, J. T., et al. 2003, PASP, 115, 362 (the SpeX instrument)
Harris, A. W. 1998, Icarus, 131, 291 (NEATM)
pyspextool: https://github.com/pyspextool/pyspextool
"""

from __future__ import annotations

PACKAGE_NAME = "spexrock"  # single place the name lives besides pyproject.toml
__version__ = "0.1.0"

from spexrock.config import ReduceConfig  # noqa: E402

__all__ = ["ReduceConfig", "PACKAGE_NAME", "__version__"]
