"""Order merging for cross-dispersed (SXD/LXD) reductions.

Cross-dispersed SpeX spectra come out of extraction as separate orders with
small wavelength overlaps.  pyspextool's merge module scales adjacent
orders in their overlap regions and combines them onto a single wavelength
grid; prism spectra are single-order and skip this stage.  Merging runs on
the telluric-corrected spectrum, so the order-scale factors are computed on
reflectance rather than on raw counts (safer where the overlaps are
telluric-contaminated, which for LXD is most of them).

References
----------
Cushing, M. C., Vacca, W. D., & Rayner, J. T. 2004, PASP, 116, 362
"""

from __future__ import annotations

from spexrock.config import ReduceConfig
from spexrock import engine  # noqa: F401  (sets MPLBACKEND before pyspextool)

import pyspextool as ps  # noqa: E402


def merged_name(config: ReduceConfig) -> str:
    """Output root of the merged spectrum."""
    return f"{config.object_stem}_reflectance_merged"


def merge(config: ReduceConfig, reflectance_root: str) -> str:
    """Merge the orders of the telluric-corrected spectrum.

    Returns the file root of the spectrum the rest of the pipeline should
    consume: the merged file for cross-dispersed modes, or the input
    unchanged for prism (nothing to merge).
    """
    if not (config.is_cross_dispersed and config.merge_orders):
        return reflectance_root

    ps.merge.merge(reflectance_root + ".fits",
                   merged_name(config),
                   verbose=config.verbose)
    return merged_name(config)
