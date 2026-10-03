"""Quality-assurance summary plotting.

pyspextool already writes detailed per-stage QA files (flat edges, arc line
identifications, spatial profiles, traces, telluric shifts) into the run's
``qa/`` directory; this module adds the one figure those don't provide --
the final normalized reflectance spectrum with its uncertainty envelope and,
when NEATM correction ran, the modeled thermal excess -- so a reduction can
be judged at a glance from a single image.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from spexrock.config import ReduceConfig


def summary_plot(spectra: np.ndarray, config: ReduceConfig,
                 excess_spectra: np.ndarray | None = None) -> Path:
    """Plot the final spectrum; returns the path of the saved PNG."""
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 5))
    for i in range(spectra.shape[0]):
        wave, flux, err = spectra[i, 0, :], spectra[i, 1, :], spectra[i, 2, :]
        good = np.isfinite(wave) & np.isfinite(flux)
        ax.plot(wave[good], flux[good], lw=0.8,
                color="C0" if spectra.shape[0] == 1 else None)
        ax.fill_between(wave[good], (flux - err)[good], (flux + err)[good],
                        alpha=0.25, lw=0, color="C0")
    if excess_spectra is not None:
        for i in range(excess_spectra.shape[0]):
            wave, excess = excess_spectra[i, 0, :], excess_spectra[i, 1, :]
            good = np.isfinite(wave) & np.isfinite(excess) & (excess > 0)
            if np.any(good):
                ax.plot(wave[good], excess[good], ls="--", color="C3",
                        lw=0.9, label="NEATM thermal excess" if i == 0 else None)

    ax.axhline(1.0, color="0.7", lw=0.5, zorder=0)
    ax.axvline(config.normalization_wavelength, color="0.7", lw=0.5, zorder=0)
    ax.set_xlabel("Wavelength ($\\mu$m)")
    ax.set_ylabel(f"Reflectance (= 1 at {config.normalization_wavelength:.2f} $\\mu$m)")
    ax.set_title(f"{config.object_name}  /  {config.analog_name}   "
                 f"[{config.instrument} {config.mode}]")
    finite = spectra[:, 1, :][np.isfinite(spectra[:, 1, :])]
    if finite.size:
        lo, hi = np.nanpercentile(finite, [1, 99])
        ax.set_ylim(max(lo - 0.2, -0.5), hi + 0.2)
    if excess_spectra is not None:
        ax.legend(frameon=False)
    fig.tight_layout()

    path = config.qa_dir / f"{config.object_stem}_final_summary.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
