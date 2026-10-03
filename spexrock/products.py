"""Final reflectance products: normalization, thermal correction, and writers.

The pipeline's upstream stages leave a telluric-corrected (and, for
cross-dispersed modes, merged) spectrum in pyspextool's FITS layout --
an array of shape (norders x napertures, 4, npix) holding wavelength [um],
flux, uncertainty, and bitmask per order.  This module turns that into the
science product an asteroid spectroscopist actually uses: a single
reflectance spectrum normalized to unity at the configured wavelength,
optionally NEATM-corrected, written as

* ``<object>_final.fits``  -- normalized spectrum in the same pyspextool
  layout (readable by pyspextool tools), provenance in the header;
* ``<object>_final.dat``   -- four-column ASCII (wavelength, reflectance,
  uncertainty, flag) with a commented provenance header;
* ``<object>_final.csv``   -- the same table for spreadsheet users.

Normalization divides flux and uncertainty by the median flux inside
+/- ``normalization_halfwidth_um`` of the normalization wavelength; if that
window is empty (e.g. LXD-only data normalized at a prism wavelength) the
median of the cleanest quartile of the spectrum is used instead, with a
warning.

References
----------
Cushing, M. C., Vacca, W. D., & Rayner, J. T. 2004, PASP, 116, 362
DeMeo, F. E., et al. 2009, Icarus, 202, 160 (1.2 um normalization convention)
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from astropy.io import fits

from spexrock.config import ReduceConfig
from spexrock import thermal as thermal_mod
from spexrock import engine


def load_spectrum(path: Path) -> tuple[np.ndarray, dict]:
    """Read a pyspextool spectra FITS; returns (spectra, info)."""
    from pyspextool.io.read_spectra_fits import read_spectra_fits
    return read_spectra_fits(str(path))


def normalize(spectra: np.ndarray, wavelength_um: float,
              halfwidth_um: float) -> tuple[np.ndarray, float]:
    """Normalize all orders by the median flux near ``wavelength_um``.

    Returns (normalized spectra, normalization factor).  Operates across
    orders so inter-order scaling from the merge stage is preserved.
    """
    spectra = np.array(spectra, dtype=float, copy=True)
    wave = spectra[:, 0, :]
    flux = spectra[:, 1, :]
    window = (np.abs(wave - wavelength_um) <= halfwidth_um) & np.isfinite(flux)

    if np.any(window):
        factor = float(np.nanmedian(flux[window]))
    else:
        finite = flux[np.isfinite(flux)]
        factor = float(np.nanmedian(
            finite[finite >= np.nanpercentile(finite, 75)]))
        print(f"WARNING: no data within {halfwidth_um} um of "
              f"{wavelength_um} um; normalized to the brightest-quartile "
              "median instead.")
    if factor == 0 or not np.isfinite(factor):
        raise ValueError("Normalization factor is zero or non-finite; "
                         "check the reflectance spectrum.")

    spectra[:, 1, :] /= factor
    spectra[:, 2, :] /= factor
    return spectra, factor


def apply_thermal_correction(spectra: np.ndarray,
                             config: ReduceConfig) -> tuple[np.ndarray, np.ndarray]:
    """NEATM thermal-excess removal on normalized spectra (per order).

    Returns ``(corrected_spectra, excess_spectra)`` where the second array
    mirrors the first but carries the modeled excess in its flux plane (for
    QA plotting)."""
    required = {"neatm_p_v": config.neatm_p_v, "neatm_r_au": config.neatm_r_au,
                "neatm_alpha_deg": config.neatm_alpha_deg}
    missing = [k for k, v in required.items() if v is None]
    if missing:
        raise ValueError("NEATM enabled but missing config fields: "
                         + ", ".join(missing))

    spectra = np.array(spectra, dtype=float, copy=True)
    excess_spectra = np.array(spectra, dtype=float, copy=True)
    excess_spectra[:, 1, :] = np.nan
    for i in range(spectra.shape[0]):
        wave, flux, err = spectra[i, 0, :], spectra[i, 1, :], spectra[i, 2, :]
        good = np.isfinite(wave) & np.isfinite(flux)
        if not np.any(good):
            continue
        corrected, excess = thermal_mod.remove_excess(
            wave[good], flux[good], err[good],
            p_v=config.neatm_p_v, g_slope=config.neatm_g_slope,
            eta=config.neatm_eta, r_au=config.neatm_r_au,
            alpha_deg=config.neatm_alpha_deg,
            emissivity=config.neatm_emissivity,
            fit_eta=config.neatm_fit_eta,
            spectrum_file=config.solar_spectrum_file)
        spectra[i, 1, good] = corrected
        excess_spectra[i, 1, good] = excess
        if np.nanmax(excess) > 0.01:
            print(f"NEATM: peak thermal excess {np.nanmax(excess):.3f} "
                  f"(reflectance units) at {wave[good][np.nanargmax(excess)]:.2f} um")
    return spectra, excess_spectra


def final_name(config: ReduceConfig) -> str:
    return f"{config.object_stem}_final"


def write_products(spectra: np.ndarray, info: dict, config: ReduceConfig,
                   source_path: Path, normalization_factor: float) -> list[Path]:
    """Write the final FITS/ASCII/CSV products; returns the paths written."""
    root = config.proc_dir / final_name(config)
    written: list[Path] = []

    # FITS: reuse the source file so pyspextool's header/layout survive.
    with fits.open(source_path) as hdul:
        hdul[0].data = spectra.astype(np.float32)
        header = hdul[0].header
        header["NORMWAVE"] = (config.normalization_wavelength,
                              "Normalization wavelength (um)")
        header["NORMFAC"] = (normalization_factor,
                             "Flux divided by this factor")
        header["THERMCOR"] = (bool(config.neatm_enabled),
                              "NEATM thermal excess removed")
        from spexrock import __version__
        header["SPEXROCK"] = (__version__, "SpexRock version")
        header["PYSPEXT"] = (engine.pyspextool_version(), "pyspextool version")
        fits_path = root.with_suffix(".fits")
        hdul.writeto(fits_path, overwrite=True)
    written.append(fits_path)

    # Flatten orders into one wavelength-sorted table.
    wave = spectra[:, 0, :].ravel()
    flux = spectra[:, 1, :].ravel()
    err = spectra[:, 2, :].ravel()
    flag = spectra[:, 3, :].ravel()
    good = np.isfinite(wave)
    order = np.argsort(wave[good])
    table = np.column_stack([wave[good][order], flux[good][order],
                             err[good][order], flag[good][order]])

    provenance = {
        "object": config.object_name,
        "analog": config.analog_name,
        "instrument": config.instrument,
        "mode": config.mode,
        "object_frames": config.object_files,
        "analog_frames": config.analog_files,
        "normalization_wavelength_um": config.normalization_wavelength,
        "neatm_thermal_correction": config.neatm_enabled,
        "spexrock_version": None,
        "pyspextool_version": engine.pyspextool_version(),
    }
    from spexrock import __version__
    provenance["spexrock_version"] = __version__

    header_lines = ["SpexRock reflectance spectrum",
                    *(f"{k}: {v}" for k, v in provenance.items()),
                    "columns: wavelength_um reflectance uncertainty flag"]
    dat_path = root.with_suffix(".dat")
    np.savetxt(dat_path, table, fmt="%12.6f %12.6f %12.6f %6d",
               header="\n".join(header_lines))
    written.append(dat_path)

    csv_path = root.with_suffix(".csv")
    np.savetxt(csv_path, table, fmt="%.6f", delimiter=",",
               header="wavelength_um,reflectance,uncertainty,flag",
               comments="")
    written.append(csv_path)

    # Config snapshot next to the products for exact reproducibility.
    config_path = config.output_dir / "config_used.json"
    record = {k: (str(v) if isinstance(v, Path) else v)
              for k, v in asdict(config).items()}
    config_path.write_text(json.dumps(record, indent=2) + "\n")
    written.append(config_path)

    return written
