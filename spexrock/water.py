"""Residual telluric-water corrections beyond standard-star division.

Dividing an asteroid spectrum by a solar-analog observed at a slightly
different airmass (and time) leaves residual water absorption or emission
structure, worst in the 2.86-3.05 um wing where the 3-um band science
lives.  This module implements two residual corrections:

1. :func:`airmass_regression` -- a model-free, night-specific method.  The
   night's repeated standard-star visits span a range of airmasses; for
   each wavelength, Beer-Lambert absorption makes ln(flux) linear in
   airmass, so regressing the (continuum-normalized) visits yields the
   night's actual optical-depth spectrum tau(lambda).  The reflectance is
   then corrected by exp(tau * dAM) for the object-minus-standard airmass
   difference.  No atmospheric model or line list is involved; the
   atmosphere is measured from the data themselves.

2. :func:`atran_power_scaling` -- the classical model route (the
   Bus/Volquardsen approach in spirit): divide by an ATRAN transmission
   model raised to a fitted power x, where x absorbs the precipitable-water
   and airmass mismatch in the weak-line limit (T^x scaling).  x is fitted
   by minimizing telluric-correlated structure in the corrected spectrum
   within water-dominated windows.

Both operate on a 1D (wavelength, reflectance, error) spectrum after the
pipeline's telluric division, so they compose with any reduction route.

References
----------
Lord, S. D. 1992, NASA Tech. Mem. 103957 (ATRAN)
Smette, A., et al. 2015, A&A 576, A77 (Molecfit; the full-fitting approach)
Ulmer-Moll, S., et al. 2019, A&A 621, A79 (method comparison)
Rivkin, A. S., et al. 2022, PSJ 3, 153 (ATRAN residual step in LMNOP)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# Windows dominated by telluric H2O (avoiding the opaque 2.5-2.86 core and
# the asteroid band minimum region as much as practical).
WATER_WINDOWS = [(2.86, 3.02), (3.10, 3.25)]
CLEAN_WINDOW = (2.10, 2.40)


def _continuum_normalize(wave: np.ndarray, flux: np.ndarray) -> np.ndarray:
    """Divide out a low-order continuum fit over the clean window region."""
    sel = (wave >= CLEAN_WINDOW[0]) & (wave <= CLEAN_WINDOW[1]) & (flux > 0)
    coef = np.polyfit(wave[sel], np.log(flux[sel]), 1)
    return flux / np.exp(np.polyval(coef, wave))


def airmass_regression(visits: list[tuple[np.ndarray, np.ndarray, float]],
                       grid: np.ndarray) -> np.ndarray:
    """Fit tau(lambda) from repeated standard visits.

    Parameters
    ----------
    visits : list of (wavelength_um, flux, airmass) per standard visit
        (same star throughout the night; 3+ visits spanning airmass).
    grid : wavelength grid (um) on which to return tau.

    Returns
    -------
    tau : ndarray
        Optical depth relative to the clean-window continuum, per `grid`
        wavelength; NaN where fewer than 3 visits contribute.
    """
    if len(visits) < 3:
        raise ValueError("airmass regression needs >= 3 standard visits")
    logf = np.full((len(visits), grid.size), np.nan)
    ams = np.array([am for _, _, am in visits], dtype=float)
    for i, (w, f, _) in enumerate(visits):
        good = np.isfinite(w) & np.isfinite(f) & (f > 0)
        fn = _continuum_normalize(w[good], f[good])
        inside = (grid >= w[good].min()) & (grid <= w[good].max())
        logf[i, inside] = np.log(np.interp(grid[inside], w[good], fn))
    tau = np.full(grid.size, np.nan)
    for j in range(grid.size):
        ok = np.isfinite(logf[:, j])
        if ok.sum() >= 3 and np.ptp(ams[ok]) > 0.05:
            slope = np.polyfit(ams[ok], logf[ok, j], 1)[0]
            tau[j] = -slope
    return tau


def apply_airmass_correction(wave: np.ndarray, reflectance: np.ndarray,
                             tau_grid: np.ndarray, tau: np.ndarray,
                             dam: float) -> np.ndarray:
    """Correct a reflectance for an object-minus-standard airmass offset.

    reflectance_true = reflectance_obs * exp(tau * dam), with
    dam = AM_object - AM_standard (positive when the object was observed
    through more atmosphere than the standard).
    """
    tau_i = np.interp(wave, tau_grid, tau, left=np.nan, right=np.nan)
    factor = np.exp(np.nan_to_num(tau_i, nan=0.0) * dam)
    return reflectance * factor


def night_airmass_correction(config, spectra: np.ndarray) -> np.ndarray:
    """Apply the empirical airmass-regression water correction to a merged
    reflectance ``spectra`` array (pyspextool layout: order x plane x pixel).

    tau(lambda) is regressed from the night's individually extracted
    standard-star frames (read back from the proc directory with their
    header airmasses) and applied at dam = <AM_object> - <AM_analog>.
    Returns the corrected array; when the night cannot support the
    regression (< 3 analog visits, or airmass span < 0.05) the input is
    returned unchanged with a console note.
    """
    from astropy.io import fits

    from spexrock import engine

    def per_frame(numbers: str):
        prefix = engine.output_prefix(config)
        out = []
        for n in engine.expand_numbers(numbers):
            for d in (5, 4):
                path = config.proc_dir / f"{prefix}{n:0{d}d}.fits"
                if path.exists():
                    with fits.open(path) as hdul:
                        data = np.asarray(hdul[0].data, dtype=float)
                        am = float(hdul[0].header["AM"])
                    wave = data[..., 0, :].ravel()
                    flux = data[..., 1, :].ravel()
                    good = np.isfinite(wave) & np.isfinite(flux)
                    order = np.argsort(wave[good])
                    out.append((wave[good][order], flux[good][order], am))
                    break
        return out

    visits = per_frame(config.analog_files)
    ams = np.array([am for _, _, am in visits])
    if len(visits) < 3 or (len(visits) and np.ptp(ams) < 0.05):
        print("water_correction=airmass: night cannot support the "
              f"regression ({len(visits)} analog frames, airmass span "
              f"{np.ptp(ams) if len(visits) else 0:.3f}); skipping.")
        return spectra

    object_ams = [am for _, _, am in per_frame(config.object_files)]
    if not object_ams:
        print("water_correction=airmass: no per-frame object spectra found; "
              "skipping.")
        return spectra
    dam = float(np.mean(object_ams) - np.mean(ams))

    grid = np.geomspace(min(v[0].min() for v in visits),
                        max(v[0].max() for v in visits), 4000)
    tau = airmass_regression(visits, grid)
    print(f"water_correction=airmass: tau from {len(visits)} analog frames "
          f"(AM {ams.min():.2f}-{ams.max():.2f}), dam = {dam:+.3f}")

    corrected = np.array(spectra, dtype=float, copy=True)
    for i in range(corrected.shape[0]):
        wave, refl = corrected[i, 0, :], corrected[i, 1, :]
        good = np.isfinite(wave) & np.isfinite(refl)
        if good.any():
            corrected[i, 1, good] = apply_airmass_correction(
                wave[good], refl[good], grid, tau, dam)
    return corrected


def load_atran(resolution: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    """Load a shipped pyspextool ATRAN transmission model (wave_um, T)."""
    from astropy.io import fits
    import pyspextool
    path = (Path(pyspextool.__file__).parent / "data"
            / f"atran{resolution}.fits")
    with fits.open(path) as hdul:
        data = hdul[0].data
    return data[0], data[1]


def atran_power_scaling(wave: np.ndarray, reflectance: np.ndarray,
                        error: np.ndarray,
                        resolution: int = 2000,
                        x_bounds: tuple[float, float] = (-0.7, 0.7),
                        ) -> tuple[np.ndarray, float]:
    """Fit and remove residual water via T_atran^x division.

    Returns (corrected reflectance, fitted exponent x).  x ~ 0 means the
    standard-star division already cancelled the water; positive x means
    residual absorption remained (object seen through more water).
    The fit minimizes the local scatter of the corrected spectrum around a
    running median inside the water windows, so a smooth real band shape is
    preserved while line-correlated structure is penalized.
    """
    from scipy.optimize import minimize_scalar

    aw, at = load_atran(resolution)
    trans = np.interp(wave, aw, at)
    trans = np.clip(trans, 1e-3, 1.0)

    sel = np.zeros(wave.size, dtype=bool)
    for lo, hi in WATER_WINDOWS:
        sel |= (wave >= lo) & (wave <= hi)
    sel &= np.isfinite(reflectance) & (error > 0) & (trans > 0.05)
    if sel.sum() < 50:
        return reflectance, 0.0

    order = np.argsort(wave[sel])
    wsel = wave[sel][order]
    # Smoothing window in WAVELENGTH (broader than telluric line blends at
    # R~2000, narrower than real band structure), robust to any sampling.
    half = 0.015  # um
    lo_idx = np.searchsorted(wsel, wsel - half)
    hi_idx = np.searchsorted(wsel, wsel + half)

    def cost(x: float) -> float:
        corr = (reflectance[sel] / trans[sel] ** x)[order]
        run = np.array([np.median(corr[lo_idx[i]:hi_idx[i]])
                        for i in range(corr.size)])
        # relative deviations: immune to the overall rescaling that T^x
        # division applies in low-transmission regions
        safe = np.abs(run) > 1e-6
        return float(np.nanmedian(np.abs(corr[safe] / run[safe] - 1.0)))

    result = minimize_scalar(cost, bounds=x_bounds, method="bounded")
    x = float(result.x)
    return reflectance / trans ** x, x
