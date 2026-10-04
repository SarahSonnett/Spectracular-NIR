"""Formal 3-um band-parameter measurement: continuum, center, depth, width.

This replaces the ad-hoc window metrics used during development with a
single documented routine, returning the quantities the 3-um literature
actually tabulates (Takir & Emery 2012; Rivkin et al. 2015, 2022):

* a linear continuum fitted to clean short-wavelength windows;
* R(2.90) and R(3.05): continuum-normalized reflectance at the standard
  comparison wavelengths;
* band depth at 2.95-3.10 um (survey convention) and at the band minimum;
* band center from a parabola fit around the smoothed minimum -- reported
  as an UPPER LIMIT (`center_is_limit`) when the minimum lies at the blue
  edge of the accessible window, the ground-based signature of sharp-type
  (phyllosilicate) bands whose true center sits inside the 2.5-2.85 um
  telluric gap;
* band width (FWHM) and integrated band area over the accessible window;
* bootstrap uncertainties for everything.

The input spectrum should already be telluric-corrected, thermal-corrected
(see :mod:`spexrock.thermal`), and masked of the 2.5-2.86 um telluric gap;
:func:`measure` applies the gap mask itself for safety.

References
----------
Takir, D., & Emery, J. P. 2012, Icarus 219, 641 (band groups; R2.90 metric)
Rivkin, A. S., et al. 2015, AJ 150, 198; 2022, PSJ 3, 153
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

TELLURIC_GAP = (2.50, 2.86)
CONTINUUM_WINDOWS = [(2.00, 2.45)]
BAND_SEARCH = (2.86, 3.35)
DEPTH_WINDOW = (2.95, 3.10)
SHARP_EDGE = 2.885  # minimum this close to the gap edge => center is a limit


@dataclass
class BandParameters:
    r290: float
    r290_err: float
    r305: float
    r305_err: float
    depth: float              # 1 - <R> over DEPTH_WINDOW
    depth_err: float
    depth_at_minimum: float
    center_um: float
    center_err: float
    center_is_limit: bool     # True => quote as "< center_um" (sharp-type)
    fwhm_um: float            # NaN when the band never recovers to half depth
    area_um: float            # integral of (1 - R) over the accessible band
    continuum_coeffs: tuple

    def as_dict(self) -> dict:
        return asdict(self)


def _bin(wave, flux, err, R=350.0, lo=1.9, hi=4.1):
    edges = np.geomspace(lo, hi, int(np.log(hi / lo) * R))
    c, m, s = [], [], []
    for i in range(len(edges) - 1):
        sel = (wave >= edges[i]) & (wave < edges[i + 1])
        if sel.sum() > 3:
            w = 1.0 / err[sel] ** 2
            c.append(np.mean(wave[sel]))
            m.append(np.sum(w * flux[sel]) / np.sum(w))
            s.append(np.sqrt(1.0 / np.sum(w)))
    return np.array(c), np.array(m), np.array(s)


def _single_measure(c, m, s):
    cont_sel = np.zeros(c.size, dtype=bool)
    for lo, hi in CONTINUUM_WINDOWS:
        cont_sel |= (c >= lo) & (c <= hi)
    coeffs = np.polyfit(c[cont_sel], m[cont_sel], 1, w=1.0 / s[cont_sel])
    ratio = m / np.polyval(coeffs, c)

    def r_at(lam, hw=0.04):
        sel = np.abs(c - lam) < hw
        return np.average(ratio[sel], weights=1.0 / s[sel] ** 2) if sel.any() else np.nan

    r290, r305 = r_at(2.90), r_at(3.05)
    ev = (c >= DEPTH_WINDOW[0]) & (c <= DEPTH_WINDOW[1])
    depth = 1.0 - np.average(ratio[ev], weights=1.0 / s[ev] ** 2)

    band = (c >= BAND_SEARCH[0]) & (c <= BAND_SEARCH[1])
    cb, rb = c[band], ratio[band]
    kernel = np.ones(7) / 7.0
    rs = np.convolve(rb, kernel, mode="same")
    i_min = int(np.argmin(rs[3:-3]) + 3)
    # parabola through +-3 points around the smoothed minimum
    sl = slice(max(i_min - 3, 0), min(i_min + 4, cb.size))
    p = np.polyfit(cb[sl], rs[sl], 2)
    center = -p[1] / (2 * p[0]) if p[0] > 0 else cb[i_min]
    center = float(np.clip(center, cb.min(), cb.max()))
    depth_min = 1.0 - float(rs[i_min])
    is_limit = center <= SHARP_EDGE

    # FWHM / area over the accessible band (blue side may be truncated)
    half = 1.0 - depth_min / 2.0
    below = rs <= half
    fwhm = np.nan
    if below.any():
        idx = np.where(below)[0]
        fwhm = float(cb[idx[-1]] - cb[idx[0]]) if idx.size > 1 else np.nan
    area = float(np.trapz(np.clip(1.0 - rb, 0, None), cb))
    return (r290, r305, depth, depth_min, center, is_limit, fwhm, area, coeffs)


def measure(wave: np.ndarray, reflectance: np.ndarray, error: np.ndarray,
            n_boot: int = 300, seed: int = 0) -> BandParameters:
    """Measure 3-um band parameters with bootstrap uncertainties."""
    rng = np.random.default_rng(seed)
    good = (np.isfinite(wave) & np.isfinite(reflectance) & (error > 0)
            & ~((wave > TELLURIC_GAP[0]) & (wave < TELLURIC_GAP[1]))
            & (np.abs(reflectance) < 20))
    c, m, s = _bin(wave[good], reflectance[good], error[good])
    base = _single_measure(c, m, s)
    boots = []
    for _ in range(n_boot):
        boots.append(_single_measure(c, m + rng.normal(0, s), s)[:5])
    boots = np.array([b[:5] for b in boots], dtype=float)
    return BandParameters(
        r290=base[0], r290_err=float(np.nanstd(boots[:, 0])),
        r305=base[1], r305_err=float(np.nanstd(boots[:, 1])),
        depth=base[2], depth_err=float(np.nanstd(boots[:, 2])),
        depth_at_minimum=base[3],
        center_um=base[4], center_err=float(np.nanstd(boots[:, 4])),
        center_is_limit=bool(base[5]),
        fwhm_um=base[6], area_um=base[7],
        continuum_coeffs=tuple(np.asarray(base[8]).tolist()))
