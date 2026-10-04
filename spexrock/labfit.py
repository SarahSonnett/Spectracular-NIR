"""Laboratory-spectra mixture fitting for compositional interpretation.

The standard practice for interpreting asteroid 3-um spectra (Vernazza et
al.; Brunetto et al.; Reddy & Sanchez; Rivkin et al.; Fornasier et al.;
Takir & Emery 2012 Appendix C) is curve matching: select laboratory
species plausibly present at the target's temperature (volatiles and
refractories), and fit linear (areal) combinations of their reflectance
spectra to the observed spectrum.  This module implements that:

* :func:`load_library` reads the library index built from the RELAB PDS4
  bundle (see ``data/lab_library/``), resampling every lab spectrum onto
  the observed wavelength grid (convolved to the data's resolution with a
  Gaussian kernel when the lab sampling is finer);
* :func:`fit_combination` performs a non-negative least-squares (NNLS) fit
  of a given set of endmembers, with an optional multiplicative linear
  slope as a nuisance term (grain-size/space-weathering proxy);
* :func:`search_mixtures` exhaustively tries all 1-, 2-, and 3-component
  combinations from a candidate list and ranks them by reduced chi-squared
  -- the Takir & Emery Appendix-C procedure.

Interpretation caveats (state these with any result):
1. Linear mixing is AREAL; intimate (grain-scale) mixing is non-linear in
   reflectance, so NNLS weights are geographic area fractions at best, not
   mass fractions (Hapke modeling would be required for that).
2. RELAB spectra are room-temperature; band shapes of volatiles (H2O ice
   especially) shift and sharpen at asteroid temperatures (~150-220 K).
   Cryogenic ice spectra (SSHADE; Mastrapa et al. optical constants) are
   the known gap in the starter library.
3. Grain size changes band contrast; prefer matching grain-size series
   where the library has them.

References
----------
Takir, D., & Emery, J. P. 2012, Icarus 219, 641 (curve-matching procedure)
Vernazza, P., et al. 2015, ApJ 806, 204; 2017, AJ 153, 72
Reddy, V., Sanchez, J. A., et al. 2015, in Asteroids IV (spectral analysis)
Fornasier, S., et al. 2014, Icarus 233, 163
"""

from __future__ import annotations

import csv
import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import nnls

TELLURIC_GAP = (2.50, 2.86)


@dataclass
class LabSpectrum:
    species_group: str
    species: str
    sample_id: str
    spectrum_id: str
    file: Path
    wave: np.ndarray
    reflectance: np.ndarray

    def resampled(self, grid: np.ndarray, data_resolution: float = 350.0) -> np.ndarray:
        """Lab reflectance on `grid`, smoothed toward the data resolution."""
        wave, refl = self.wave, self.reflectance
        # smooth with a Gaussian whose FWHM matches the grid's R, if the lab
        # sampling is appreciably finer
        med_dl = np.median(np.diff(wave))
        target_dl = np.median(grid) / data_resolution
        if med_dl < target_dl / 3:
            sigma_pix = target_dl / 2.355 / med_dl
            n = max(int(4 * sigma_pix) | 1, 3)
            x = np.arange(n) - n // 2
            k = np.exp(-0.5 * (x / sigma_pix) ** 2)
            k /= k.sum()
            # reflect-pad so the convolution has no zero-padding edge droop
            padded = np.pad(refl, n // 2, mode="reflect")
            refl = np.convolve(padded, k, mode="valid")
        out = np.interp(grid, wave, refl, left=np.nan, right=np.nan)
        return out


def load_library(library_dir: str | Path) -> list[LabSpectrum]:
    """Load every spectrum listed in ``library_index.csv``."""
    library_dir = Path(library_dir)
    out = []
    with open(library_dir / "library_index.csv") as f:
        for row in csv.DictReader(f):
            path = library_dir / "spectra" / Path(row["file"]).name
            if not path.exists():
                continue
            table = _read_relab_csv(path)
            if table is None:
                continue
            wave, refl = table
            out.append(LabSpectrum(row.get("species_group", "?"),
                                   row.get("species", "?"),
                                   row.get("sample_id", "?"),
                                   row.get("spectrum_id", "?"),
                                   path, wave, refl))
    return out


def _read_relab_csv(path: Path):
    """RELAB PDS4 reflectance CSVs: header line(s) then wavelength,reflectance
    (wavelength may be nm or um; standard deviation column optional)."""
    rows = []
    with open(path, errors="replace") as f:
        for line in f:
            parts = line.replace(",", " ").split()
            if len(parts) >= 2:
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
    if len(rows) < 20:
        return None
    arr = np.array(rows)
    wave, refl = arr[:, 0], arr[:, 1]
    if wave.max() > 100:          # nanometers -> microns
        wave = wave / 1000.0
    order = np.argsort(wave)
    return wave[order], refl[order]


def fit_combination(grid: np.ndarray, data: np.ndarray, err: np.ndarray,
                    endmembers: list[LabSpectrum],
                    slope_nuisance: bool = True):
    """NNLS fit of data ~ sum_i w_i * lab_i(grid) [* (1 + a*(lam-lam0))].

    Returns (weights, model, reduced_chi2).  Weights are normalized areal
    fractions over the fitted endmembers; NaN regions are excluded.
    """
    columns = [em.resampled(grid) for em in endmembers]
    A = np.array(columns).T
    good = np.isfinite(data) & (err > 0) & np.all(np.isfinite(A), axis=1)
    if slope_nuisance:
        lam0 = np.nanmedian(grid[good])
        # fold a +/- linear slope into the design matrix via two
        # non-negative half-slopes applied to the mean endmember
        base = np.nanmean(A[:, :], axis=1)
        A = np.column_stack([A, base * (grid - lam0), -base * (grid - lam0)])
    Aw = A[good] / err[good, None]
    bw = data[good] / err[good]
    w, _ = nnls(Aw, bw)
    model = A @ w
    ndof = good.sum() - np.count_nonzero(w)
    chi2 = float(np.sum(((data[good] - model[good]) / err[good]) ** 2) / max(ndof, 1))
    n_em = len(endmembers)
    weights = w[:n_em]
    total = weights.sum()
    weights = weights / total if total > 0 else weights
    return weights, model, chi2


def search_mixtures(grid: np.ndarray, data: np.ndarray, err: np.ndarray,
                    library: list[LabSpectrum], max_components: int = 3,
                    top_n: int = 15, slope_nuisance: bool = True):
    """Exhaustive 1..max_components mixture search, ranked by reduced chi2.

    Returns a list of (chi2, [(species, sample_id, fraction), ...], model).
    The telluric gap is masked automatically.  Every library spectrum is
    resampled onto `grid` exactly once.
    """
    gap = (grid > TELLURIC_GAP[0]) & (grid < TELLURIC_GAP[1])
    err = np.where(gap, np.inf, err)
    columns = np.array([em.resampled(grid) for em in library])   # (nlib, n)
    good_base = np.isfinite(data) & np.isfinite(err) & (err > 0) & (err < np.inf)
    lam0 = np.nanmedian(grid[good_base])
    results = []
    for k in range(1, max_components + 1):
        for combo in itertools.combinations(range(len(library)), k):
            A = columns[list(combo)].T
            good = good_base & np.all(np.isfinite(A), axis=1)
            if good.sum() < 50:
                continue
            if slope_nuisance:
                base = np.nanmean(A, axis=1)
                A_fit = np.column_stack([A, base * (grid - lam0),
                                         -base * (grid - lam0)])
            else:
                A_fit = A
            try:
                w, _ = nnls(A_fit[good] / err[good, None], data[good] / err[good])
            except Exception:
                continue
            model = A_fit @ w
            ndof = good.sum() - np.count_nonzero(w)
            chi2 = float(np.sum(((data[good] - model[good]) / err[good]) ** 2)
                         / max(ndof, 1))
            if not np.isfinite(chi2):
                continue
            weights = w[:k]
            total = weights.sum()
            if total <= 0:
                continue
            weights = weights / total
            parts = [(library[i].species, library[i].sample_id, float(frac))
                     for i, frac in zip(combo, weights) if frac > 0.01]
            results.append((chi2, parts, model))
    results.sort(key=lambda r: r[0])
    return results[:top_n]
