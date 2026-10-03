"""NEATM thermal-excess estimation and removal.

Beyond ~2.5 um an asteroid's observed flux is reflected sunlight plus its
own thermal emission; a reflectance spectrum built as object/solar-analog
therefore turns up at the long-wavelength end unless the thermal
contribution is removed.  This module implements the Near-Earth Asteroid
Thermal Model (Harris 1998) and expresses its emission as an additive
"thermal excess" in normalized-reflectance units,

    R_obs(lambda) = R_true(lambda) + E(lambda),

where E is the ratio of the NEATM thermal flux to the reflected-light flux
that a unit-reflectance surface would produce.  The reflected continuum is
modeled with the observing geometry (r, Delta, alpha), the IAU H-G phase
function, and a solar spectrum (a 5772 K Planck curve scaled to the solar
constant by default -- good to a few percent in the 2-5 um continuum -- or
a user-supplied two-column file, wavelength [um] vs irradiance at 1 AU
[W m-2 um-1]).

In NEATM the dayside temperature is T(mu) = T_ss * mu^(1/4) with mu the
cosine of the angular distance from the subsolar point, zero flux from the
nightside, and the subsolar temperature

    T_ss = [ (1 - A) * S / (eta * eps * sigma * r^2) ]^(1/4),

with A = q(G) * p_v the Bond albedo (q = 0.290 + 0.684 G), eta the beaming
parameter, eps the emissivity.  Fluxes come from numerically integrating
the Planck radiance over the illuminated-and-visible portion of the sphere
seen at phase angle alpha.

The key assumption for spectra normalized in the near-IR: the geometric
albedo at the normalization wavelength approximates p_v.  For red or very
blue objects this misestimates the reflected continuum by the corresponding
color factor; the beaming parameter (or `fit_eta`) absorbs much of the
error in practice.

References
----------
Harris, A. W. 1998, Icarus, 131, 291 (NEATM)
Bowell, E., et al. 1989, in Asteroids II, 524 (H-G system, q(G))
Rivkin, A. S., et al. 2005, Icarus, 175, 175 (thermal excess in reflectance)
Delbo, M., & Harris, A. W. 2002, M&PS, 37, 1929 (NEATM practice)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# Physical constants (SI)
H_PLANCK = 6.62607015e-34   # J s
C_LIGHT = 2.99792458e8      # m s-1
K_BOLTZ = 1.380649e-23      # J K-1
SIGMA_SB = 5.670374419e-8   # W m-2 K-4
SOLAR_CONSTANT = 1361.0     # W m-2 at 1 AU
T_SUN = 5772.0              # K, IAU nominal
R_SUN_M = 6.957e8           # m
AU_M = 1.495978707e11       # m


def planck_lambda(wavelength_um: np.ndarray | float, temperature: float) -> np.ndarray:
    """Planck spectral radiance B_lambda in W m-2 um-1 sr-1."""
    wavelength_m = np.asarray(wavelength_um, dtype=float) * 1e-6
    with np.errstate(over="ignore", divide="ignore"):
        exponent = H_PLANCK * C_LIGHT / (wavelength_m * K_BOLTZ * max(temperature, 1e-6))
        radiance_m = (2.0 * H_PLANCK * C_LIGHT**2 / wavelength_m**5
                      / np.expm1(np.clip(exponent, None, 700.0)))
    return radiance_m * 1e-6  # per meter -> per micron


def solar_irradiance(wavelength_um: np.ndarray,
                     spectrum_file: Path | None = None) -> np.ndarray:
    """Solar spectral irradiance at 1 AU, W m-2 um-1.

    Default: pi * B_lambda(5772 K) * (R_sun/AU)^2, whose integral is the
    solar constant by construction of the IAU nominal values.
    """
    wavelength_um = np.asarray(wavelength_um, dtype=float)
    if spectrum_file is not None:
        table = np.loadtxt(spectrum_file)
        return np.interp(wavelength_um, table[:, 0], table[:, 1])
    return np.pi * planck_lambda(wavelength_um, T_SUN) * (R_SUN_M / AU_M) ** 2


def diameter_km(h_mag: float, p_v: float) -> float:
    """Effective diameter from absolute magnitude and geometric albedo,
    D = 1329 km * 10^(-H/5) / sqrt(p_v) (Fowler & Chillemi 1992 convention)."""
    return 1329.0 * 10.0 ** (-h_mag / 5.0) / np.sqrt(p_v)


def bond_albedo(p_v: float, g_slope: float) -> float:
    """Bond albedo A = q * p_v with the H-G phase integral q = 0.290 + 0.684 G."""
    return (0.290 + 0.684 * g_slope) * p_v


def subsolar_temperature(p_v: float, g_slope: float, eta: float,
                         r_au: float, emissivity: float = 0.9) -> float:
    """NEATM subsolar temperature in K."""
    albedo = bond_albedo(p_v, g_slope)
    return ((1.0 - albedo) * SOLAR_CONSTANT
            / (eta * emissivity * SIGMA_SB * r_au**2)) ** 0.25


def hg_phase_function(alpha_deg: float, g_slope: float) -> float:
    """IAU H-G phase function Phi(alpha), unity at zero phase (Bowell 1989)."""
    alpha = np.radians(abs(alpha_deg))
    tan_half = np.tan(alpha / 2.0)
    phi1 = np.exp(-3.33 * tan_half**0.63)
    phi2 = np.exp(-1.87 * tan_half**1.22)
    return float((1.0 - g_slope) * phi1 + g_slope * phi2)


def neatm_disk_integral(wavelength_um: np.ndarray, t_ss: float,
                        alpha_deg: float, n_grid: int = 128) -> np.ndarray:
    """The NEATM angular integral I(lambda), in W m-2 um-1 sr-1.

    I = integral over the illuminated-and-visible surface of
    B_lambda(T_ss * mu^(1/4)) * cos^2(phi) * cos(theta - alpha) dtheta dphi,
    with theta longitude from the subsolar meridian, phi latitude, and the
    observer in the equatorial plane at longitude alpha.  The thermal flux
    density at Earth is then  F = emissivity * (D/2)^2 / Delta^2 * I.
    A blackbody disk sanity check: at alpha = 0 and uniform temperature the
    integral reduces to pi * B_lambda(T).
    """
    wavelength_um = np.asarray(wavelength_um, dtype=float)
    alpha = np.radians(alpha_deg)
    # Illuminated: |theta| < pi/2. Visible: |theta - alpha| < pi/2.
    theta_lo = alpha - np.pi / 2.0 if alpha >= 0 else -np.pi / 2.0
    theta_hi = np.pi / 2.0 if alpha >= 0 else alpha + np.pi / 2.0
    theta = np.linspace(theta_lo, theta_hi, n_grid)
    phi = np.linspace(-np.pi / 2.0, np.pi / 2.0, n_grid)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")

    mu_sun = np.clip(np.cos(phi_grid) * np.cos(theta_grid), 0.0, None)
    mu_obs = np.clip(np.cos(phi_grid) * np.cos(theta_grid - alpha), 0.0, None)
    temperature = t_ss * mu_sun**0.25

    integral = np.empty_like(wavelength_um)
    weight = np.cos(phi_grid) * mu_obs  # cos(phi) Jacobian * observer projection
    dtheta = theta[1] - theta[0]
    dphi = phi[1] - phi[0]
    for i, wl in enumerate(wavelength_um):
        radiance = _planck_grid(wl, temperature)
        integral[i] = np.sum(radiance * weight) * dtheta * dphi
    return integral


def _planck_grid(wavelength_um: float, temperature: np.ndarray) -> np.ndarray:
    """Planck radiance (W m-2 um-1 sr-1) for one wavelength over a temperature grid."""
    wavelength_m = wavelength_um * 1e-6
    out = np.zeros_like(temperature)
    warm = temperature > 1.0
    with np.errstate(over="ignore"):
        exponent = (H_PLANCK * C_LIGHT
                    / (wavelength_m * K_BOLTZ * temperature[warm]))
        out[warm] = (2.0 * H_PLANCK * C_LIGHT**2 / wavelength_m**5
                     / np.expm1(np.clip(exponent, None, 700.0))) * 1e-6
    return out


def thermal_excess(wavelength_um: np.ndarray,
                   p_v: float,
                   g_slope: float,
                   eta: float,
                   r_au: float,
                   alpha_deg: float,
                   emissivity: float = 0.9,
                   spectrum_file: Path | None = None) -> np.ndarray:
    """Thermal excess E(lambda) in normalized-reflectance units.

    E = F_thermal / F_reflected(R=1), where the unit-reflectance reflected
    flux is F_ref = p_v * Phi(alpha) * S_sun(lambda)/r^2 * (D/2)^2/Delta^2.
    The diameter and Delta cancel in the ratio, so E needs no H or Delta:

        E(lambda) = emissivity * I(lambda) * r^2 / (p_v * Phi(alpha) * S_sun(lambda))
    """
    t_ss = subsolar_temperature(p_v, g_slope, eta, r_au, emissivity)
    integral = neatm_disk_integral(wavelength_um, t_ss, alpha_deg)
    reflected = (p_v * hg_phase_function(alpha_deg, g_slope)
                 * solar_irradiance(wavelength_um, spectrum_file) / r_au**2)
    return emissivity * integral / reflected


def remove_excess(wavelength_um: np.ndarray,
                  reflectance: np.ndarray,
                  uncertainty: np.ndarray | None,
                  p_v: float,
                  g_slope: float,
                  eta: float,
                  r_au: float,
                  alpha_deg: float,
                  emissivity: float = 0.9,
                  fit_eta: bool = False,
                  spectrum_file: Path | None = None,
                  ) -> tuple[np.ndarray, np.ndarray]:
    """Subtract the NEATM thermal excess from a normalized reflectance spectrum.

    Returns ``(corrected_reflectance, excess)``.  With ``fit_eta`` the
    beaming parameter is adjusted (bounded 0.6-3.0) so that the corrected
    reflectance beyond 3.5 um is consistent with a linear continuum
    extrapolated from the 2.0-3.5 um region -- the usual practice when eta
    is not known from radiometry.  The uncertainty array is unchanged by the
    subtraction (the model is treated as exact) but is used to weight the
    eta fit.
    """
    wavelength_um = np.asarray(wavelength_um, dtype=float)
    reflectance = np.asarray(reflectance, dtype=float)

    def excess_for(eta_value: float) -> np.ndarray:
        return thermal_excess(wavelength_um, p_v, g_slope, eta_value, r_au,
                              alpha_deg, emissivity, spectrum_file)

    if fit_eta:
        from scipy.optimize import minimize_scalar

        short = (wavelength_um > 2.0) & (wavelength_um < 3.5)
        long_ = wavelength_um >= 3.5
        good = np.isfinite(reflectance)
        if uncertainty is not None:
            weight = np.where(np.asarray(uncertainty) > 0,
                              1.0 / np.square(uncertainty), 0.0)
        else:
            weight = np.ones_like(reflectance)

        if np.sum(short & good) > 10 and np.sum(long_ & good) > 10:

            def cost(eta_value: float) -> float:
                # Refit the continuum on the *corrected* spectrum each trial,
                # otherwise the thermal contamination of the 2-3.5 um region
                # biases the continuum high and the fit under-subtracts.
                corrected = reflectance - excess_for(eta_value)
                coeffs = np.polyfit(wavelength_um[short & good],
                                    corrected[short & good], 1)
                continuum = np.polyval(coeffs, wavelength_um)
                resid = corrected - continuum
                mask = long_ & good
                return float(np.sum(weight[mask] * resid[mask] ** 2))

            result = minimize_scalar(cost, bounds=(0.6, 3.0), method="bounded")
            eta = float(result.x)
            print(f"NEATM: fitted beaming parameter eta = {eta:.2f}")
        else:
            print("NEATM: too few points beyond 3.5 um to fit eta; "
                  f"using eta = {eta:.2f}")

    excess = excess_for(eta)
    return reflectance - excess, excess
