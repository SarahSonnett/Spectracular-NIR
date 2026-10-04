"""Tests for SpexRock.

Runnable either way:

    pytest tests/test_spexrock.py
    python tests/test_spexrock.py

The fast tests exercise the config, the NEATM physics, and the pipeline
plumbing without touching pyspextool or any data.  The end-to-end test
reduces the pyspextool sample data and is skipped automatically when the
cloned ``data/sample/test_data`` directory (or pyspextool itself) is
missing.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # noqa: E402

from spexrock.config import ReduceConfig  # noqa: E402
from spexrock import thermal  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
SAMPLE = REPO / "data" / "sample" / "test_data"


# ----------------------------------------------------------------------
# config
# ----------------------------------------------------------------------

def test_config_json_roundtrip(tmp_path):
    config = ReduceConfig(object_name="16 Psyche", analog_name="HD 28099",
                          mode="LXD", raw_dir=tmp_path,
                          object_files="10-13", neatm_enabled=True,
                          neatm_p_v=0.12, neatm_r_au=2.6,
                          neatm_delta_au=1.7, neatm_alpha_deg=12.0)
    path = config.save(tmp_path / "config.json")
    loaded = ReduceConfig.load(path)
    assert loaded == config
    # file is plain JSON with string paths
    record = json.loads(path.read_text())
    assert record["mode"] == "LXD"
    assert isinstance(record["raw_dir"], str)


def test_config_validation():
    with pytest.raises(ValueError):
        ReduceConfig(mode="LXD_long")   # mode must be prism|SXD|LXD
    with pytest.raises(ValueError):
        ReduceConfig(instrument="ishell")


def test_config_derived_defaults():
    assert ReduceConfig(mode="prism").normalization_wavelength == 1.20
    assert ReduceConfig(mode="LXD").normalization_wavelength == 3.55
    assert ReduceConfig(mode="LXD").is_cross_dispersed
    assert not ReduceConfig(mode="prism").is_cross_dispersed
    # analog_info: name unless the full offline dict is supplied
    assert ReduceConfig(analog_name="HD 28099").analog_info == "HD 28099"
    info = ReduceConfig(analog_name="HD 28099", analog_sptype="G2V",
                        analog_bmag=8.7, analog_vmag=8.1).analog_info
    assert info == {"id": "HD 28099", "sptype": "G2V",
                    "bmag": 8.7, "vmag": 8.1}


# ----------------------------------------------------------------------
# engine plumbing (no pyspextool calls)
# ----------------------------------------------------------------------

def test_output_prefix_matches_notebook_convention():
    from spexrock import engine
    config = ReduceConfig()
    assert engine.output_prefix(config) == "A-B_bc_lc_flat_wc_fxbdpx-"
    config = ReduceConfig(correct_bias=False, fix_badpixels=False)
    assert engine.output_prefix(config) == "A-B_lc_flat_wc-"


def test_expand_and_pending_numbers(tmp_path):
    from spexrock import engine
    assert engine.expand_numbers("20-22,60,62-63") == [20, 21, 22, 60, 62, 63]
    config = ReduceConfig(raw_dir=tmp_path, output_dir=tmp_path / "out")
    config.proc_dir.mkdir(parents=True)
    prefix = engine.output_prefix(config)
    # mark the first block fully extracted, second partially
    for n in (20, 21, 22):
        (config.proc_dir / f"{prefix}{n:05d}.fits").touch()
    (config.proc_dir / f"{prefix}00060.fits").touch()
    assert engine.pending_numbers(config, "20-22,60-61") == "60-61"
    assert engine.pending_numbers(config, "20-22") == ""


def test_calibration_names():
    from spexrock import engine
    config = ReduceConfig(flat_files="45-49", arc_files="43-44")
    assert engine.flat_name(config) == "flat45-49.fits"
    assert engine.wavecal_name(config) == "wavecal43-44.fits"


# ----------------------------------------------------------------------
# NEATM physics
# ----------------------------------------------------------------------

def test_solar_constant_recovered():
    wavelengths = np.geomspace(0.05, 1000.0, 20000)
    irradiance = thermal.solar_irradiance(wavelengths)
    total = np.trapz(irradiance, wavelengths)
    assert abs(total - thermal.SOLAR_CONSTANT) / thermal.SOLAR_CONSTANT < 0.02


def test_planck_wien_peak():
    wavelengths = np.linspace(1, 50, 20000)
    for temperature in (200.0, 300.0, 400.0):
        radiance = thermal.planck_lambda(wavelengths, temperature)
        peak = wavelengths[np.argmax(radiance)]
        assert abs(peak - 2898.0 / temperature) < 0.05 * (2898.0 / temperature)


def test_diameter_from_h():
    assert abs(thermal.diameter_km(15.0, 0.2) - 2.972) < 0.01
    # brighter (smaller H) => bigger; higher albedo => smaller
    assert thermal.diameter_km(10.0, 0.2) > thermal.diameter_km(15.0, 0.2)
    assert thermal.diameter_km(15.0, 0.4) < thermal.diameter_km(15.0, 0.2)


def test_subsolar_temperature_scalings():
    t1 = thermal.subsolar_temperature(0.1, 0.15, 1.0, 1.0)
    t2 = thermal.subsolar_temperature(0.1, 0.15, 1.0, 4.0)
    assert abs(t1 / t2 - 2.0) < 1e-6          # T ~ r^-1/2
    hot = thermal.subsolar_temperature(0.1, 0.15, 0.8, 1.0)
    assert hot > t1                            # lower beaming => hotter
    assert 380 < t1 < 420                      # ~394 K for these parameters


def test_phase_function_limits():
    assert thermal.hg_phase_function(0.0, 0.15) == pytest.approx(1.0)
    assert 0.0 < thermal.hg_phase_function(30.0, 0.15) < 1.0


def test_disk_integral_blackbody_limit():
    # With mu^(1/4) replaced by 1 the integral at alpha=0 must be pi*B.
    # Emulate by evaluating at a wavelength deep in the Rayleigh-Jeans tail
    # where B varies slowly with T, so T ~ T_ss over most of the disk.
    wavelength = np.array([2000.0])            # um, extreme RJ regime
    t_ss = 400.0
    integral = thermal.neatm_disk_integral(wavelength, t_ss, 0.0, n_grid=256)
    limit = np.pi * thermal.planck_lambda(wavelength, t_ss)
    ratio = float(integral[0] / limit[0])
    # RJ: B ~ T, disk-average of mu^(1/4) with projection weight is 8/9
    assert abs(ratio - 8.0 / 9.0) < 0.02


def test_thermal_excess_behavior():
    wavelengths = np.linspace(1.5, 5.0, 100)
    excess = thermal.thermal_excess(wavelengths, p_v=0.15, g_slope=0.15,
                                    eta=1.0, r_au=1.5, alpha_deg=20.0)
    assert np.all(np.isfinite(excess))
    assert excess[0] < 1e-3                    # negligible at 1.5 um
    assert excess[-1] > excess[0]              # grows to the red
    assert np.all(np.diff(excess) >= -1e-12)   # monotonic
    cooler = thermal.thermal_excess(wavelengths, p_v=0.15, g_slope=0.15,
                                    eta=1.5, r_au=1.5, alpha_deg=20.0)
    assert cooler[-1] < excess[-1]             # higher eta => less excess


def test_remove_excess_roundtrip():
    wavelengths = np.linspace(0.8, 5.0, 400)
    truth = np.ones_like(wavelengths)
    excess = thermal.thermal_excess(wavelengths, p_v=0.2, g_slope=0.15,
                                    eta=1.0, r_au=1.1, alpha_deg=30.0)
    observed = truth + excess
    corrected, model = thermal.remove_excess(
        wavelengths, observed, None, p_v=0.2, g_slope=0.15, eta=1.0,
        r_au=1.1, alpha_deg=30.0)
    assert np.allclose(corrected, truth, atol=1e-10)
    assert np.allclose(model, excess)


def test_remove_excess_eta_fit_recovers():
    wavelengths = np.linspace(0.8, 5.0, 500)
    truth = np.ones_like(wavelengths)
    excess = thermal.thermal_excess(wavelengths, p_v=0.2, g_slope=0.15,
                                    eta=1.3, r_au=1.1, alpha_deg=30.0)
    observed = truth + excess
    corrected, _ = thermal.remove_excess(
        wavelengths, observed, None, p_v=0.2, g_slope=0.15,
        eta=1.0, fit_eta=True, r_au=1.1, alpha_deg=30.0)
    residual = np.abs(corrected[wavelengths > 4.0] - 1.0)
    assert np.max(residual) < 0.02


# ----------------------------------------------------------------------
# water residual corrections
# ----------------------------------------------------------------------

def test_airmass_regression_recovers_tau():
    from spexrock import water
    rng = np.random.default_rng(5)
    grid = np.linspace(2.0, 3.4, 400)
    tau_true = 0.3 * np.exp(-0.5 * ((grid - 3.0) / 0.05) ** 2)  # fake line
    visits = []
    for am in (1.05, 1.2, 1.4, 1.6):
        flux = np.exp(-tau_true * am) * (1 + rng.normal(0, 0.003, grid.size))
        visits.append((grid, flux, am))
    tau = water.airmass_regression(visits, grid)
    good = np.isfinite(tau)
    assert np.nanmax(np.abs(tau[good] - tau_true[good])) < 0.03
    # applying the correction undoes a known mismatch
    refl = np.exp(-tau_true * 0.3)          # object 0.3 airmass deeper
    fixed = water.apply_airmass_correction(grid, refl, grid, tau, 0.3)
    assert np.nanmax(np.abs(fixed[good] - 1.0)) < 0.05


def test_atran_power_scaling_recovers_exponent():
    from spexrock import water
    rng = np.random.default_rng(6)
    aw, at = water.load_atran(2000)
    # resample to a realistic R~2500-like pixel grid (the shipped ATRAN
    # files are heavily oversampled)
    wave = np.arange(1.8, 3.6, 5e-4)
    trans = np.clip(np.interp(wave, aw, at), 1e-3, 1)
    x_true = 0.25
    refl = trans ** x_true * (1 + rng.normal(0, 0.01, wave.size))
    err = np.full(wave.size, 0.01)
    corrected, x_fit = water.atran_power_scaling(wave, refl, err)
    assert abs(x_fit - x_true) < 0.08
    band = (wave > 2.86) & (wave < 3.25)
    assert np.nanstd(corrected[band]) < np.nanstd(refl[band])


# ----------------------------------------------------------------------
# products
# ----------------------------------------------------------------------

def _fake_spectra(norders=2, npix=50):
    spectra = np.full((norders, 4, npix), np.nan)
    for i in range(norders):
        wave = np.linspace(1.0 + i, 2.0 + i, npix)
        spectra[i, 0] = wave
        spectra[i, 1] = 2.0 * (1 + 0.1 * (wave - 1.5))   # sloped, level 2
        spectra[i, 2] = 0.02
        spectra[i, 3] = 0
    return spectra


def test_normalize_at_wavelength():
    from spexrock import products
    spectra, factor = products.normalize(_fake_spectra(), 1.5, 0.05)
    window = np.abs(spectra[0, 0] - 1.5) <= 0.05
    assert np.nanmedian(spectra[0, 1][window]) == pytest.approx(1.0, abs=1e-6)
    assert factor == pytest.approx(2.0, rel=0.02)
    # uncertainties scale with the flux
    assert np.nanmedian(spectra[0, 2]) == pytest.approx(0.02 / factor, rel=1e-6)


def test_normalize_empty_window_falls_back():
    from spexrock import products
    spectra, factor = products.normalize(_fake_spectra(norders=1), 5.0, 0.01)
    assert np.isfinite(factor) and factor > 0


# ----------------------------------------------------------------------
# band parameters + lab fitting
# ----------------------------------------------------------------------

def _synthetic_band(depth=0.25, center=3.00, width=0.18, n=6000, noise=0.01,
                    seed=4):
    rng = np.random.default_rng(seed)
    wave = np.linspace(1.9, 4.1, n)
    refl = 1.0 - depth * np.exp(-0.5 * ((wave - center) / width) ** 2)
    refl *= 1.0 + 0.05 * (wave - 2.2)          # mild red slope
    err = np.full(n, noise)
    return wave, refl + rng.normal(0, noise, n), err


def test_band_measure_recovers_synthetic():
    from spexrock import bands
    wave, refl, err = _synthetic_band()
    p = bands.measure(wave, refl, err, n_boot=120)
    assert abs(p.center_um - 3.00) < 0.03
    assert not p.center_is_limit
    assert abs(p.depth_at_minimum - 0.25) < 0.03
    assert abs(p.r290 - (1 - 0.25 * np.exp(-0.5 * ((2.90 - 3.0) / 0.18) ** 2))) < 0.03
    assert p.depth_err < 0.02
    assert np.isfinite(p.area_um) and p.area_um > 0


def test_band_center_limit_flag():
    from spexrock import bands
    # band centered inside the telluric gap -> accessible minimum at the
    # blue edge -> center flagged as a limit (sharp-type convention)
    wave, refl, err = _synthetic_band(center=2.75, width=0.25)
    p = bands.measure(wave, refl, err, n_boot=60)
    assert p.center_is_limit


def test_labfit_recovers_mixture(tmp_path):
    from spexrock import labfit
    grid = np.linspace(2.0, 4.0, 900)
    em_a = 1.0 - 0.5 * np.exp(-0.5 * ((grid - 2.95) / 0.12) ** 2)   # "serpentine"
    em_b = np.full_like(grid, 0.6)                                   # "opaque"
    specs = []
    for name, refl in (("serp", em_a), ("opaque", em_b)):
        specs.append(labfit.LabSpectrum("grp", name, name, name,
                                        tmp_path / f"{name}.csv",
                                        grid.copy(), refl.copy()))
    rng = np.random.default_rng(8)
    truth = 0.7 * em_a + 0.3 * em_b
    data = truth + rng.normal(0, 0.005, grid.size)
    err = np.full(grid.size, 0.005)
    weights, model, chi2 = labfit.fit_combination(grid, data, err, specs,
                                                  slope_nuisance=False)
    assert abs(weights[0] - 0.7) < 0.05
    assert abs(weights[1] - 0.3) < 0.05
    assert chi2 < 2.0
    ranked = labfit.search_mixtures(grid, data, err, specs, max_components=2,
                                    top_n=3)
    assert len(ranked[0][1]) == 2      # best model uses both endmembers


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def test_cli_template_init(tmp_path):
    from spexrock.cli.run import main
    config_path = tmp_path / "template.json"
    assert main([str(config_path), "--init"]) == 0
    loaded = ReduceConfig.load(config_path)
    assert loaded == ReduceConfig()
    # refuses to clobber
    assert main([str(config_path), "--init"]) == 1


# ----------------------------------------------------------------------
# end-to-end on the pyspextool sample data (slow; needs data + network-free)
# ----------------------------------------------------------------------

LXD_RAW = SAMPLE / "raw" / "uspex-LXD" / "data"
PRISM_RAW = SAMPLE / "raw" / "uspex-prism" / "data"


@pytest.mark.slow
@pytest.mark.skipif(not LXD_RAW.is_dir(),
                    reason="sample data not present (data/sample/test_data)")
def test_end_to_end_uspex_lxd(tmp_path):
    pytest.importorskip("pyspextool")
    from spexrock import pipeline

    config = ReduceConfig(
        object_name="CD-2317771A", analog_name="HD223352",
        instrument="uspex", mode="LXD",
        raw_dir=LXD_RAW, output_dir=tmp_path / "out",
        source_prefix="spc-", flat_prefix="flat-", arc_prefix="arc-",
        object_files="23-24", analog_files="33-34",
        flat_files="45-49", arc_files="43-44",
        analog_sptype="B8V", analog_bmag=4.53, analog_vmag=4.57,
        normalization_wavelength_um=2.2)
    written = pipeline.run(config)
    finals = [p for p in written if p.suffix == ".dat"]
    assert finals, "no final .dat product written"
    table = np.loadtxt(finals[0])
    wave, refl = table[:, 0], table[:, 1]
    assert wave.min() < 2.0 and wave.max() > 4.0
    near = np.abs(wave - 2.2) < 0.05
    assert np.nanmedian(refl[near]) == pytest.approx(1.0, abs=0.05)


@pytest.mark.slow
@pytest.mark.skipif(not PRISM_RAW.is_dir(),
                    reason="sample data not present (data/sample/test_data)")
def test_end_to_end_uspex_prism(tmp_path):
    pytest.importorskip("pyspextool")
    from spexrock import pipeline

    config = ReduceConfig(
        object_name="J2010-1707", analog_name="HD193689",
        instrument="uspex", mode="prism",
        raw_dir=PRISM_RAW, output_dir=tmp_path / "out",
        source_prefix="sbd.2022B046.221019.spc.",
        flat_prefix="sbd.2022B046.221019.flat.",
        arc_prefix="sbd.2022B046.221019.arc.",
        object_files="1-2", analog_files="7-8",
        flat_files="15-19", arc_files="20",
        analog_sptype="A0V", analog_bmag=7.1, analog_vmag=7.1)
    written = pipeline.run(config)
    finals = [p for p in written if p.suffix == ".dat"]
    assert finals, "no final .dat product written"
    table = np.loadtxt(finals[0])
    wave, refl = table[:, 0], table[:, 1]
    assert wave.min() < 0.9 and wave.max() > 2.4
    near = np.abs(wave - 1.20) < 0.05
    assert np.nanmedian(refl[near]) == pytest.approx(1.0, abs=0.05)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "-m", "not slow"]))
