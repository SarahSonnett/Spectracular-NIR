"""The SpexRock pipeline orchestrator.

:func:`run` executes the full raw-to-reflectance sequence described in the
package docstring, printing a stage banner as it goes.  With
``resume=True`` (the default) stages whose output files already exist in
the run's output directory are skipped, so an interrupted or tweaked
reduction picks up where it left off; pass ``resume=False`` (or delete the
output directory) for a clean rerun.

References
----------
Cushing, M. C., Vacca, W. D., & Rayner, J. T. 2004, PASP, 116, 362
Vacca, W. D., Cushing, M. C., & Rayner, J. T. 2003, PASP, 115, 389
"""

from __future__ import annotations

from pathlib import Path

from spexrock.config import ReduceConfig
from spexrock import engine, merging, products, qa, telluric


def _banner(text: str) -> None:
    print(f"\n=== SpexRock: {text} ===")


def run(config: ReduceConfig, resume: bool = True) -> list[Path]:
    """Run the pipeline; returns the list of final product paths."""
    _banner(f"{config.object_name} [{config.instrument} {config.mode}] "
            f"analog {config.analog_name}")
    engine.setup(config)

    # --- calibrations ----------------------------------------------------
    flat = config.cal_dir / engine.flat_name(config)
    wavecal = config.cal_dir / engine.wavecal_name(config)
    if resume and flat.exists() and wavecal.exists():
        _banner("calibrations exist; skipping (resume)")
    else:
        _banner("calibrations (flat + wavecal)")
        engine.make_calibrations(config)

    # --- extraction + combination ---------------------------------------
    for which in ("analog", "object"):
        name = {"object": config.object_stem,
                "analog": config.analog_stem}[which]
        combined = config.proc_dir / f"{name}.fits"
        if resume and combined.exists():
            _banner(f"{which} already combined; skipping (resume)")
            continue
        if which == "object" and config.stack_object_images:
            stack_file = config.proc_dir / (engine.stack_name(config) + ".fits")
            if resume and stack_file.exists():
                _banner("object image stack exists; skipping (resume)")
            else:
                _banner("stacking object images (faint mode)")
                engine.combine_object_images(config)
        _banner(f"extracting {which}")
        engine.extract_set(config, which)
        _banner(f"combining {which}")
        engine.combine_set(config, which)

    # --- telluric correction / reflectance -------------------------------
    reflectance_root = telluric.reflectance_name(config)
    reflectance_file = config.proc_dir / (reflectance_root + ".fits")
    if resume and reflectance_file.exists():
        _banner("reflectance exists; skipping telluric (resume)")
    else:
        _banner("telluric correction (solar-analog division)")
        telluric.correct(config)

    # --- order merging ----------------------------------------------------
    if config.is_cross_dispersed and config.merge_orders:
        merged_file = config.proc_dir / (merging.merged_name(config) + ".fits")
        if resume and merged_file.exists():
            _banner("merged spectrum exists; skipping merge (resume)")
            final_root = merging.merged_name(config)
        else:
            _banner("merging orders")
            final_root = merging.merge(config, reflectance_root)
    else:
        final_root = reflectance_root
    spectrum_path = config.proc_dir / (final_root + ".fits")

    # --- normalization, thermal correction, products ---------------------
    _banner("building final products")
    spectra, info = products.load_spectrum(spectrum_path)
    spectra, factor = products.normalize(spectra,
                                         config.normalization_wavelength,
                                         config.normalization_halfwidth_um)
    print(f"Normalized at {config.normalization_wavelength:.3f} um "
          f"(factor {factor:.5g})")

    excess_spectra = None
    if config.neatm_enabled:
        _banner("NEATM thermal-excess removal")
        spectra, excess_spectra = products.apply_thermal_correction(spectra, config)

    written = products.write_products(spectra, info, config,
                                      spectrum_path, factor)
    if config.qa_plots:
        written.append(qa.summary_plot(spectra, config, excess_spectra))

    _banner("done")
    for path in written:
        print(f"  wrote {path}")
    return written
