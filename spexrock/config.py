"""Run configuration for a SpexRock reduction.

One :class:`ReduceConfig` instance fully describes a reduction: which raw
frames make up the flats, arcs, object, and solar-analog observations, how
to extract them, how to build the reflectance spectrum, and (optionally) the
physical/observing circumstances needed for NEATM thermal-excess removal.
It replaces the "top cell" of a reduction notebook.  Configs round-trip
through JSON (:meth:`ReduceConfig.save` / :meth:`ReduceConfig.load`) so a
reduction is reproducible from a single small text file that can live next
to the data.

File-number strings use pyspextool's index syntax, e.g. ``"23-24"`` or
``"1,3,5-9"``.  Prefixes are everything before the zero-padded frame number
in the raw file names (``spc-`` for ``spc-00023.a.fits``,
``sbd.2022B046.221019.spc.`` for archive-style names).

References
----------
Cushing, M. C., Vacca, W. D., & Rayner, J. T. 2004, PASP, 116, 362
Harris, A. W. 1998, Icarus, 131, 291
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path

# Default normalization wavelength (micron) per observing mode.  1.20 um is
# the Bus-DeMeo convention for NIR asteroid spectra and sits in clean
# atmosphere for prism/SXD; LXD-only spectra normalize in the clean window
# near 3.55 um instead.
NORMALIZATION_DEFAULTS_UM = {"prism": 1.20, "SXD": 1.20, "LXD": 3.55}

CROSS_DISPERSED_MODES = ("SXD", "LXD")


@dataclass
class ReduceConfig:
    """Everything needed to reduce one asteroid (+ solar analog) data set.

    Only the names, directories, and frame-number fields are required; the
    remaining fields default to the choices used in the pyspextool tutorial
    reductions and can be tuned per data set.
    """

    # --- identities -----------------------------------------------------
    object_name: str = "object"          # asteroid, e.g. "1 Ceres"
    analog_name: str = "analog"          # solar-analog star, e.g. "HD 28099"
    instrument: str = "uspex"            # 'uspex' (post-2014) | 'spex' (pre-upgrade)
    mode: str = "prism"                  # 'prism' | 'SXD' | 'LXD'

    # --- paths ----------------------------------------------------------
    raw_dir: Path = Path("raw")          # directory holding the raw FITS frames
    output_dir: Path = Path("spexrock_out")  # cals/, proc/, qa/ created inside

    # --- raw frame selection --------------------------------------------
    source_prefix: str = "spc-"          # object AND analog frames
    # Some observers name object and standard files differently
    # (e.g. 'asteroid-' vs 'hd12345-'); these override source_prefix.
    object_prefix: str | None = None
    analog_prefix: str | None = None
    flat_prefix: str = "flat-"
    arc_prefix: str = "arc-"
    object_files: str = ""               # e.g. "23-24" (A-B pairs)
    analog_files: str = ""               # e.g. "33-34"
    flat_files: str = ""                 # e.g. "45-49"
    arc_files: str = ""                  # e.g. "43-44"
    sky_files: str | None = None         # LXD: sky frames for the long orders
    #                                      (None = use the object frames, the
    #                                      standard IRTF practice)

    # --- extraction (passed through to pyspextool.extract.extract) ------
    reduction_mode: str = "A-B"          # 'A-B' | 'A' | 'A-Sky'
    # Faint-object mode: median-stack all pair-subtracted object images into
    # one high-S/N frame and extract once, instead of extracting pair by
    # pair.  Essential when the target is invisible in a single A-B pair
    # (typical for 3-5 um asteroid work).  The analog is always extracted
    # pair by pair.
    stack_object_images: bool = False
    aperture_find_method: str = "auto"   # 'auto' | 'guess' | 'fixed'
    # For 'auto': the NUMBER of apertures to find (2 for an A-B pair, 1 for
    # staring).  For 'guess'/'fixed': the aperture position(s) in arcsec.
    aperture_find_parameter: float | list = 2
    # Fix the A/B trace signs (e.g. ["+", "-"], ordered like the aperture
    # positions) instead of letting pyspextool infer them from the profile —
    # essential for faint targets, where noisy profiles flip the inference.
    aperture_signs: list[str] | None = None
    aperture_radius_arcsec: float = 1.5
    psf_radius_arcsec: float | None = None  # None = aperture (sum) extraction
    bg_annulus_arcsec: list[float] = field(default_factory=lambda: [2.0, 2.5])
    include_orders: str | None = None    # None = all orders
    exclude_orders: str | None = None    # e.g. "4" to drop a saturated M-band order
    # LXD_long sky frames often saturate the longest-wavelength order; allow
    # the wavecal to proceed (that order should then be excluded above).
    ignore_wavecal_saturation: bool = False
    correct_bias: bool = True
    correct_linearity: bool = True
    flat_field: bool = True
    fix_badpixels: bool = True
    use_stored_wavecal_solution: bool = False  # fall-back when too few arc lines

    # --- solar analog SIMBAD bypass -------------------------------------
    # If all three are given, no network query is made; the analog is
    # described by this dictionary instead (only matters for correction
    # types other than 'reflectance', but a full bypass keeps runs offline).
    analog_sptype: str | None = None     # e.g. "G2V"
    analog_bmag: float | None = None
    analog_vmag: float | None = None

    # --- telluric / reflectance -----------------------------------------
    correction_type: str = "reflectance"  # plain object/analog ratio
    normalization_wavelength_um: float | None = None  # None = mode default
    normalization_halfwidth_um: float = 0.05
    airmass_warn_threshold: float = 0.15
    # Residual telluric-water correction applied to the merged reflectance:
    # 'none', or 'airmass' = empirical Beer-Lambert tau(lambda) regressed
    # from the night's standard visits, applied at the object-minus-standard
    # airmass offset (the method that won the comparison in the
    # water-removal investigation; it is fit on the high-S/N standards, so
    # it is safe at any target S/N).  Needs >= 3 standard visits spanning
    # >= 0.05 airmass, otherwise it is skipped with a console note.
    water_correction: str = "none"
    # Extraction weighting for the ANALOG. The analog is bright, so its
    # noise budget is dominated by profile systematics, not photons --
    # the regime where optimal extraction's assumptions fail (we measured
    # order-envelope-shaped reflectance corruption when analogs were
    # optimally extracted). Default None = sum extraction for the analog
    # regardless of psf_radius_arcsec, which then applies to the object
    # only.
    analog_psf_radius_arcsec: float | None = None

    # --- merging (cross-dispersed modes only) ---------------------------
    merge_orders: bool = True

    # --- NEATM thermal-excess removal (optional, LXD-relevant) ----------
    neatm_enabled: bool = False
    neatm_h_mag: float | None = None     # absolute magnitude H
    neatm_g_slope: float = 0.15          # phase slope G
    neatm_p_v: float | None = None       # geometric visible albedo
    neatm_eta: float = 1.0               # beaming parameter
    neatm_fit_eta: bool = False          # fit eta to the >4.5 um excess instead
    neatm_r_au: float | None = None      # heliocentric distance at epoch
    neatm_delta_au: float | None = None  # geocentric distance at epoch
    neatm_alpha_deg: float | None = None  # solar phase angle at epoch
    neatm_emissivity: float = 0.9
    solar_spectrum_file: Path | None = None  # optional; default 5772 K Planck

    # --- QA / console ---------------------------------------------------
    qa_plots: bool = True
    verbose: bool = True

    def __post_init__(self):
        self.raw_dir = Path(self.raw_dir).expanduser()
        self.output_dir = Path(self.output_dir).expanduser()
        if self.solar_spectrum_file is not None:
            self.solar_spectrum_file = Path(self.solar_spectrum_file).expanduser()
        if self.mode not in NORMALIZATION_DEFAULTS_UM:
            raise ValueError(
                f"mode must be one of {sorted(NORMALIZATION_DEFAULTS_UM)}, "
                f"got {self.mode!r}")
        if self.instrument not in ("uspex", "spex"):
            raise ValueError(f"instrument must be 'uspex' or 'spex', got {self.instrument!r}")
        if self.water_correction not in ("none", "airmass"):
            raise ValueError("water_correction must be 'none' or 'airmass', "
                             f"got {self.water_correction!r}")

    # --- derived --------------------------------------------------------
    @property
    def object_stem(self) -> str:
        """Space-free version of the object name, used in every filename
        (pyspextool strips spaces internally, so spaces break round trips)."""
        return self.object_name.replace(" ", "")

    @property
    def analog_stem(self) -> str:
        """Space-free version of the analog name, used in every filename."""
        return self.analog_name.replace(" ", "")

    @property
    def cal_dir(self) -> Path:
        return self.output_dir / "cals"

    @property
    def proc_dir(self) -> Path:
        return self.output_dir / "proc"

    @property
    def qa_dir(self) -> Path:
        return self.output_dir / "qa"

    @property
    def is_cross_dispersed(self) -> bool:
        return self.mode in CROSS_DISPERSED_MODES

    @property
    def normalization_wavelength(self) -> float:
        if self.normalization_wavelength_um is not None:
            return self.normalization_wavelength_um
        return NORMALIZATION_DEFAULTS_UM[self.mode]

    @property
    def object_file_prefix(self) -> str:
        return self.object_prefix or self.source_prefix

    @property
    def analog_file_prefix(self) -> str:
        return self.analog_prefix or self.source_prefix

    @property
    def analog_info(self) -> str | dict:
        """What to hand pyspextool as ``standard_info``: an offline dict when
        the analog's type and magnitudes were supplied, else its name (SIMBAD)."""
        if None not in (self.analog_sptype, self.analog_bmag, self.analog_vmag):
            return {"id": self.analog_name, "sptype": self.analog_sptype,
                    "bmag": self.analog_bmag, "vmag": self.analog_vmag}
        return self.analog_name

    # --- JSON round trip ------------------------------------------------
    def save(self, path: str | Path) -> Path:
        path = Path(path)
        record = asdict(self)
        for key, value in record.items():
            if isinstance(value, Path):
                record[key] = str(value)
        path.write_text(json.dumps(record, indent=2) + "\n")
        return path

    @classmethod
    def load(cls, path: str | Path) -> "ReduceConfig":
        record = json.loads(Path(path).read_text())
        return cls(**record)
