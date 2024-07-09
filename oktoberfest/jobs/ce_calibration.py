import logging
from pathlib import Path
from typing import List, Optional, Union

import numpy as np
from sklearn.linear_model import LinearRegression, RANSACRegressor

from oktoberfest import __copyright__, __version__
from oktoberfest import plotting as pl
from oktoberfest import predict as pr
from oktoberfest import preprocessing as pp
from oktoberfest import rescore as re

from ..data.spectra import Spectra
from ..utils import Config, JobPool, ProcessStep, group_iterator

logger = logging.getLogger(__name__)


def preprocess(spectra_files: List[Path], config: Config) -> List[Path]:
    preprocess_search_step = ProcessStep(config.output, "preprocessing_search")
    if not preprocess_search_step.is_done():
        # load search results
        if not config.search_results_type == "internal":
            logger.info(f"Converting search results from {config.search_results} to internal search result.")

            msms_output = config.output / "msms"
            msms_output.mkdir(exist_ok=True)
            internal_search_file = msms_output / "msms.prosit"
            tmt_label = config.tag

            search_results = pp.convert_search(
                input_path=config.search_results,
                search_engine=config.search_results_type,
                tmt_label=tmt_label,
                output_file=internal_search_file,
            )
            if config.spectra_type.lower() in ["d", "hdf"]:
                timstof_metadata = pp.convert_timstof_metadata(
                    input_path=config.search_results,
                    search_engine=config.search_results_type,
                    output_file=msms_output / "tims_meta.csv",
                )
        else:
            internal_search_file = config.search_results
            search_results = pp.load_search(internal_search_file)
            # TODO add support for internal timstof metadata
        logger.info(f"Read {len(search_results)} PSMs from {internal_search_file}")

        # filter search results
        search_results = pp.filter_peptides_for_model(peptides=search_results, model=config.models["intensity"])

        # split search results
        searchfiles_found = pp.split_search(
            search_results=search_results,
            output_dir=config.output / "msms",
            filenames=[spectra_file.stem for spectra_file in spectra_files],
        )
        # split timstof metadata
        if config.spectra_type.lower() in ["d", "hdf"]:
            _ = pp.split_timstof_metadata(
                timstof_metadata=timstof_metadata,
                output_dir=config.output / "msms",
                filenames=searchfiles_found,
            )
        preprocess_search_step.mark_done()
    else:
        searchfiles_found = [msms_file.stem for msms_file in (config.output / "msms").glob("*rescore")]
    spectra_files_to_return = []
    for spectra_file in spectra_files:
        if spectra_file.stem in searchfiles_found:
            spectra_files_to_return.append(spectra_file)

    return spectra_files_to_return


def _annotate_and_get_library(spectra_file: Path, config: Config, tims_meta_file: Optional[Path] = None) -> Spectra:
    data_dir = config.output / "data"
    data_dir.mkdir(exist_ok=True)
    hdf5_path = data_dir / spectra_file.with_suffix(".mzml.hdf5").name
    if hdf5_path.is_file():
        aspec = Spectra.from_hdf5(hdf5_path)
        instrument_type = config.instrument_type
        if instrument_type is not None and aspec.obs["INSTRUMENT_TYPES"].values[0] != instrument_type:
            aspec.obs["INSTRUMENT_TYPES"] = instrument_type
            aspec.write_as_hdf5(hdf5_path)
    else:
        spectra_dir = config.output / "spectra"
        spectra_dir.mkdir(exist_ok=True)
        format_ = spectra_file.suffix.lower()
        if format_ == ".raw":
            file_to_load = spectra_dir / spectra_file.with_suffix(".mzML").name
            pp.convert_raw_to_mzml(spectra_file, file_to_load, thermo_exe=config.thermo_exe)
        elif format_ in [".mzml", ".hdf"]:
            file_to_load = spectra_file
        elif format_ == ".d":
            file_to_load = spectra_dir / spectra_file.with_suffix(".hdf").name
            pp.convert_d_to_hdf(spectra_file, file_to_load)
        spectra = pp.load_spectra(file_to_load, tims_meta_file=tims_meta_file)
        config_instrument_type = config.instrument_type
        if config_instrument_type is not None:
            spectra["INSTRUMENT_TYPES"] = config_instrument_type
        search = pp.load_search(config.output / "msms" / spectra_file.with_suffix(".rescore").name)
        library = pp.merge_spectra_and_peptides(spectra, search)
        aspec = pp.annotate_spectral_library(
            library, mass_tol=config.mass_tolerance, unit_mass_tol=config.unit_mass_tolerance
        )
        aspec.write_as_hdf5(hdf5_path)  # write_metadata_annotation

    return aspec


def _get_best_ce(library: Spectra, spectra_file: Path, config: Config):
    results_dir = config.output / "results"
    results_dir.mkdir(exist_ok=True)
    if (library.obs["FRAGMENTATION"] == "HCD").any():
        server_kwargs = {
            "server_url": config.prediction_server,
            "ssl": config.ssl,
            "model_name": config.models["intensity"],
        }
        use_ransac_model = config.use_ransac_model
        alignment_library = pr.ce_calibration(library, config.ce_range, use_ransac_model, **server_kwargs)

        if use_ransac_model:
            logger.info("Performing RANSAC regression")
            calib_group = (
                alignment_library.obs.groupby(
                    by=["PRECURSOR_CHARGE", "ORIG_COLLISION_ENERGY", "COLLISION_ENERGY", "MASS"], as_index=False
                )["SPECTRAL_ANGLE"]
                .mean()
                .groupby(["PRECURSOR_CHARGE", "ORIG_COLLISION_ENERGY", "MASS"], as_index=False)
                .apply(lambda x: x.loc[x["SPECTRAL_ANGLE"].idxmax()])
            )
            calib_group["delta_collision_energy"] = (
                calib_group["COLLISION_ENERGY"] - calib_group["ORIG_COLLISION_ENERGY"]
            )
            x = calib_group[["MASS", "PRECURSOR_CHARGE"]]  # input feature
            y = calib_group["delta_collision_energy"]  # target variable
            ransac = RANSACRegressor(LinearRegression(), residual_threshold=1.5, random_state=42)
            ransac.fit(x, y)

            for charge, df in calib_group.groupby("PRECURSOR_CHARGE"):
                r2_score = ransac.score(df[["MASS", "PRECURSOR_CHARGE"]], df["COLLISION_ENERGY"])
                title = f"Scatter Plot with RANSAC Model \nSlope: {ransac.estimator_.coef_[0]:.2f}, "
                title += f"Intercept: {ransac.estimator_.intercept_:.2f}, R2: {r2_score:.2f}"
                pl.plot_ce_ransac_model(
                    sa_ce_df=df,
                    filename=results_dir / f"{spectra_file.stem}_ce_ransac_model_{charge}.svg",
                    title=title,
                )

            delta_ce = ransac.predict(library.obs[["MASS", "PRECURSOR_CHARGE"]])
            library.obs["COLLISION_ENERGY"] = np.maximum(0, library.obs["COLLISION_ENERGY"] + delta_ce)

        else:
            ce_alignment = alignment_library.obs.groupby(by=["COLLISION_ENERGY"])["SPECTRAL_ANGLE"].mean()

            best_ce = ce_alignment.idxmax()
            pl.plot_mean_sa_ce(
                sa_ce_df=ce_alignment.to_frame().reset_index(),
                filename=results_dir / f"{spectra_file.stem}_mean_spectral_angle_ce.svg",
            )
            pl.plot_violin_sa_ce(
                sa_ce_df=alignment_library.obs[["COLLISION_ENERGY", "SPECTRAL_ANGLE"]],
                filename=results_dir / f"{spectra_file.stem}_violin_spectral_angle_ce.svg",
            )
            library.obs["COLLISION_ENERGY"] = best_ce
            with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
                f.write(str(best_ce))
                f.close()
    else:
        best_ce = 35
        library.obs["COLLISION_ENERGY"] = best_ce

        with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
            f.write(str(best_ce))
            f.close()\


def ce_calib(spectra_file: Path, config: Config) -> Spectra:
    ce_calib_step = ProcessStep(config.output, "ce_calib." + spectra_file.stem)
    if ce_calib_step.is_done():
        hdf5_path = config.output / "data" / spectra_file.with_suffix(".mzml.hdf5").name
        if hdf5_path.is_file():
            library = Spectra.from_hdf5(hdf5_path)
            return library
        else:
            raise FileNotFoundError(f"{hdf5_path} not found but ce_calib.{spectra_file.stem} found. Please check.")
    tims_meta_file = None
    if config.spectra_type.lower() in ["hdf", "d"]:  # if it is timstof
        tims_meta_file = config.output / "msms" / spectra_file.with_suffix(".timsmeta").name
    aspec = _annotate_and_get_library(spectra_file, config, tims_meta_file=tims_meta_file)
    _get_best_ce(aspec, spectra_file, config)

    aspec.write_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.hdf5").name)

    ce_calib_step.mark_done()

    return aspec

def run_ce_calibration(
    config_path: Union[str, Path],
):
    """
    Create a CeCalibration object and run the CE calibration.

    # TODO full description
    :param config_path: path to config file
    """
    config = Config()
    config.read(config_path)

    # load spectra file names
    spectra_files = pp.list_spectra(input_dir=config.spectra, input_format=config.spectra_type)

    proc_dir = config.output / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)

    spectra_files = preprocess(spectra_files, config)

    if config.num_threads > 1:
        processing_pool = JobPool(processes=config.num_threads)
        for spectra_file in spectra_files:
            processing_pool.apply_async(ce_calib, [spectra_file, config])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            ce_calib(spectra_file, config)