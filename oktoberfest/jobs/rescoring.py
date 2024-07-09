import logging
from pathlib import Path
from typing import Union

from oktoberfest import plotting as pl
from oktoberfest import predict as pr
from oktoberfest import preprocessing as pp
from oktoberfest import rescore as re
from oktoberfest.jobs import ce_calibration as ce

from ..utils import Config, JobPool, ProcessStep, group_iterator

logger = logging.getLogger(__name__)


def _calculate_features(spectra_file: Path, config: Config):
    library = ce.ce_calib(spectra_file, config)

    calc_feature_step = ProcessStep(config.output, "calculate_features." + spectra_file.stem)
    if calc_feature_step.is_done():
        return

    predict_step = ProcessStep(config.output, "predict." + spectra_file.stem)
    if not predict_step.is_done():

        predict_kwargs = {
            "server_url": config.prediction_server,
            "ssl": config.ssl,
        }

        if "alphapept" in config.models["intensity"].lower():
            chunk_idx = list(group_iterator(df=library.obs, group_by_column="PEPTIDE_LENGTH"))
        else:
            chunk_idx = None
        pr.predict_intensities(
            data=library, chunk_idx=chunk_idx, model_name=config.models["intensity"], **predict_kwargs
        )

        pr.predict_rt(data=library, model_name=config.models["irt"], **predict_kwargs)

        library.write_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.hdf5").name)
        predict_step.mark_done()

    # produce percolator tab files
    fdr_dir = config.output / "results" / config.fdr_estimation_method
    fdr_dir.mkdir(exist_ok=True)

    re.generate_features(
        library=library,
        search_type="original",
        output_file=fdr_dir / spectra_file.with_suffix(".original.tab").name,
        all_features=config.all_features,
        regression_method=config.curve_fitting_method,
    )
    re.generate_features(
        library=library,
        search_type="rescore",
        output_file=fdr_dir / spectra_file.with_suffix(".rescore.tab").name,
        all_features=config.all_features,
        regression_method=config.curve_fitting_method,
    )

    calc_feature_step.mark_done()


def _rescore(fdr_dir: Path, config: Config):
    """
    High level rescore function for original and rescore.

    :param fdr_dir: the output directory
    :param config: the configuration object
    :raises ValueError: if the provided fdr estimation method in the config is not recognized
    """
    rescore_original_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_original")
    rescore_prosit_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_prosit")

    if config.fdr_estimation_method == "percolator":
        if not rescore_original_step.is_done():
            re.rescore_with_percolator(input_file=fdr_dir / "original.tab", output_folder=fdr_dir)
            rescore_original_step.mark_done()
        if not rescore_prosit_step.is_done():
            re.rescore_with_percolator(input_file=fdr_dir / "rescore.tab", output_folder=fdr_dir)
            rescore_prosit_step.mark_done()
    elif config.fdr_estimation_method == "mokapot":
        if not rescore_original_step.is_done():
            re.rescore_with_mokapot(input_file=fdr_dir / "original.tab", output_folder=fdr_dir)
            rescore_original_step.mark_done()
        if not rescore_prosit_step.is_done():
            re.rescore_with_mokapot(input_file=fdr_dir / "rescore.tab", output_folder=fdr_dir)
            rescore_prosit_step.mark_done()
    else:
        raise ValueError(
            'f{config.fdr_estimation_method} is not a valid rescoring tool, use either "percolator" or "mokapot"'
        )


def run_rescoring(config_path: Union[str, Path]):
    """
    Create a ReScore object and run the rescoring.

    # TODO full description
    :param config_path: path to config file
    """
    config = Config()
    config.read(config_path)

    # load spectra file names
    spectra_files = pp.list_spectra(input_dir=config.spectra, input_format=config.spectra_type)

    proc_dir = config.output / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)

    spectra_files = ce.preprocess(spectra_files, config)

    if config.num_threads > 1:
        processing_pool = JobPool(processes=config.num_threads)
        for spectra_file in spectra_files:
            processing_pool.apply_async(_calculate_features, [spectra_file, config])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            _calculate_features(spectra_file, config)

    # prepare rescoring

    fdr_dir = config.output / "results" / config.fdr_estimation_method

    original_tab_files = [fdr_dir / spectra_file.with_suffix(".original.tab").name for spectra_file in spectra_files]
    rescore_tab_files = [fdr_dir / spectra_file.with_suffix(".rescore.tab").name for spectra_file in spectra_files]

    prepare_tab_original_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_prepare_tab_original")
    prepare_tab_rescore_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_prepare_tab_prosit")

    if not prepare_tab_original_step.is_done():
        logger.info("Merging input tab files for rescoring without peptide property prediction")
        re.merge_input(tab_files=original_tab_files, output_file=fdr_dir / "original.tab")
        prepare_tab_original_step.mark_done()

    if not prepare_tab_rescore_step.is_done():
        logger.info("Merging input tab files for rescoring with peptide property prediction")
        re.merge_input(tab_files=rescore_tab_files, output_file=fdr_dir / "rescore.tab")
        prepare_tab_rescore_step.mark_done()

    # rescoring
    _rescore(fdr_dir, config)

    # plotting
    logger.info("Generating summary plots...")
    pl.plot_all(fdr_dir)

    logger.info("Finished rescoring.")