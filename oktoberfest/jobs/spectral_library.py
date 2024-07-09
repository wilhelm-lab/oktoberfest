import logging
import pickle
import sys
import time
from functools import partial
from math import ceil
from multiprocessing import Manager, Process, pool
from pathlib import Path
from typing import Dict, List, Tuple, Type, Union

import pandas as pd
from spectrum_io.spectral_library import MSP, DLib, SpectralLibrary, Spectronaut
from tqdm.auto import tqdm

from oktoberfest import __copyright__, __version__
from oktoberfest import predict as pr
from oktoberfest import preprocessing as pp

from ..data.spectra import Spectra
from ..utils import Config, ProcessStep, group_iterator

logger = logging.getLogger(__name__)


def _make_predictions_error_callback(failure_progress_tracker, failure_lock, error):
    logger.error(
        f"Prediction failed due to: {error} Batch will be missing from output file. "
        "DO NOT STOP THIS RUN: The index of the batch is stored and your output file will be appended "
        "by the missing batches if you rerun without changing your config file after this run is completed."
    )
    with failure_lock:
        failure_progress_tracker.value += 1

def _make_predictions(int_model, irt_model, predict_kwargs, queue_out, progress, lock, batch_df):
    predictions = {
        **pr.predict_at_once(batch_df, model_name=int_model, **predict_kwargs),
        **pr.predict_at_once(batch_df, model_name=irt_model, **predict_kwargs),
    }
    queue_out.put((predictions, batch_df))
    with lock:
        progress.value += 1


def _get_writer_and_output(results_path: Path, output_format: str) -> Tuple[Type[SpectralLibrary], Path]:
    if output_format == "msp":
        return MSP, results_path / "myPrositLib.msp"
    elif output_format == "spectronaut":
        return Spectronaut, results_path / "myPrositLib.csv"
    elif output_format == "dlib":
        return DLib, results_path / "myPrositLib.dlib"
    else:
        raise ValueError(f"{output_format} is not supported as spectral library type")


def _get_batches_and_mode(out_file: Path, failed_batch_file: Path, obs: pd.DataFrame, batchsize: int, model: str):
    if out_file.is_file():
        if failed_batch_file.is_file():
            with open(failed_batch_file, "rb") as fh:
                batch_iterator = pickle.load(fh)
            mode = "a"
            logger.warning(
                f"Found existing spectral library {out_file}. "
                "Attempting to append missing batches from previous run..."
            )
        else:
            logger.error(
                f"A file {out_file} already exists but no information about missing batches "
                "from a previous run could be found. Stopping to prevent corruption / data loss. "
                "If this is intended, delete the file and rerun."
            )
            sys.exit(1)
    else:
        if "alphapept" in model.lower():
            batch_iterator = group_iterator(df=obs, group_by_column="PEPTIDE_LENGTH", max_batch_size=batchsize)
        else:
            batch_iterator = (
                obs.index[i * batchsize : (i + 1) * batchsize].to_numpy() for i in range(ceil(len(obs) / batchsize))
            )
        mode = "w"

    return list(batch_iterator), mode


def _speclib_from_digestion(config: Config) -> Spectra:
    library_input_type = config.library_input_type
    if library_input_type == "fasta":
        digest_step = ProcessStep(config.output, "speclib_digested")
        library_file = config.output / "prosit_input.csv"
        if not digest_step.is_done():
            peptide_dict = pp.digest(
                fasta=config.library_input,
                digestion=config.digestion,
                missed_cleavages=config.missed_cleavages,
                db=config.db,
                enzyme=config.enzyme,
                special_aas=config.special_aas,
                min_length=config.min_length,
                max_length=config.max_length,
            )
            metadata = pp.generate_metadata(
                peptides=list(peptide_dict.keys()),
                collision_energy=config.collision_energy,
                precursor_charge=config.precursor_charge,
                fragmentation=config.fragmentation,
                instrument_type=config.instrument_type,
                proteins=list(peptide_dict.values()),
            )
            library_file = config.output / "prosit_input.csv"
            metadata.to_csv(library_file, sep=",", index=None)
            digest_step.mark_done()
    elif library_input_type == "peptides":
        library_file = config.library_input
    else:
        raise ValueError(f'Library input type {library_input_type} not understood. Can only be "fasta" or "peptides".')
    spec_library = pp.gen_lib(library_file)

    pp_and_filter_step = ProcessStep(config.output, "speclib_filtered")

    data_dir = config.output / "data"
    if not pp_and_filter_step.is_done():
        data_dir.mkdir(exist_ok=True)
        spec_library = pp.process_and_filter_spectra_data(
            library=spec_library, model=config.models["intensity"], tmt_label=config.tag
        )
        spec_library.write_as_hdf5(data_dir / f"{library_file.stem}_filtered.hdf5")
        pp_and_filter_step.mark_done()
    else:
        spec_library = Spectra.from_hdf5(data_dir / f"{library_file.stem}_filtered.hdf5")

    return spec_library


def _update(pbar: tqdm, postfix_values: Dict[str, int]):
    total_val = sum(postfix_values.values())
    if total_val > pbar.n:
        pbar.set_postfix(**postfix_values)
        pbar.n = total_val
        pbar.refresh()


def _check_write_failed_batch_file(failed_batch_file: Path, n_failed: int, results: List[pool.AsyncResult]):
    if n_failed > 0:
        failed_batches = []
        for i, result in enumerate(results):
            try:
                result.get()
            except Exception:
                failed_batches.append(i)
        logger.error(
            f"Prediction for {n_failed} / {i+1} batches failed. Check the log to find out why. "
            "Then rerun without changing the config file to append only the missing batches to your output file."
        )
        with open(failed_batch_file, "wb") as fh:
            pickle.dump(failed_batches, fh)
        sys.exit(1)


def generate_spectral_lib(config_path: Union[str, Path]):
    """
    Create a SpectralLibrary object and generate the spectral library.

    # TODO full description
    :param config_path: path to config file
    """
    config = Config()
    config.read(config_path)

    spec_library = _speclib_from_digestion(config)

    server_kwargs = {
        "server_url": config.prediction_server,
        "ssl": config.ssl,
        "disable_progress_bar": True,
    }

    speclib_written_step = ProcessStep(config.output, "speclib_written")
    if not speclib_written_step.is_done():
        results_path = config.output / "results"
        results_path.mkdir(exist_ok=True)

        batchsize = config.batch_size
        failed_batch_file = config.output / "data" / "speclib_failed_batches.pkl"
        writer, out_file = _get_writer_and_output(results_path, config.output_format)
        batches, mode = _get_batches_and_mode(
            out_file, failed_batch_file, spec_library.obs, batchsize, config.models["intensity"]
        )
        speclib = writer(out_file, mode=mode, min_intensity_threshold=config.min_intensity)
        n_batches = len(batches)

        with Manager() as manager:
            # setup
            shared_queue = manager.Queue(maxsize=config.num_threads)
            prediction_progress = manager.Value("i", 0)
            prediction_failure_progress = manager.Value("i", 0)
            writing_progress = manager.Value("i", 0)

            lock = manager.Lock()
            lock_failure = manager.Lock()

            # Create a pool for producer processes
            predictor_pool = pool.Pool(config.num_threads)

            consumer_process = Process(
                target=speclib.async_write,
                args=(
                    shared_queue,
                    writing_progress,
                ),
            )

            try:
                results = []
                for batch in batches:
                    result = predictor_pool.apply_async(
                        _make_predictions,
                        (
                            config.models["intensity"],
                            config.models["irt"],
                            server_kwargs,
                            shared_queue,
                            prediction_progress,
                            lock,
                            spec_library.obs.loc[batch],
                        ),
                        error_callback=partial(
                            _make_predictions_error_callback, prediction_failure_progress, lock_failure
                        ),
                    )
                    results.append(result)
                predictor_pool.close()

                with tqdm(
                    total=n_batches, desc="Writing library", postfix={"successful": 0, "missing": 0}
                ) as writer_pbar:
                    consumer_process.start()
                    with tqdm(
                        total=n_batches, desc="Getting predictions", postfix={"successful": 0, "failed": 0}
                    ) as predictor_pbar:
                        while predictor_pbar.n < n_batches:
                            time.sleep(1)
                            pr_fail_val = prediction_failure_progress.value
                            _update(predictor_pbar, {"failed": pr_fail_val, "successful": prediction_progress.value})
                            _update(writer_pbar, {"successful": writing_progress.value, "missing": pr_fail_val})
                        shared_queue.put(None)  # signal the writer process, that it is done
                        predictor_pool.join()  # properly await and terminate the pool

                    while writer_pbar.n < n_batches:  # have to keep updating the writer pbar
                        time.sleep(1)
                        _update(
                            writer_pbar,
                            {"successful": writing_progress.value, "missing": prediction_failure_progress.value},
                        )
                    consumer_process.join()  # properly await the termination of the writer process

                _check_write_failed_batch_file(failed_batch_file, prediction_failure_progress.value, results)

            finally:
                predictor_pool.terminate()
                predictor_pool.join()
                consumer_process.terminate()
                consumer_process.join()

        logger.info("Finished writing the library to disk")
        failed_batch_file.unlink(missing_ok=True)
        speclib_written_step.mark_done()