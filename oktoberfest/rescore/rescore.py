import logging
import subprocess
from pathlib import Path
from typing import List, Optional, Union

import mokapot
import numpy as np
import pandas as pd
from spectrum_fundamentals.metrics.percolator import Percolator

from ..data import FragmentType, Spectra

logger = logging.getLogger(__name__)


def generate_features(
    library: Spectra,
    search_type: str,
    output_file: Union[str, Path],
    all_features: bool = False,
    regression_method: str = "spline",
):
    """
    Generate features to be used for percolator or mokapot target decoy separation.

    The function calculates a range of metrics and features on the provided library for the chosen
    fdr estimation method, then writes the input tab file to the chosen output file.

    :param library: the library to perform feature generation on
    :param search_type: One of "original" and "rescore", which determines the generated features
    :param output_file: the location to the generated tab file to be used for percolator / mokapot
    :param all_features: whether to use all features or only the standard set TODO
    :param regression_method: The regression method to use for iRT alignment

    :Example:

    .. code-block:: python

        from oktoberfest import rescore as re
        from oktoberfest import predict as pr
        import pandas as pd
        import numpy as np
        # Required columns: RAW_FILE, MODIFIED_SEQUENCE, SEQUENCE, CALCULATED_MASS, SCAN_NUMBER,
        # COLLISION_ENERGY, PRECURSOR_CHARGE, REVERSE and SCORE
        meta_df = pd.DataFrame({"RAW_FILE": ["File1","File1"],
                                "MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]RFVQ","RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"],
                                "SEQUENCE": ["AAACRFVQ","RMPCHKPYL"],
                                "CALCULATED_MASS": [1000,4000],
                                "SCAN_NUMBER": [1,2],
                                "COLLISION_ENERGY": [30,35],
                                "PRECURSOR_CHARGE": [1,2],
                                "FRAGMENTATION": ["HCD","HCD"],
                                "REVERSE": [False,False],
                                "SCORE": [0,0]})
        var = Spectra._gen_vars_df()
        library = Spectra(obs=meta_df, var=var)
        raw_intensities = np.random.rand(2,174)
        mzs = np.random.rand(2,174)*1000
        annotation = np.array([var.index,var.index])
        library.add_intensities(raw_intensities, annotation, FragmentType.RAW)
        library.add_mzs(mzs, FragmentType.MZ)
        library.strings_to_categoricals()
        pr.predict_intensities(data=library,
                            model_name="Prosit_2020_intensity_HCD",
                            server_url="koina.wilhelmlab.org:443",
                            ssl=True,
                            targets=["intensities", "annotation"])
        re.generate_features(library, search_type="original", regression_method="spline", output_file="original.tab")
    """
    perc_features = Percolator(
        metadata=library.get_meta_data().reset_index(drop=True),
        pred_intensities=library.get_matrix(FragmentType.PRED)[0],
        true_intensities=library.get_matrix(FragmentType.RAW)[0],
        mz=library.get_matrix(FragmentType.MZ)[0],
        input_type=search_type,
        all_features_flag=all_features,
        regression_method=regression_method,
    )
    perc_features.calc()
    perc_features.write_to_file(str(output_file))


def merge_input(
    tab_files: List[Path],
    output_file=Union[str, Path],
):
    """
    Merge spectra file identifier specific tab files into one large file for combined percolation.

    The function takes a list of tab files and concatenates them before writing a combined tab file back to
    the chosen output file location.

    Fastest solution according to:
    https://stackoverflow.com/questions/44211461/what-is-the-fastest-way-to-combine-100-csv-files-with-headers-into-one

    :param tab_files: list of paths pointing to the individual tab files to be concatenated
    :param output_file: path to the generated output tab file

    :Example:

    .. code-block:: python

        from oktoberfest import rescore as re
        from pathlib import Path
        tabfile1 = Path("out/results/percolator/rescore1.tab")
        tabfile2 = Path("out/results/percolator/rescore2.tab")
        filelist = [tabfile1,tabfile2]
        re.merge_input(filelist, "out/results/percolator/rescore2.tab")
    """
    with open(output_file, "wb") as fout:
        first = True
        for tab_file in tab_files:
            with open(tab_file, "rb") as f:
                if not first:
                    next(f)  # skip the header
                else:
                    first = False
                fout.write(f.read())
            tab_file.unlink()

    # TODO make this more efficient
    df_prosit = pd.read_csv(output_file, sep="\t")
    df_prosit = df_prosit.fillna(0)

    # We exploit the expmass column here by assigning a unique id per filename+scannr group.
    # This ensures percolator will not deduplicate scannrs between filenames, as it uses only
    # filename+ExpMass for TDC.
    df_prosit.insert(loc=4, column="ExpMass", value=df_prosit.groupby(["filename", "ScanNr"]).ngroup())

    df_prosit.to_csv(output_file, sep="\t", index=False)


def rescore_with_percolator(
    input_file: Union[str, Path],
    output_folder: Optional[Union[str, Path]] = None,
    num_threads: int = 3,
    test_fdr: float = 0.01,
    train_fdr: float = 0.01,
):
    """
    Rescore using percolator.

    The function takes an input file location, as well as an output folder location for percolator outputs and
    executes percolator with additional optional parameters.

    :param input_file: Path to percolator tab file
    :param output_folder: An optional output folder for all percolator files, default is the parent directory of the input_file
    :param num_threads: The number of threads used in parallel for percolator
    :param test_fdr: the fdr cutoff for the test set
    :param train_fdr: the fdr cutoff for the train set
    :raises FileNotFoundError: if the input file does not exist

    :Example:

    .. code-block:: python

        from oktoberfest import rescore as re
        from pathlib import Path
        file = Path("out/results/percolator/rescore.tab")
        re.rescore_with_percolator(file, num_threads=3, test_fdr=0.01, train_fdr=0.01)
    """
    if isinstance(input_file, str):
        input_file = Path(input_file)

    if not input_file.is_file():
        raise FileNotFoundError(f"{input_file} does not exist.")

    if output_folder is None:
        output_folder = input_file.parent
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)

    file_prefix = input_file.stem
    weights_file = output_folder / f"{file_prefix}.percolator.weights.csv"
    target_psms = output_folder / f"{file_prefix}.percolator.psms.txt"
    decoy_psms = output_folder / f"{file_prefix}.percolator.decoy.psms.txt"
    target_peptides = output_folder / f"{file_prefix}.percolator.peptides.txt"
    decoy_peptides = output_folder / f"{file_prefix}.percolator.decoy.peptides.txt"
    log_file = output_folder / f"{file_prefix}.log"

    cmd = f"percolator --weights {weights_file} \
                    --num-threads {num_threads} \
                    --subset-max-train 500000 \
                    --post-processing-tdc \
                    --testFDR {test_fdr} \
                    --trainFDR {train_fdr} \
                    --results-psms {target_psms} \
                    --decoy-results-psms {decoy_psms} \
                    --results-peptides {target_peptides} \
                    --decoy-results-peptides {decoy_peptides} \
                    {input_file} 2> {log_file}"

    logger.info(f"Starting percolator with command {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    logger.info("Finished rescoring using percolator.")


def rescore_with_mokapot(
    input_file: Union[str, Path],
    output_folder: Optional[Union[str, Path]] = None,
    test_fdr: float = 0.01,
):
    """
    Rescore using mokapot.

    The function takes an input file location, as well as an output folder location for mokapot outputs and
    executes mokapot with additional optional parameters.

    :param input_file: Path to percolator tab file
    :param output_folder: An optional output folder for all percolator files, default is the parent directory of the input_file
    :param test_fdr: the fdr cutoff for the test set
    :raises FileNotFoundError: if the input file does not exist

    :Example:

    .. code-block:: python

        from oktoberfest import rescore as re
        from pathlib import Path
        file = Path("out/results/percolator/rescore.tab")
        re.rescore_with_mokapot(file, test_fdr=0.01)
    """
    if isinstance(input_file, str):
        input_file = Path(input_file)

    if not input_file.is_file():
        raise FileNotFoundError(f"{input_file} does not exist.")

    if output_folder is None:
        output_folder = input_file.parent
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)

    mokapot_logger = logging.getLogger("mokapot")
    mokapot_logger.setLevel(logging.INFO)
    log_file = output_folder / f"{input_file.stem}.log"

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    log_formatter = logging.Formatter("%(levelname)s: %(message)s")
    file_handler.setFormatter(log_formatter)
    mokapot_logger.addHandler(file_handler)

    np.random.seed(123)

    psms = mokapot.read_pin(input_file)
    logger.info(f"Number of PSMs used for training: {len(psms)}")

    logger.info("Running mokapot...")
    results, models = mokapot.brew(psms, test_fdr=test_fdr)
    results.to_txt(dest_dir=output_folder, file_root=f"{input_file.stem}", decoys=True)
    logger.info("Finished rescoring using mokapot.")
