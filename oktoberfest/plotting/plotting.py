import os
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from scipy import stats
from spectrum_io.raw import ThermoRaw

from oktoberfest.data.spectra import Spectra
from oktoberfest.utils import Config

# Set the default fontsize and linewidth
plt.rcParams.update({"font.size": 14, "axes.linewidth": 1.5, "xtick.major.width": 1.5, "ytick.major.width": 1.5})


def _check_columns(df: pd.DataFrame):
    """
    Check columns to make plotting work for mokapot.

    :param df: a pandas dataframe to check columns for
    :raises AssertionError: if expected columns do not correspond to either mokapot or percolator output

    :return: Tuple of column identifiers: (score, q-value, Proteins, SpecId)
    """
    mokapot_mapping = {
        "mokapot score": "score",
        "mokapot q-value": "q-value",
        "Peptide": "peptide",
        "SpecId": "PSMId",
    }
    if set(mokapot_mapping.keys()).issubset(df.columns):
        column_names = list(mokapot_mapping.keys())
    elif set(mokapot_mapping.values()).issubset(df.columns):
        column_names = list(mokapot_mapping.values())
    else:
        raise AssertionError("Missing columns in one of the input files. Please check.")

    return column_names


def plot_score_distribution(target: pd.DataFrame, decoy: pd.DataFrame, level: str, filename: Union[str, Path]):
    """
    Generate histogram of the score distribution for targets and decoys.

    :param target: mokapot / percolator target output
    :param decoy: mokapot / percolator decoy output
    :param level: The level on which to produce the comparison. Can be either "peptide" or "psm"
    :param filename: the path to the location used for storing the plot

    :Example:

    .. code-block:: python

        >>> from oktoberfest import plotting as pl
        >>> import pandas as pd
        >>> # Required columns: PSMId, score, q-value and peptide
        >>> target_df = pd.DataFrame({"PSMId": ["F1-15-TAIASPEK-1-5","F1-5-HARPQTTLR-2-6","F2-14-RVYDPASPQRR-2-5",
        >>>                             "F1-12-FSTQDHAAAAIAK-2-2","F2-63-ISDPTSPLRTR-2-9","F1-16-ADHPLRTR-1-5"],
        >>>                             "score": [-0.1,-0.5,-0.5,0.7,0.4,0.7],
        >>>                             "q-value": [0.005,0.008,0.002,0.006,0.004,0.001],
        >>>                             "peptide": ["TAIASPEK","HARPQTTLR","RVYDPASPQRR",
        >>>                             "FSTQDHAAAAIAK","ISDPTSPLRTR","ADHPLRTR"]})
        >>> decoy_df = pd.DataFrame({"PSMId": ["F1-11-KLYNANYIK-3-7","F2-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4"],
        >>>                         "score": [-0.1,-0.5,-0.5],
        >>>                         "q-value": [0.006,0.004,0.003],
        >>>                         "peptide": ["KLYNANYIK","LGLTKLQLH","EFAVEVLK"]})
        >>> pl.plot_score_distribution(target=target_df,
        >>>                             decoy=decoy_df,
        >>>                             level="psm",
        >>>                             filename="./tests/doctests/output/score_distribution_plot.svg")
    """
    score_col, _, _, _ = _check_columns(target)

    plt.figure(figsize=(8, 6))
    bins = np.linspace(-3, 2, 15).tolist()
    plt.hist(target[score_col], bins, label="Targets", rwidth=0.5, color="#48AF00", alpha=1.0)
    plt.hist(decoy[score_col], bins, label="Decoys", rwidth=0.5, color="#FE7312", alpha=0.7)
    plt.xlabel("Score")
    plt.legend(loc="upper right")
    plt.title(f"Score Distribution ({level.capitalize()})")
    plt.savefig(filename, dpi=300)
    plt.close()


def joint_plot(
    prosit_target: pd.DataFrame,
    prosit_decoy: pd.DataFrame,
    andromeda_target: pd.DataFrame,
    andromeda_decoy: pd.DataFrame,
    level: str,
    filename: Union[str, Path],
):
    """
    Generate joint plot to compare rescoring with and without peptide property predictions.

    :param prosit_target: mokapot / percolator target output for rescoring with peptide property prediction
    :param prosit_decoy: mokapot / percolator decoy output for rescoring with peptide property prediction
    :param andromeda_target: mokapot / percolator target output for rescoring without peptide property prediction
    :param andromeda_decoy: mokapot / percolator decoy output for rescoring without peptide property prediction
    :param level: The level on which to produce the comparison. Can be either "peptide" or "psm"
    :param filename: the path to the location used for storing the plot

    :raises ValueError: if a wrong level is provided

    :Example:

    .. code-block:: python

        >>> from oktoberfest import plotting as pl
        >>> import pandas as pd
        >>> # Required columns: PSMId, score, q-value and peptide
        >>> target_df = pd.DataFrame({"PSMId": ["F1-15-TAIASPEK-1-5","F2-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4",
        >>>                         "F2-63-ISDPTSPLRTR-2-9","F1-16-ADHPLRTR-1-5","F2-4-YLNPLRTK-1-5"],
        >>>                         "q-value": [0.005,0.008,0.002,0.006,0.004,0.001],
        >>>                         "score": [-0.1,-0.5,-0.5,0.7,0.4,0.5],
        >>>                         "peptide": ["TAIASPEK","LGLTKLQLH","EFAVEVLK","ISDPTSPLRTSR","ADHPLRTR","YLNPLRTK"]})
        >>> decoy_df = pd.DataFrame({"PSMId": ["F1-11-KLYNANYIK-3-7","F2-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4"],
        >>>                         "q-value": [0.006,0.004,0.003],
        >>>                         "score": [-0.1,-0.5,-0.5],
        >>>                         "peptide": ["KLYNANYIK","LGLTKLQLH","EFAVEVLK"]})
        >>> andromeda_target_df = pd.DataFrame({"PSMId": ["F1-15-TAIASPEK-1-5","F2-59-LGLTKLQLH-3-9","F2-63-ISDPTSPLRTR-2-9",
        >>>                                     "F1-16-ADHPLRTR-1-5","F2-4-YLNPLRTK-1-5"],
        >>>                                     "q-value": [0.005,0.008,0.002,0.006,0.001],
        >>>                                     "score": [-0.2,-0.8,0.5,0.3,0.4],
        >>>                                     "peptide": ["TAIASPEK","LGLTKLQLH","ISDPTSPLRTSR","ADHPLRTR","YLNPLRTK"]})
        >>> andromeda_decoy_df = pd.DataFrame({"PSMId": ["F1-11-KLYNANYIK-3-7","F2-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4"],
        >>>                                 "q-value": [0.007,0.005,0.002],
        >>>                                 "score": [-0.2,-0.7,-0.8],
        >>>                                 "peptide": ["KLYNANYIK","LGLTKLQLH","EFAVEVLK"]})
        >>> pl.joint_plot(prosit_target=target_df,
        >>>                 prosit_decoy=decoy_df,
        >>>                 andromeda_target=andromeda_target_df,
        >>>                 andromeda_decoy=andromeda_decoy_df,
        >>>                 level="psm",
        >>>                 filename="./tests/doctests/output/joint_plot.svg")
    """
    score_col, _, peptide_col, psm_col = _check_columns(prosit_target)
    if level.lower() == "peptide":
        join_col = peptide_col
    elif level.lower() == "psm":
        join_col = psm_col
    else:
        raise ValueError(f"level can only be peptide or psm. Given: {level}")

    targets = prosit_target.merge(
        andromeda_target, on=join_col, how="outer", suffixes=["_prosit", "_andromeda"], indicator=True
    )
    decoys = prosit_decoy.merge(
        andromeda_decoy, on=join_col, how="outer", suffixes=["_prosit", "_andromeda"], indicator=True
    )
    df_targets = pd.DataFrame()
    df_targets["prosit_score"] = targets[f"{score_col}_prosit"]
    df_targets["and_score"] = targets[f"{score_col}_andromeda"]
    df_targets["type"] = "targets"
    df_targets["color"] = "#48AF00"

    df_decoys = pd.DataFrame()
    df_decoys["prosit_score"] = decoys[f"{score_col}_prosit"]
    df_decoys["and_score"] = decoys[f"{score_col}_andromeda"]
    df_decoys["type"] = "decoys"
    df_decoys["color"] = "#FE7312"

    df_all = pd.concat([df_targets, df_decoys], axis=0)
    df_all.dropna(inplace=True)
    jplot = sns.jointplot(
        data=df_all,
        x="and_score",
        y="prosit_score",
        marginal_ticks=True,
        hue=df_all["type"],
        palette=["#15853B", "#E2700E"],
        ratio=2,
        height=10,
        joint_kws={"rasterized": True, "edgecolor": "none", "s": 10},
    )
    jplot.ax_joint.axhline(y=0, c="red")
    jplot.ax_joint.axvline(x=0, c="red")
    jplot.ax_marg_y.axhline(y=0, c="red")
    jplot.ax_marg_x.axvline(x=0, c="red")

    jplot.ax_joint.set_ylabel("Score\n(peptide property prediction)")
    jplot.ax_joint.set_xlabel("Score\n(search engine)")
    jplot.figure.suptitle(f"Score distribution ({level.capitalize()})", y=0.99)
    plt.savefig(filename, dpi=300)
    plt.close()


def plot_gain_loss(prosit_target: pd.DataFrame, andromeda_target: pd.DataFrame, level: str, filename: Union[str, Path]):
    """
    Generate venn barplots to show lost, common and shared targets below 1% FDR attributed to peptide property predictions.

    :param prosit_target: mokapot / percolator target output for rescoring with peptide property prediction
    :param andromeda_target: mokapot / percolator target output for rescoring without peptide property prediction
    :param level: The level on which to produce the comparison. Can be either "peptide" or "psm"
    :param filename: the path to the location used for storing the plot

    :raises ValueError: if a wrong level is provided

    :Example:

    .. code-block:: python

        >>> from oktoberfest import plotting as pl
        >>> import pandas as pd
        >>> # Required columns: PSMId, score, q-value and peptide
        >>> prosit_df = pd.DataFrame({"PSMId": ["F1-15-TAIASPEK-1-5","F2-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4",
        >>>                         "F2-63-ISDPTSPLRTR-2-9","F1-16-ADHPLRTR-1-5"],
        >>>                         "q-value": [0.005,0.008,0.002,0.006,0.004],
        >>>                         "score": [-0.1,-0.5,-0.5,0.7,0.4],
        >>>                         "peptide": ["TAIASPEK","LGLTKLQLH","EFAVEVLK","ISDPTSPLRTSR","ADHPLRTR"]})
        >>> andromeda_df = pd.DataFrame({"PSMId": ["F1-11-KLYNANYIK-3-7","F2-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4"],
        >>>                             "q-value": [0.006,0.004,0.003],
        >>>                             "score": [-0.1,-0.5,-0.5],
        >>>                             "peptide": ["KLYNANYIK","LGLTKLQLH","EFAVEVLK"]})
        >>> pl.plot_gain_loss(prosit_target=prosit_df,
        >>>                     andromeda_target=andromeda_df,
        >>>                     level="psm",
        >>>                     filename="./tests/doctests/output/gain_loss_psm_plot.svg")
    """
    _, qval_col, peptide_col, psm_col = _check_columns(prosit_target)

    if level.lower() == "peptide":
        join_col = peptide_col
    elif level.lower() == "psm":
        join_col = psm_col
    else:
        raise ValueError(f"level can only be peptide or psm. Given: {level}")

    andromeda_target = andromeda_target[andromeda_target[qval_col] < 0.01]
    prosit_target = prosit_target[prosit_target[qval_col] < 0.01]
    merged_df = prosit_target.merge(andromeda_target, how="inner", on=join_col)

    shared = len(merged_df.index)
    gained = len(prosit_target.index) - shared
    lost = len(andromeda_target.index) - shared

    fig, ax = plt.subplots(1, figsize=(1.5, 10))
    labels = [""]
    ax1 = ax.bar(labels, shared, width=0.5, color="#115795")
    ax2 = ax.bar(labels, gained, width=0.5, bottom=shared, color="#007D3E")
    ax3 = ax.bar(labels, -lost, color="#E17224", width=0.5)

    for r1, r2, r3, v1, v2, v3 in zip(ax1, ax2, ax3, [shared], [gained], [lost]):
        h1 = r1.get_height()
        h2 = r2.get_height()
        plt.text(
            r1.get_x() + r1.get_width() / 2.0,
            h1 / 2.0,
            "%d" % v1,
            ha="center",
            va="bottom",
            color="white",
            fontweight="bold",
        )
        plt.text(
            r2.get_x() + r2.get_width() / 2.0,
            h1 + h2 / 2,
            "%d" % v2,
            ha="center",
            va="bottom",
            color="black",
        )
        plt.text(
            r3.get_x() + r3.get_width() / 2.0,
            -0.025 * gained,
            "%d" % -v3,
            ha="center",
            va="bottom",
            color="black",
        )

    plt.ylim(-lost - 100, h1 + h2 + 30)
    # remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    # grid
    ax.set_ylabel(f"number of target {level.lower()}s below 1% FDR")
    ax.set_axisbelow(True)
    ax.yaxis.grid(color="black")
    ax.tick_params(axis="y", which="major")
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)

    legend_label = ["Common", "Gained", "Lost"]
    plt.legend(legend_label, ncol=1, bbox_to_anchor=([1.2, 0.5, 0, 0]), frameon=False)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def plot_mean_sa_ce(sa_ce_df: pd.DataFrame, filename: Union[str, Path]):
    """
    Generate dotplot for spectral angle distribution over range of collision energies used for fragment intensity prediction.

    :param sa_ce_df: a dataframe containing the two columns "COLLISION_ENERGY", "SPECTRAL_ANGLE".
    :param filename: the path to the location used for storing the plot

    :Example:

    .. code-block:: python

        >>> from oktoberfest import plotting as pl
        >>> import pandas as pd
        >>> # Required columns: SPECTRA_ANGLE and COLLISION_ENERGY
        >>> sa_ce_df = pd.DataFrame({"SPECTRAL_ANGLE": [0.7,0.5,0.3,0.8], "COLLISION_ENERGY": [34,31,34,31]})
        >>> pl.plot_mean_sa_ce(sa_ce_df=sa_ce_df, filename="./tests/doctests/output/mean_sa_ce_plot.svg")
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    sns.scatterplot(data=sa_ce_df, x="COLLISION_ENERGY", y="SPECTRAL_ANGLE", ax=ax)
    ax.axvline(x=sa_ce_df["COLLISION_ENERGY"][sa_ce_df["SPECTRAL_ANGLE"].idxmax()], color="red")
    plt.grid()
    plt.savefig(filename, dpi=300)
    plt.close()


def plot_violin_sa_ce(sa_ce_df: pd.DataFrame, filename: Union[str, Path]):
    """
    Generate violinplot for spectral angle distribution over range of collision energies used for fragment intensity prediction.

    :param sa_ce_df: a dataframe containing the two columns "COLLISION_ENERGY", "SPECTRAL_ANGLE".
    :param filename: the path to the location used for storing the plot

    :Example:

    .. code-block:: python

        >>> from oktoberfest import plotting as pl
        >>> import pandas as pd
        >>> # Required columns: SPECTRA_ANGLE and COLLISION_ENERGY
        >>> sa_ce_df = pd.DataFrame({"SPECTRAL_ANGLE": [0.7,0.5,0.3,0.8], "COLLISION_ENERGY": [34,31,34,31]})
        >>> pl.plot_violin_sa_ce(sa_ce_df=sa_ce_df, filename="./tests/doctests/output/violin_sa_ce_plot.svg")
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    sns.violinplot(data=sa_ce_df, x="COLLISION_ENERGY", y="SPECTRAL_ANGLE", ax=ax, color="#1f77b4")
    ax.axvline(
        x=sa_ce_df["COLLISION_ENERGY"][sa_ce_df["SPECTRAL_ANGLE"].idxmax()] - sa_ce_df["COLLISION_ENERGY"].min(),
        color="red",
    )
    plt.xticks(rotation=90)
    plt.grid()
    plt.savefig(filename, dpi=300)
    plt.close()


def plot_pred_rt_vs_irt(
    prosit_df: pd.DataFrame, prosit_target: pd.DataFrame, outpath: Union[str, Path], suffix: Union[str, Path]
):
    """
    Generate scatterplot to compare predicted indexed retention time against (aligned) experimentally observed retention time.

    :param prosit_df: mokapot / percolator input tab for rescoring with peptide property prediction
    :param prosit_target: mokapot / percolator target output for rescoring with peptide property prediction
    :param outpath: the path to the location used for storing the plot without the filename
    :param suffix: the suffix of the filename, which will be prepended by the rawfile name

    :Example:

    .. code-block:: python

        >>> from oktoberfest import plotting as pl
        >>> import pandas as pd
        >>> # Required columns: SpecId, RT, iRT and pred_RT
        >>> prosit_df = pd.DataFrame({"SpecId": ["F1-15-TAIASPEK-1-5","F1-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4"],
        >>>                         "RT": np.array([28.1,56.54,83.7]),
        >>>                         "iRT": np.array([25.21,54.76,82.88]),
        >>>                         "pred_RT": np.array([32.23,60.01,99.11])})
        >>> # Required columns: PSMId, score, q-value and peptide
        >>> target_df = pd.DataFrame({"PSMId": ["F1-15-TAIASPEK-1-5","F1-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4"],
        >>>                         "q-value": [0.005,0.003,0.002],
        >>>                         "score": [0.7,0.4,0.5],
        >>>                         "peptide": ["TAIASPEK","LGLTKLQLH","EFAVEVLK"]})
        >>> pl.plot_pred_rt_vs_irt(prosit_df=prosit_df,
        >>>                         prosit_target=target_df,
        >>>                         outpath="./tests/doctests/output/",
        >>>                         suffix="pred_irt_vs_irt")
    """
    _, qval_col, _, psm_col = _check_columns(prosit_target)

    targets = prosit_df.merge(prosit_target, how="inner", left_on="SpecId", right_on=psm_col)
    targets = targets.loc[targets[qval_col] < 0.01, ["SpecId", "RT", "iRT", "pred_RT"]]
    targets.columns = ["SpecId", "experimental RT", "aligned RT", "predicted iRT"]
    targets["rawfile"] = targets["SpecId"].str.split("-", n=1).str[0]
    targets.sort_values("predicted iRT")
    for rawfile in targets["rawfile"].unique():
        fig, ax = plt.subplots(figsize=(8, 8))
        data = targets[targets["rawfile"] == rawfile]

        if (data["predicted iRT"] == data["experimental RT"]).all():
            # this can now happen if we skip the iRT prediction as we then use a perfect prediction as replacement
            kernel = np.zeros(data.shape[0])
        else:
            both = np.vstack([data["predicted iRT"], data["experimental RT"]])
            kernel = stats.gaussian_kde(both)(both)

        sns.scatterplot(
            data=data, x="predicted iRT", y="experimental RT", label="predicted iRT", c=kernel, ax=ax, linewidth=0
        )
        sns.lineplot(data=data, x="predicted iRT", y="aligned RT", label="alignment", ax=ax, c="k")
        ax.set_ylabel("(aligned) experimental RT")
        ax.set_xlabel("predicted iRT")
        plt.title("RT alignment")
        plt.legend(labels=("predicted iRT", "alignment"), loc="best", fancybox=True, shadow=True)
        plt.grid()
        plt.savefig(Path(outpath) / f"{rawfile}_{suffix}", dpi=300)
        plt.close()


def plot_sa_distribution(prosit_df: pd.DataFrame, target_df: pd.DataFrame, decoy_df: pd.DataFrame, filename: Path):
    """Generate spectral angle distribution for targets and decoys.

    :param prosit_df: mokapot / percolator input tab for rescoring with peptide property prediction
    :param target_df: mokapot / percolator target output for rescoring with peptide property prediction on the psm level
    :param decoy_df: mokapot / percolator decoy output for rescoring with peptide property prediction on the psm level
    :param filename: the path to the location used for storing the plot

    :Example:

    .. code-block:: python

        >>> from oktoberfest import plotting as pl
        >>> import pandas as pd
        >>> # Required columns: SpecId and spectral_angle
        >>> prosit_df = pd.DataFrame({"SpecId": ["F1-15-TAIASPEK-1-5","F2-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4",
        >>>                         "F2-63-ISDPTSPLRTR-2-9","F1-16-ADHPLRTR-1-5","F1-11-KLYNANYIK-3-7","F2-4-YLNPLRTK-1-5"],
        >>>                         "spectral_angle": [0.6,0.2,0.3,0.6,0.4,0.2,0.5]})
        >>> # Required columns: PSMId, score, q-value and peptide
        >>> target_df = pd.DataFrame({"PSMId": ["F1-15-TAIASPEK-1-5","F2-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4",
        >>>                         "F2-63-ISDPTSPLRTR-2-9","F1-16-ADHPLRTR-1-5","F2-4-YLNPLRTK-1-5"],
        >>>                         "q-value": [0.005,0.008,0.002,0.006,0.004,0.001],
        >>>                         "score": [-0.1,-0.5,-0.5,0.7,0.4,0.5],
        >>>                         "peptide": ["TAIASPEK","LGLTKLQLH","EFAVEVLK","ISDPTSPLRTSR","ADHPLRTR","YLNPLRTK"]})
        >>> decoy_df = pd.DataFrame({"PSMId": ["F1-11-KLYNANYIK-3-7","F2-59-LGLTKLQLH-3-9","F1-24-EFAVEVLK-2-4"],
        >>>                         "q-value": [0.006,0.004,0.003],
        >>>                         "score": [-0.1,-0.5,-0.5],
        >>>                         "peptide": ["KLYNANYIK","LGLTKLQLH","EFAVEVLK"]})
        >>> pl.plot_sa_distribution(prosit_df=prosit_df,
        >>>                         target_df=target_df,
        >>>                         decoy_df=decoy_df,
        >>>                         filename="./tests/doctests/output/sa_distribution_plot.svg")
    """
    _, _, _, psm_col = _check_columns(target_df)
    target = prosit_df.merge(target_df, how="inner", left_on="SpecId", right_on=psm_col)
    decoy = prosit_df.merge(decoy_df, how="inner", left_on="SpecId", right_on=psm_col)
    plt.figure(figsize=(8, 6))
    bins = np.linspace(0, 1, 15).tolist()
    plt.hist(target.spectral_angle, bins, label="Targets", rwidth=0.5, color="#48AF00", alpha=1.0)
    plt.hist(decoy.spectral_angle, bins, label="Decoys", rwidth=0.5, color="#FE7312", alpha=0.7)
    plt.xlabel("Spectral angle", size=14)
    plt.title("Target vs Decoys Spectral Angle Distribution")
    plt.legend(loc="upper right")
    plt.savefig(filename, dpi=300)
    plt.close()


def plot_mirror_spectrum(
    spec_pred: Spectra,
    mzml: pd.DataFrame,
    raw_file: str,
    scan_number: int,
    config: Config,
    prosit_df: pd.DataFrame,
    target_df: pd.DataFrame,
    decoy_df: pd.DataFrame,
    pdf: PdfPages,
):
    """
    Generate a mirror plot comparing an experimental and predicted MS/MS spectrum.

    :param spec_pred: Spectra object containing predicted MS/MS spectra
    :param mzml: ThermoRaw object containing experimental MS/MS spectra from an mzML file
    :param raw_file: The name of the raw file being processed
    :param scan_number: The scan number of the spectrum to be plotted
    :param config: the configuration object
    :param prosit_df: mokapot / percolator input tab for rescoring with peptide property prediction
    :param target_df: mokapot / percolator target output for rescoring with peptide property prediction on the psm level
    :param decoy_df: mokapot / percolator decoy output for rescoring with peptide property prediction on the psm level
    :param pdf: PDF file object for saving mirror plots

    :raises ValueError: If the mass analyzer type is unknown.
    """
    score_col, _, _, spec_col = _check_columns(target_df)
    target_df["target"] = True
    decoy_df["target"] = False
    concat_target_decoy = pd.concat([target_df, decoy_df])
    if spec_col == "PSMId":
        concat_target_decoy["ScanNr"] = concat_target_decoy["PSMId"].str.split("-").str[-1].astype(int)
    filtered_obs = spec_pred.obs[spec_pred.obs["SCAN_NUMBER"] == scan_number]
    if filtered_obs.empty:
        print(f"Warning: Scan number {scan_number} for {raw_file} not found in the prediction file.")
        return
    obs = filtered_obs.iloc[0]

    mod_sequence = obs["MODIFIED_SEQUENCE"].replace("[]-", "").replace("-[]", "")
    charge = obs["PRECURSOR_CHARGE"]
    mass = obs["MASS"]
    mass_analyzer = obs["MASS_ANALYZER"]
    fragm = obs["FRAGMENTATION"]
    raw_file = obs["RAW_FILE"]
    rt = obs["RETENTION_TIME"].round(2)
    ce = obs["COLLISION_ENERGY"]
    model = config.models["intensity"]
    ion_types = config.ion_types
    abs_rt_diff = prosit_df[(prosit_df["ScanNr"] == scan_number) & (prosit_df["filename"] == raw_file)][
        "abs_rt_diff"
    ].iloc[0]
    sa = prosit_df[(prosit_df["ScanNr"] == scan_number) & (prosit_df["filename"] == raw_file)]["spectral_angle"].iloc[0]
    ce = prosit_df[(prosit_df["ScanNr"] == scan_number) & (prosit_df["filename"] == raw_file)][
        "collision_energy_aligned"
    ].iloc[0]
    score = concat_target_decoy[
        (concat_target_decoy["ScanNr"] == scan_number) & (concat_target_decoy["filename"] == raw_file)
    ][score_col].iloc[0]

    # Set tolerance based on mass analyzer
    if mass_analyzer == "FTMS":
        fragment_tol_mass = 20.0
        fragment_tol_mode = "ppm"
    elif mass_analyzer == "ITMS":
        fragment_tol_mass = 0.4
        fragment_tol_mode = "Da"
    else:
        raise ValueError(f"Unknown mass analyzer: {mass_analyzer}")

    # Get experimental spectrum
    mz_exp = np.array(mzml[mzml["SCAN_NUMBER"] == scan_number]["MZ"].iloc[0])
    intensity_exp = np.array(mzml[mzml["SCAN_NUMBER"] == scan_number]["INTENSITIES"].iloc[0])

    top_spectrum = sus.MsmsSpectrum("", mass, charge, mz=mz_exp, intensity=intensity_exp)
    top_spectrum = top_spectrum.annotate_proforma(
        mod_sequence, fragment_tol_mass, fragment_tol_mode, ion_types="byrI", max_ion_charge=charge
    )

    # Get predicted spectrum
    idx = spec_pred.obs.index.get_loc(spec_pred.obs[spec_pred.obs["SCAN_NUMBER"] == scan_number].index[0])
    mz_pred = spec_pred.layers["mz"].toarray()[idx]
    intensity_pred = spec_pred.layers["pred_int"].toarray()[idx]

    bot_spectrum = sus.MsmsSpectrum("", mass, charge, mz=mz_pred, intensity=intensity_pred)
    bot_spectrum = bot_spectrum.annotate_proforma(
        mod_sequence, fragment_tol_mass, fragment_tol_mode, ion_types=ion_types, max_ion_charge=charge
    )
    fig = plt.figure(figsize=(12, 7))  # Adjust total figure size

    gs = GridSpec(2, 2, width_ratios=[2, 1], wspace=0.2)

    ax_mirror = fig.add_subplot(gs[:, 0])
    ax_kde = fig.add_subplot(gs[3])
    ax_rt = fig.add_subplot(gs[1])

    # Mirror spectrum plot
    ax_mirror.set_xlabel("m/z")
    title = f"Modified sequence: {mod_sequence}, charge: {charge}, retention time: {rt}"
    title_2 = f"Fragmentation: {fragm}, mass analyzer: {mass_analyzer}, collision energy aligned: {ce}"
    title_top = f"Top: experimental, raw file: {raw_file}, scan number: {scan_number}"
    title_bottom = f"Bottom: prediction, model: {model}, spectral angle: {sa}"
    ax_mirror.set_title(f"{title}\n{title_2}\n{title_top}\n{title_bottom}", fontsize=10)
    sup.mirror(top_spectrum, bot_spectrum, ax=ax_mirror)

    # KDE score plot
    sns.kdeplot(data=concat_target_decoy, x=score_col, hue="target", ax=ax_kde)

    ax_kde.axvline(score, color="black", linestyle="--", linewidth=1, label=f"score: {score:.2f}")
    ax_kde.legend(loc="upper right", fontsize=8)
    ax_kde.set_xlabel(score_col, fontsize=10)
    ax_kde.set_ylabel("Density", fontsize=10)

    # KDE abs rt diff plot
    sns.kdeplot(data=prosit_df, x="abs_rt_diff", ax=ax_rt)
    ax_rt.axvline(abs_rt_diff, color="black", linestyle="--", linewidth=1, label=f"abs rt diff: {abs_rt_diff:.2f}")
    ax_rt.legend(loc="upper right", fontsize=8)
    ax_rt.set_xlabel("abs rt diff", fontsize=10)
    ax_rt.set_ylabel("Density", fontsize=10)

    pdf.savefig(fig)


def plot_all(data_dir: Path, config: Config):
    """
    Generate all plots after a rescoring run.

    :param data_dir: the directory containing all inputs / outputs from either percolator or mokapot
    :param config: the configuration object
    """
    fdr_method = data_dir.stem
    prosit_df = pd.read_csv(data_dir / "rescore.tab", delimiter="\t")
    prosit_pep_target = pd.read_csv(data_dir / f"rescore.{fdr_method}.peptides.txt", delimiter="\t")
    prosit_pep_decoy = pd.read_csv(data_dir / f"rescore.{fdr_method}.decoy.peptides.txt", delimiter="\t")
    prosit_psms_target = pd.read_csv(data_dir / f"rescore.{fdr_method}.psms.txt", delimiter="\t")
    prosit_psms_decoy = pd.read_csv(data_dir / f"rescore.{fdr_method}.decoy.psms.txt", delimiter="\t")

    andromeda_pep_target = pd.read_csv(data_dir / f"original.{fdr_method}.peptides.txt", delimiter="\t")
    andromeda_pep_decoy = pd.read_csv(data_dir / f"original.{fdr_method}.decoy.peptides.txt", delimiter="\t")
    andromeda_psms_target = pd.read_csv(data_dir / f"original.{fdr_method}.psms.txt", delimiter="\t")
    andromeda_psms_decoy = pd.read_csv(data_dir / f"original.{fdr_method}.decoy.psms.txt", delimiter="\t")

    plot_score_distribution(
        prosit_pep_target,
        prosit_pep_decoy,
        "peptide",
        data_dir / "rescore_target_vs_decoys_peptide_bins.svg",
    )
    plot_score_distribution(
        prosit_psms_target,
        prosit_psms_decoy,
        "psm",
        data_dir / "rescore_target_vs_decoys_psm_bins.svg",
    )
    plot_score_distribution(
        andromeda_pep_target,
        andromeda_pep_decoy,
        "peptide",
        data_dir / "original_target_vs_decoys_peptide_bins.svg",
    )
    plot_score_distribution(
        andromeda_psms_target,
        andromeda_psms_decoy,
        "psm",
        data_dir / "original_target_vs_decoys_psm_bins.svg",
    )
    plot_sa_distribution(
        prosit_df,
        prosit_psms_target,
        prosit_psms_decoy,
        data_dir / "target_vs_decoys_sa_distribution.svg",
    )
    joint_plot(
        prosit_pep_target,
        prosit_pep_decoy,
        andromeda_pep_target,
        andromeda_pep_decoy,
        "peptide",
        data_dir / "rescore_original_joint_plot_peptide.svg",
    )
    joint_plot(
        prosit_psms_target,
        prosit_psms_decoy,
        andromeda_psms_target,
        andromeda_psms_decoy,
        "psm",
        data_dir / "rescore_original_joint_plot_psm.svg",
    )

    plot_gain_loss(prosit_pep_target, andromeda_pep_target, "peptide", data_dir / "peptide_1%_FDR.svg")
    plot_gain_loss(prosit_psms_target, andromeda_psms_target, "psm", data_dir / "psm_1%_FDR.svg")

    plot_pred_rt_vs_irt(prosit_df, prosit_psms_target, data_dir, "irt_vs_pred_rt.svg")

    base_mzml_path = os.path.abspath(os.path.join(data_dir, "../../spectra"))
    base_hdf5_path = os.path.abspath(os.path.join(data_dir, "../../data"))
    mirror_plots_dict = config.mirror_plots

    fdr_method = data_dir.stem

    prosit_psms_target = pd.read_csv(data_dir / f"rescore.{fdr_method}.psms.txt", delimiter="\t")
    prosit_psms_decoy = pd.read_csv(data_dir / f"rescore.{fdr_method}.decoy.psms.txt", delimiter="\t")

    pdf_path = data_dir / "mirror_plots.pdf"

    with PdfPages(pdf_path) as pdf:
        for raw_file, scan_numbers in mirror_plots_dict.items():
            mzml_path = os.path.join(base_mzml_path, f"{raw_file}.mzML")
            hdf5_path = os.path.join(base_hdf5_path, f"{raw_file}.mzml.pred.hdf5")

            mzml = ThermoRaw.read_mzml(source=mzml_path)
            spec_pred = Spectra.from_hdf5(hdf5_path)

            for scan_number in scan_numbers:
                plot_mirror_spectrum(
                    spec_pred,
                    mzml,
                    raw_file,
                    scan_number,
                    config,
                    prosit_df,
                    prosit_psms_target,
                    prosit_psms_decoy,
                    pdf,
                )


def plot_ce_ransac_model(
    sa_ce_df: pd.DataFrame, filename: Path, xlabel: str = "MASS", ylabel: str = "delta collision enegery", **kwargs
):
    """Generate plot (mass vs ce difference)."""
    df = sa_ce_df.reset_index()
    df = df[["MASS", "delta_collision_energy", "SPECTRAL_ANGLE"]]
    fig, ax = plt.subplots()
    sns.scatterplot(data=df, x="MASS", y="delta_collision_energy", hue="SPECTRAL_ANGLE", alpha=0.4, ax=ax)
    sns.regplot(
        data=df, x="MASS", y="delta_collision_energy", scatter=False, ci=None, line_kws={"linestyle": "--"}, ax=ax
    )
    ax.set(xlabel=xlabel, ylabel=ylabel, **kwargs)
    plt.savefig(filename, dpi=300)
    plt.close()
