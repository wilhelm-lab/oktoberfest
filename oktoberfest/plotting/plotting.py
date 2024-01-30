from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

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
    """
    score_col, _, _, _ = _check_columns(target)

    plt.figure(figsize=(8, 6))
    bins = np.linspace(-3, 2, 15)
    plt.hist(target[score_col], bins, label="Targets", rwidth=0.5, color="#48AF00", alpha=1.0)
    plt.hist(decoy[score_col], bins, label="Decoys", rwidth=0.5, color="#FE7312", alpha=0.7)
    plt.xlabel("Score")
    plt.legend(loc="upper right")
    plt.title(f"Score Distribution ({level.capitalize()})")
    plt.savefig(filename, dpi=300)
    plt.plot()
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
    jplot.ax_joint.set_ylabel("Score\n(peptide property prediction)")
    jplot.ax_joint.set_xlabel("Score\n(search engine)")
    jplot.fig.suptitle(f"Score distribution ({level.capitalize()})", y=0.99)
    plt.savefig(filename, dpi=300)
    plt.plot()
    plt.close()


def plot_gain_loss(prosit_target: pd.DataFrame, andromeda_target: pd.DataFrame, level: str, filename: Union[str, Path]):
    """
    Generate venn barplots to show lost, common and shared targets below 1% FDR attributed to peptide property predictions.

    :param prosit_target: mokapot / percolator target output for rescoring with peptide property prediction
    :param andromeda_target: mokapot / percolator target output for rescoring without peptide property prediction
    :param level: The level on which to produce the comparison. Can be either "peptide" or "psm"
    :param filename: the path to the location used for storing the plot

    :raises ValueError: if a wrong level is provided
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
    ax.set_ylabel("number of lost-common-shared targets below 1% FDR")
    ax.set_axisbelow(True)
    ax.yaxis.grid(color="black")
    ax.tick_params(axis="y", which="major")
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)

    legend_label = ["Common", "Gained", "Lost"]
    plt.legend(legend_label, ncol=1, bbox_to_anchor=([1.2, 0.5, 0, 0]), frameon=False)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.plot()
    plt.close()


def plot_mean_sa_ce(sa_ce_df: pd.DataFrame, filename: Union[str, Path]):
    """
    Generate dotplot for spectral angle distribution over range of collision energies used for fragment intensity prediction.

    :param sa_ce_df: a dataframe containing the two columns "COLLISION_ENERGY", "SPECTRAL_ANGLE".
    :param filename: the path to the location used for storing the plot
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    sns.scatterplot(data=sa_ce_df, x="COLLISION_ENERGY", y="SPECTRAL_ANGLE", ax=ax)
    ax.axvline(x=sa_ce_df["COLLISION_ENERGY"][sa_ce_df["SPECTRAL_ANGLE"].idxmax()], color="red")
    plt.grid()
    plt.savefig(filename, dpi=300)
    plt.plot()
    plt.close()


def plot_violin_sa_ce(sa_ce_df: pd.DataFrame, filename: Union[str, Path]):
    """
    Generate violinplot for spectral angle distribution over range of collision energies used for fragment intensity prediction.

    :param sa_ce_df: a dataframe containing the two columns "COLLISION_ENERGY", "SPECTRAL_ANGLE".
    :param filename: the path to the location used for storing the plot
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    sns.violinplot(data=sa_ce_df, x="COLLISION_ENERGY", y="SPECTRAL_ANGLE", ax=ax, color="#1f77b4")
    ax.axvline(x=sa_ce_df["COLLISION_ENERGY"][sa_ce_df["SPECTRAL_ANGLE"].idxmax()], color="red")
    plt.grid()
    plt.savefig(filename, dpi=300)
    plt.plot()
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
        plt.plot()
        plt.close()


def plot_sa_distribution(prosit_df: pd.DataFrame, target_df: pd.DataFrame, decoy_df: pd.DataFrame, filename: Path):
    """Generate spectral angle distribution for targets and decoys.

    :param prosit_df: mokapot / percolator input tab for rescoring with peptide property prediction
    :param target_df: mokapot / percolator target output for rescoring with peptide property prediction on the psm level
    :param decoy_df: mokapot / percolator decoy output for rescoring with peptide property prediction on the psm level
    :param filename: the path to the location used for storing the plot
    """
    _, _, _, psm_col = _check_columns(target_df)
    target = prosit_df.merge(target_df, how="inner", left_on="SpecId", right_on=psm_col)
    decoy = prosit_df.merge(decoy_df, how="inner", left_on="SpecId", right_on=psm_col)
    plt.figure(figsize=(8, 6))
    bins = np.linspace(0, 1, 15)
    plt.hist(target.spectral_angle, bins, label="Targets", rwidth=0.5, color="#48AF00", alpha=1.0)
    plt.hist(decoy.spectral_angle, bins, label="Decoys", rwidth=0.5, color="#FE7312", alpha=0.7)
    plt.xlabel("Spectral angle", size=14)
    plt.title("Target vs Decoys Spectral Angle Distribution")
    plt.legend(loc="upper right")
    plt.savefig(filename, dpi=300)
    plt.plot()
    plt.close()


def plot_all(data_dir: Path):
    """
    Generate all plots after a rescoring run.

    :param data_dir: the directory containing all inputs / outputs from either percolator or mokapot.

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


def plot_ce_ransac_model(
    sa_ce_df: pd.DataFrame,
    filename: Path,
    xlabel: str = "MASS",
    ylabel: str = "delta collision enegery",
    *args,
    **kwargs,
):
    """Generate plot (mass vs ce difference)."""
    df = sa_ce_df.reset_index()
    df = df[["MASS", "delta_collision_energy", "SPECTRAL_ANGLE"]]
    fig, ax = plt.subplots()
    sns.scatterplot(data=df, x="MASS", y="delta_collision_energy", hue="SPECTRAL_ANGLE", alpha=0.4, ax=ax)
    sns.regplot(
        data=df, x="MASS", y="delta_collision_energy", scatter=False, ci=None, line_kws={"linestyle": "--"}, ax=ax
    )
    ax.set(xlabel=xlabel, ylabel=ylabel, *args, **kwargs)
    plt.savefig(filename, dpi=300)
