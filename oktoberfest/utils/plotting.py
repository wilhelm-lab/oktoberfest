import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import Config


def plot_target_decoy(target: pd.DataFrame, decoy: pd.DataFrame, type: str, search_type: str, directory: str):
    """Generate target-decoy distribution of the score."""
    plt.figure(figsize=(8, 6))
    bins = np.linspace(-3, 2, 15)
    plt.hist(target.score, bins, label="Targets", rwidth=0.5, color="#48AF00")
    plt.hist(decoy.score, bins, label="Decoys", rwidth=0.5, color="#FE7312")
    plt.xlabel("Score", size=14)
    plt.title(f"{search_type} Target vs Decoys ({type})")
    plt.legend(loc="upper right")
    plt.savefig(directory + f"/{search_type}_Target_vs_Decoys_{type}_bins.png", dpi=300)


def joint_plot(
    prosit_target: pd.DataFrame,
    prosit_decoy: pd.DataFrame,
    andromeda_target: pd.DataFrame,
    andromeda_decoy: pd.DataFrame,
    type: str,
    directory: str,
):
    """Generate joint plot (correlation between Prosit and Andromeda score)."""
    if type == "Peptides":
        join_col = "proteinIds"
    else:
        join_col = "PSMId"

    targets = prosit_target.merge(andromeda_target, on=join_col, how="outer", suffixes=["", "_"], indicator=True)
    decoys = prosit_decoy.merge(andromeda_decoy, on=join_col, how="outer", suffixes=["", "_"], indicator=True)
    df_targets = pd.DataFrame()
    df_targets["prosit_score"] = targets.score
    df_targets["and_score"] = targets.score_
    df_targets["type"] = "targets"
    df_targets["color"] = "#48AF00"

    df_decoys = pd.DataFrame()
    df_decoys["prosit_score"] = decoys.score
    df_decoys["and_score"] = decoys.score_
    df_decoys["type"] = "decoys"
    df_decoys["color"] = "#FE7312"

    df_all = pd.concat([df_targets, df_decoys], axis=0)
    df_all = df_all[~df_all.index.duplicated()]
    jplot = sns.jointplot(
        data=df_all,
        x="and_score",
        y="prosit_score",
        marginal_ticks=True,
        hue=df_all.type,
        palette=["#48AF00", "#FE7312"],
        ratio=2,
        height=10,
    )
    jplot.ax_joint.set_ylabel("rescored_score")
    jplot.ax_joint.set_xlabel("original_score")
    plt.savefig(directory + f"/Rescored_Original_joint_plot_{type}.png", dpi=300)


def plot_gain_loss(
    prosit_target: pd.DataFrame,
    andromeda_target: pd.DataFrame,
    type: str,
    directory: str,
):
    """Generate gain-loss plot (peptides/PSMs 1% FDR)."""
    if type == "Peptides":
        join_col = "proteinIds"
    else:
        join_col = "PSMId"

    andromeda_target = andromeda_target[andromeda_target["q-value"] < 0.01]
    prosit_target = prosit_target[prosit_target["q-value"] < 0.01]
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
            fontsize=12,
            fontweight="bold",
        )
        plt.text(
            r2.get_x() + r2.get_width() / 2.0,
            h1 + h2 / 2,
            "%d" % v2,
            ha="center",
            va="bottom",
            color="black",
            fontsize=12,
        )
        plt.text(
            r3.get_x() + r3.get_width() / 2.0,
            -0.025 * gained,
            "%d" % -v3,
            ha="center",
            va="bottom",
            color="black",
            fontsize=12,
        )

    plt.ylim(-lost - 100, h1 + h2 + 30)
    # remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    # grid
    ax.set_ylabel("Percentage", fontsize=14)
    ax.set_axisbelow(True)
    ax.yaxis.grid(color="gray")
    ax.tick_params(axis="y", which="major", labelsize=13)
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)

    legend_label = ["Common", "Gained", "Lost"]
    plt.legend(legend_label, ncol=1, bbox_to_anchor=([1.2, 0.5, 0, 0]), frameon=False)
    plt.title(f"{type} 1% FDR\n", fontsize=14)
    plt.savefig(directory + f"/{type}_1%_FDR.png", dpi=300, bbox_inches="tight")


def plot_mean_sa_ce(sa_ce_df: pd.DataFrame, directory: str, raw_file_name: str):
    """Generate plot (ce vs mean sa)."""
    directory = directory + ""
    directory = directory.replace("/mzML", "")
    directory = directory.replace("/percolator", "")
    df = sa_ce_df.to_frame()
    df = df.reset_index()
    df = df[["COLLISION_ENERGY", "SPECTRAL_ANGLE"]]
    sns.lmplot(data=df, x="COLLISION_ENERGY", y="SPECTRAL_ANGLE", ci=None, order=2, truncate=False)
    plt.savefig(directory + "/" + raw_file_name + "mean_spectral_angle_ce.png", dpi=300)


def plot_all(percolator_path: str, config: Config):
    """Generate all plots and save them as png in the percolator folder."""
    fdr_estimation_method = config.fdr_estimation_method
    if fdr_estimation_method == "mokapot":
        files = [
            "/rescore.mokapot.peptides.txt",
            "/rescore.mokapot.decoy.peptides.txt",
            "/rescore.mokapot.psms.txt",
            "/rescore.mokapot.decoy.psms.txt",
            "/original.mokapot.peptides.txt",
            "/original.mokapot.decoy.peptides.txt",
            "/original.mokapot.psms.txt",
            "/original.mokapot.decoy.psms.txt",
        ]

        for f in files:
            prefix = "rescore" if f.startswith("/rescore") else "original"
            target_decoy = "decoy." if "decoy" in f else "target."

            file = percolator_path + f
            file_renamed = percolator_path + "/" + prefix + "_" + target_decoy + f.split(".")[-2]
            df = pd.read_csv(file, delimiter="\t")
            df.rename(
                columns=(
                    {
                        "mokapot score": "score",
                        "mokapot q-value": "q-value",
                        "Proteins": "proteinIds",
                        "SpecId": "PSMId",
                    }
                ),
                inplace=True,
            )
            df.to_csv(file, sep="\t", index=False)
            os.rename(file, file_renamed)

    prosit_pep_target = pd.read_csv(percolator_path + "/rescore_target.peptides", delimiter="\t")
    prosit_pep_decoy = pd.read_csv(percolator_path + "/rescore_decoy.peptides", delimiter="\t")
    prosit_psms_target = pd.read_csv(percolator_path + "/rescore_target.psms", delimiter="\t")
    prosit_psms_decoy = pd.read_csv(percolator_path + "/rescore_decoy.psms", delimiter="\t")

    andromeda_pep_target = pd.read_csv(percolator_path + "/original_target.peptides", delimiter="\t")
    andromeda_pep_decoy = pd.read_csv(percolator_path + "/original_decoy.peptides", delimiter="\t")
    andromeda_psms_target = pd.read_csv(percolator_path + "/original_target.psms", delimiter="\t")
    andromeda_psms_decoy = pd.read_csv(percolator_path + "/original_decoy.psms", delimiter="\t")

    plot_target_decoy(prosit_pep_target, prosit_pep_decoy, "Peptides", "Rescore", percolator_path)
    plot_target_decoy(prosit_psms_target, prosit_psms_decoy, "PSMs", "Rescore", percolator_path)
    plot_target_decoy(andromeda_pep_target, andromeda_pep_decoy, "Peptides", "Original", percolator_path)
    plot_target_decoy(andromeda_psms_target, andromeda_psms_decoy, "PSMs", "Original", percolator_path)

    joint_plot(
        prosit_pep_target, prosit_pep_decoy, andromeda_pep_target, andromeda_pep_decoy, "Peptides", percolator_path
    )
    joint_plot(
        prosit_psms_target, prosit_psms_decoy, andromeda_psms_target, andromeda_psms_decoy, "PSMs", percolator_path
    )
    plot_gain_loss(prosit_pep_target, andromeda_pep_target, "Peptides", percolator_path)
    plot_gain_loss(prosit_psms_target, andromeda_psms_target, "PSMs", percolator_path)
