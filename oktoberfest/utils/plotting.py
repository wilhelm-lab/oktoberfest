from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Set the default fontsize and linewidth
plt.rcParams.update({"font.size": 14, "axes.linewidth": 1.5, "xtick.major.width": 1.5, "ytick.major.width": 1.5})


def plot_target_decoy(target: pd.DataFrame, decoy: pd.DataFrame, type: str, search_type: str, directory: Path):
    """Generate target-decoy distribution of the score."""
    plt.figure(figsize=(8, 6))
    bins = np.linspace(-3, 2, 15)
    plt.hist(target.score, bins, label="Targets", rwidth=0.5, color="#48AF00")
    plt.hist(decoy.score, bins, label="Decoys", rwidth=0.5, color="#FE7312")
    plt.xlabel("Score", size=14)
    plt.title(f"{search_type} Target vs Decoys ({type})")
    plt.legend(loc="upper right")
    plt.savefig(directory / f"{search_type}_Target_vs_Decoys_{type}_bins.png", dpi=300)
    plt.plot()
    plt.close()


def joint_plot(
    prosit_target: pd.DataFrame,
    prosit_decoy: pd.DataFrame,
    andromeda_target: pd.DataFrame,
    andromeda_decoy: pd.DataFrame,
    type: str,
    directory: Path,
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
    df_all.dropna(inplace=True)
    jplot = sns.jointplot(
        data=df_all,
        x="and_score",
        y="prosit_score",
        marginal_ticks=True,
        hue=df_all.type,
        palette=["#15853B", "#E2700E"],
        ratio=2,
        height=10,
        joint_kws={"rasterized": True, "edgecolor": "none", "s": 10},
    )
    jplot.ax_joint.set_ylabel("rescored_score")
    jplot.ax_joint.set_xlabel("original_score")
    plt.savefig(directory / f"Rescored_Original_joint_plot_{type}.svg", dpi=300)
    plt.plot()
    plt.close()


def plot_gain_loss(prosit_target: pd.DataFrame, andromeda_target: pd.DataFrame, type: str, directory: Path):
    """Generate gain-loss plot (peptides/PSMs 1% FDR)."""
    if type == "Peptides":
        join_col = "peptide"
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
    ax.set_ylabel("Percentage")
    ax.set_axisbelow(True)
    ax.yaxis.grid(color="black")
    ax.tick_params(axis="y", which="major")
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)

    legend_label = ["Common", "Gained", "Lost"]
    plt.legend(legend_label, ncol=1, bbox_to_anchor=([1.2, 0.5, 0, 0]), frameon=False)
    plt.title(f"{type} 1% FDR\n")
    plt.savefig(directory / f"{type}_1%_FDR.svg", dpi=300, bbox_inches="tight")
    plt.plot()
    plt.close()


def plot_mean_sa_ce(sa_ce_df: pd.DataFrame, filename: Union[str, Path], best_ce: float):
    """Generate plot (ce vs mean sa)."""
    df = sa_ce_df.to_frame()
    df = df.reset_index()
    df = df[["COLLISION_ENERGY", "SPECTRAL_ANGLE"]]
    fig, ax = plt.subplots(figsize=(8, 8))
    sns.scatterplot(data=df, x="COLLISION_ENERGY", y="SPECTRAL_ANGLE", ax=ax)
    ax.axvline(x=best_ce, color="red")
    plt.grid()
    plt.savefig(filename, dpi=300)
    plt.plot()
    plt.close()


def plot_violin_sa_ce(df: pd.DataFrame, filename: Union[str, Path], best_ce: float):
    """Generate plot (ce vs mean sa)."""
    fig, ax = plt.subplots(figsize=(8, 8))
    sns.violinplot(data=df, x="COLLISION_ENERGY", y="SPECTRAL_ANGLE", ax=ax, color="#1f77b4")
    ax.axvline(x=best_ce, color="red")
    plt.grid()
    plt.savefig(filename, dpi=300)
    plt.plot()
    plt.close()


def plot_pred_rt_vs_irt(prosit_df: pd.DataFrame, prosit_target: pd.DataFrame, directory: Path):
    """Generate pred rt vs irt plot."""
    targets = prosit_df.merge(prosit_target, how="inner", left_on="SpecId", right_on="PSMId")
    targets = targets.loc[targets["q-value"] < 0.01, ["SpecId", "RT", "iRT", "pred_RT"]]
    targets.columns = ["SpecId", "experimental RT", "aligned RT", "predicted iRT"]
    targets["rawfile"] = targets["SpecId"].str.split("-", n=1).str[0] + " aligned RT"
    targets.sort_values("predicted iRT")
    for rawfile in targets["rawfile"].unique():
        fig, ax = plt.subplots(figsize=(8, 8))
        sns.scatterplot(
            data=targets[targets["rawfile"] == rawfile],
            x="predicted iRT",
            y="experimental RT",
            label="predicted iRT",
            ax=ax,
        )
        sns.lineplot(
            data=targets[targets["rawfile"] == rawfile], x="predicted iRT", y="aligned RT", label="alignment", ax=ax
        )
        # plot.plot(targets["RT"], targets["pred_RT"], ".", c="b", label="original")
        # plt.plot(targets["iRT"], targets["pred_RT"], ".", c="r", label="smoothed")
        ax.set_ylabel("(aligned) experimental RT")
        ax.set_xlabel("predicted iRT")
        plt.title("RT alignment")
        plt.legend(labels=("predicted iRT", "alignment"), loc="best", fancybox=True, shadow=True)
        plt.grid()
        plt.savefig(directory / f"{rawfile}_pred_rt_vs_irt.svg", dpi=300)
        plt.plot()
        plt.close()


def plot_all(percolator_path: Path):
    """Generate all plots and save them as png in the percolator folder."""
    prosit_pep_target = pd.read_csv(percolator_path / "rescore.target.peptides", delimiter="\t")
    prosit_pep_decoy = pd.read_csv(percolator_path / "rescore.decoy.peptides", delimiter="\t")
    prosit_psms_target = pd.read_csv(percolator_path / "rescore.target.psms", delimiter="\t")
    prosit_psms_decoy = pd.read_csv(percolator_path / "rescore.decoy.psms", delimiter="\t")

    andromeda_pep_target = pd.read_csv(percolator_path / "original.target.peptides", delimiter="\t")
    andromeda_pep_decoy = pd.read_csv(percolator_path / "original.decoy.peptides", delimiter="\t")
    andromeda_psms_target = pd.read_csv(percolator_path / "original.target.psms", delimiter="\t")
    andromeda_psms_decoy = pd.read_csv(percolator_path / "original.decoy.psms", delimiter="\t")

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

    prosit_df = pd.read_csv(percolator_path / "rescore.tab", delimiter="\t")
    plot_pred_rt_vs_irt(prosit_df, prosit_psms_target, percolator_path)
