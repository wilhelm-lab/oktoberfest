import logging
import subprocess

from picked_group_fdr.digestion_params import DigestionParams
from picked_group_fdr.pipeline import pipeline

from . import Config

logger = logging.getLogger(__name__)


def apply_quant(config: Config) -> None:
    """
    Call picked-group-FDR for means of quantification on search results with rescored Oktoberfest output.

    Integration of rescoring results only possible for MaxQuant at the moment. For msfragger and sage
    pciked_group_fdr is applied to the prior search results.

    :param config: config object containing all Oktoberfest parameters
    """
    if config.search_results.is_file():
        path = config.search_results.parent
    else:
        path = config.search_results

    fdr_dir = config.output / "results" / config.fdr_estimation_method

    # uses Percolator (or Mokapot) output as starting point
    okt_target_psms = str(fdr_dir / f"rescore.{config.fdr_estimation_method}.psms.txt")
    okt_decoy_psms = str(fdr_dir / f"rescore.{config.fdr_estimation_method}.decoy.psms.txt")

    # folder for picked group FDR output
    pg_fdr_dir = config.output / "picked_group_fdr/"
    pg_fdr_dir.mkdir(exist_ok=True, parents=True)

    # reading fasta file paths from config
    fasta_file = [str(config.inputs["library_input"])]

    # read in digestion parameters from config
    digest_params_list = [
        DigestionParams(
            enzyme=config.fasta_digest_options["enzyme"],
            digestion=config.fasta_digest_options["digestion"],
            min_length=config.fasta_digest_options["minLength"],
            max_length=config.fasta_digest_options["maxLength"],
            cleavages=config.fasta_digest_options["missedCleavages"],
            special_aas=config.fasta_digest_options["specialAas"],
            methionine_cleavage=True,  # default
            db=config.fasta_digest_options["db"],
            use_hash_key=True,  # default
        )
    ]

    if config.search_results_type == "maxquant":
        evi_in = str(path / "evidence.txt")

        # saves updated MQ evidence file under same name, but different directory
        evi_out = str(pg_fdr_dir / "evidence.txt")

        logger.info("Update the MaxQuant evidence files with the rescored PSMs from Oktoberfest")
        pipeline.run_update_evidence(
            [evi_in], [okt_target_psms, okt_decoy_psms], [evi_out], "prosit", suppress_missing_peptide_warning=True
        )

        logger.info("Generate MQ proteinGroup-like file from updated evidence file")
        protein_groups_out = str(pg_fdr_dir / "rescore.proteinGroups.txt")

        lfq_min_peptide_ratios = 1  # 1 per default
        pipeline.run_picked_group_fdr(
            [evi_out],
            protein_groups_out,
            fasta_file,
            digest_params_list,
            config.quantification,
            lfq_min_peptide_ratios,
            suppress_missing_peptide_warning=True,
        )

        logger.info("FDR-filtering the MQ proteinGroup-like file")
        fdr_cut = 0.01
        pipeline.run_filter_fdr_maxquant(
            [protein_groups_out],
            str(pg_fdr_dir / f"rescore.proteinGroups.fdr{int(fdr_cut * 100)}.txt"),
            fdr_cutoff=fdr_cut,
        )

    if config.search_results_type == "sage":
        sage_results = str(path / "results.sage.tsv")
        sage_lfq = str(path / "lfq.tsv")

        cmd = f"python3 -u -m picked_group_fdr \
                --fasta {fasta_file} \
                --sage_results {sage_results} \
                --sage_lfq_tsv {sage_lfq} \
                --protein_groups_out {pg_fdr_dir}/combined_protein.tsv \
                --output_format fragpipe \
                --do_quant \
                --lfq_min_peptide_ratios 1 \
                --methods sage"

        logger.info(f"Starting picked-group-fdr with command {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    if config.search_results_type == "msfragger":
        fragpipe_psms = str(path / "*" / "psm.tsv")
        # reusing IonQuant results by default
        combined_ion = str(path / "combined_ion.tsv")

        cmd = f"python3 -u -m picked_group_fdr \
                --fasta {fasta_file} \
                --fragpipe_psm {fragpipe_psms} \
                --combined_ion {combined_ion} \
                --protein_groups_out {pg_fdr_dir}/combined_protein.tsv \
                --do_quant \
                --lfq_min_peptide_ratios 1 \
                --methods fragpipe"

        logger.info(f"Starting picked-group-fdr with command {cmd}")
        subprocess.run(cmd, shell=True, check=True)
