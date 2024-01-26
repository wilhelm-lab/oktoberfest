"""quantification module, based on picked-group-fdr

"""

import logging
from pathlib import Path

from sklearn.linear_model import LinearRegression, RANSACRegressor

from oktoberfest import __copyright__, __version__
from oktoberfest import rescore as re

from . import Config, JobPool, ProcessStep

import picked_group_fdr.pipeline as picked_group_fdr
from picked_group_fdr.digestion_params import DigestionParams


logger = logging.getLogger(__name__)


def apply_quant(config: Config):
    fdr_dir = config.output / "results" / config.fdr_estimation_method

    # currently only works for: MQ evidence file - TODO integrate MQ msms.txt ? 
    if config.search_results_type == "maxquant":
        # config.search is already a Path, but picked-group-fdr can't handle it
        evi_in = str(config.search_results.parent / "evidence.txt")
            
        # uses Percolator (or Mokapot) output as starting point
        okt_target_psms = str(fdr_dir / f"rescore.{config.fdr_estimation_method}.psms.txt")
        okt_decoy_psms = str(fdr_dir / f"rescore.{config.fdr_estimation_method}.decoy.psms.txt")
        # folder for picked group FDR output
        pg_fdr_dir = config.output / "picked_group_fdr/"
        pg_fdr_dir.mkdir(exist_ok=True, parents=True)
        # saves updated MQ evidence file under same name, but diff. dir.
        evi_rescor_integr = str(pg_fdr_dir / "evidence.txt")

        logger.info("Update the MaxQuant evidence files with the rescored PSMs from Oktoberfest") 
        picked_group_fdr.run_update_evidence(
            [evi_in],
            [okt_target_psms, okt_decoy_psms],
            [evi_rescor_integr],
            "prosit",
        )

        logger.info("Generate MQ proteinGroup-like file from updated evidence file")    
        protein_groups_out = str(pg_fdr_dir / "rescore.proteinGroups.txt")
        digest_params_list = [
            DigestionParams(
                config.fasta_digest_options["enzyme"], 
                config.fasta_digest_options["digestion"],
                config.fasta_digest_options["minLength"],
                config.fasta_digest_options["maxLength"],
                config.fasta_digest_options["missedCleavages"],
                config.fasta_digest_options["specialAas"],
                #True, # methionine_cleavage hardcode this way in picked group fdr
                config.fasta_digest_options["db"],   
                #config.fasta_digest_options["digestion"] == True, # use_hash_key defined that way in picked group fdr
            ),
        ]
        lfq_min_peptide_ratios = 1   # TODO 1 per default? / param?
        picked_group_fdr.run_picked_group_fdr(
            [evi_rescor_integr],   # Oktoberfest works on a single one of these
            protein_groups_out,
            [str(config.inputs["library_input"])],
            digest_params_list,
            config.quantification,
            lfq_min_peptide_ratios,
        )

        logger.info("FDR-filtering the MQ proteinGroup-like file")    
        fdr_cut = 0.01   # TODO default or param?
        picked_group_fdr.run_filter_fdr_maxquant(
            [protein_groups_out], 
            str(pg_fdr_dir / f"rescore.proteinGroups.fdr{int(fdr_cut * 100)}.txt"), 
            fdr_cutoff=fdr_cut
        )

    # TODO integrate MSFragger 
    # if config.search_results_type == "msfragger":  
    #     # Option where existing IonQuant results are reused

    #     # TODO adjust below
        
    #     # config.search is already a Path, but picked-group-fdr can't handle it
    #     evi_in = str(config.search_results.parent / "evidence.txt")
            
    #     # uses Percolator (or Mokapot) output as starting point
    #     okt_target_psms = str(fdr_dir / f"rescore.{config.fdr_estimation_method}.psms.txt")
    #     okt_decoy_psms = str(fdr_dir / f"rescore.{config.fdr_estimation_method}.decoy.psms.txt")
    #     # folder for picked group FDR output
    #     pg_fdr_dir = config.output / "picked_group_fdr/"
    #     pg_fdr_dir.mkdir(exist_ok=True, parents=True)
    #     # saves updated MQ evidence file under same name, but diff. dir.
    #     evi_rescor_integr = str(pg_fdr_dir / "evidence.txt")

    #     logger.info("Update the MaxQuant evidence files with the rescored PSMs from Oktoberfest") 

    # Your best bet are probably the functions in https://github.com/kusterlab/picked_group_fdr/tree/develop/picked_group_fdr/quant
    # the closest thing [to the MQ approach] are probably these two scripts:
    # https://github.com/kusterlab/picked_group_fdr/blob/develop/picked_group_fdr/pipeline/sage_quantification.py
    # https://github.com/kusterlab/picked_group_fdr/blob/develop/picked_group_fdr/pipeline/update_fragpipe_results.py