Features for target / decoy separation
======================================

Oktoberfest creates tab files for target / decoy separation and fdr estimation using Percolator and mokapot. This provides an overview of the features that are used for rescoring with and without peptide property prediction.


Always present
--------------

.. table::
    :class: fixed-table

    +----------------+--------------------------------------------------------------------------------------+
    | Feature        | Description                                                                          |
    +================+======================================================================================+
    | SpecId         | internal spectrum ID                                                                 |
    +----------------+--------------------------------------------------------------------------------------+
    | ScanNr         | scan number hash                                                                     |
    +----------------+--------------------------------------------------------------------------------------+
    | Peptide        | peptide sequence                                                                     |
    +----------------+--------------------------------------------------------------------------------------+
    | missedCleavages| number of missed cleavages                                                           |
    +----------------+--------------------------------------------------------------------------------------+
    | sequence_length| peptide sequence length                                                              |
    +----------------+--------------------------------------------------------------------------------------+
    | Mass           | experimental mass                                                                    |
    +----------------+--------------------------------------------------------------------------------------+
    | ExpMass        | fixed constant to allow PSMs with different theoretical mass to compete              |
    |                | for the same PSM in target-decoy competition                                         |
    +----------------+--------------------------------------------------------------------------------------+
    | Label          | 1 indicates target, -1 decoy                                                         |
    +----------------+--------------------------------------------------------------------------------------+
    | Protein        | protein description                                                                  |
    +----------------+--------------------------------------------------------------------------------------+
    | deltaM_ppm     | delta mass in ppm                                                                    |
    +----------------+--------------------------------------------------------------------------------------+
    | absDeltaM_ppm  | absolute delta mass in ppm                                                           |
    +----------------+--------------------------------------------------------------------------------------+
    | deltaM_da      | delta mass in Da                                                                     |
    +----------------+--------------------------------------------------------------------------------------+
    | absDeltaM_da   | absolute delta mass in Da                                                            |
    +----------------+--------------------------------------------------------------------------------------+
    | Charge2        | boolean for charge state 2                                                           |
    +----------------+--------------------------------------------------------------------------------------+
    | Charge3        | boolean for charge state 3                                                           |
    +----------------+--------------------------------------------------------------------------------------+
    | KR             | number of K or R amino acids in sequence                                             |
    +----------------+--------------------------------------------------------------------------------------+

Features related to Prosit
--------------------------

.. table::
    :class: fixed-table

    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | Feature                                          | Description                                                                                                               |
    +==================================================+===========================================================================================================================+
    | spectral_angle                                   | Normalized spectral contrast angle (SA) on all potential y- and b-ions                                                    |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_seen_nonzero                                | number of observed non-zero ions with non-zero predicted intensity                                                        |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | not_pred_seen                                    | number of observed non-zero ions not predicted                                                                            |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | not_pred_seen_b                                  | number of observed non-zero b-ions not predicted                                                                          |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | not_pred_seen_y                                  | number of observed non-zero y-ions not predicted                                                                          |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_nonZero_fragments                           | number of ions with non-zero predicted intensity                                                                          |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_nonZero_b                                   | number of b-ions with non-zero predicted intensity                                                                        |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_nonZero_y                                   | number of y-ions with non-zero predicted intensity                                                                        |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_not_seen                                    | number of ions with non-zero predicted intensity not observed                                                             |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_not_seen_b                                  | number of b-ions with non-zero predicted intensity not observed                                                           |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_not_seen_y                                  | number of y-ions with non-zero predicted intensity not observed                                                           |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_seen_nonzero_b                              | number of observed b-ions with non-zero predicted intensity                                                               |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_seen_nonzero_y                              | number of observed y-ions with non-zero predicted intensity                                                               |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_seen_zero                                   | number of predicted and observed zero-intensity ions                                                                      |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_seen_zero_b                                 | number of predicted and observed zero-intensity b-ions                                                                    |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | pred_seen_zero_y                                 | number of predicted and observed zero-intensity y-ions                                                                    |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | raw_nonZero_fragments                            | number of observed ions                                                                                                   |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | raw_nonZero_b                                    | number of observed b-ions                                                                                                 |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | raw_nonZero_y                                    | number of observed y-ions                                                                                                 |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_not_pred_seen                                | (number of observed ions not predicted) / (number of theoretically observable fragments [length*charge^2])                |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_not_pred_seen_b                              | (number of observed b-ions not predicted) / (number of theoretically observable fragments [length*charge])                |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_not_pred_seen_y                              | (number of observed y-ions not predicted) / (number of theoretically observable fragments [length*charge])                |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_nonZero_b                               | (number of observed b-ions not predicted) / (number of theoretically observable fragments [length*charge])                |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_nonZero_y                               | (number of observed y-ions not predicted) / (number of theoretically observable fragments [length*charge^2])              |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_not_seen                                | (number of ions predicted, but not observed) / (number of theoretically observable fragments [length*charge^2])           |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_not_seen_b                              | (number of b-ions ions predicted, but not observed) / (number of theoretically observable fragments [length*charge])      |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_not_seen_y                              | (number of y-ions ions predicted, but not observed) / (number of theoretically observable fragments [length*charge])      |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_seen_nonzero                            | (number of observed ions predicted ) / (number of theoretically observable fragments [length*charge^2])                   |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_seen_nonzero_b                          | (number of observed b-ions predicted) / (number of theoretically observable fragments [length*charge])                    |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_seen_nonzero_y                          | (number of observed y-ions predicted) / (number of theoretically observable fragments [length*charge])                    |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_seen_zero                               | (number of predicted and observed zero-intensity ions) / (number of theoretically observable fragments [length*charge^2]) |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_seen_zero_b                             | (number of predicted and observed zero-intensity b-ions) / (number of theoretically observable fragments [length*charge]) |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_pred_seen_zero_y                             | (number of predicted and observed zero-intensity y-ions) / (number of theoretically observable fragments [length*charge]) |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_raw_nonZero_fragments                        | (number of observed ions) / (number of theoretically observable fragments [length*charge^2])                              |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_raw_nonZero_b                                | (number of observed b-ions) / (number of theoretically observable fragments [length*charge])                              |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | rel_raw_nonZero_y                                | (number of observed y-ions) / (number of theoretically observable fragments [length*charge])                              |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_not_pred_seen2pred_nonZero_fragments     | (number of observed non-zero ions not predicted) / (number of ions with non-zero predicted intensity)                     |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_not_pred_seen_b2pred_nonZero_b           | (number of observed non-zero b-ions not predicted) / (number of b-ions with non-zero predicted intensity)                 |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_not_pred_seen_y2pred_nonZero_y           | (number of observed non-zero y-ions not predicted) / (number of y-ions with non-zero predicted intensity)                 |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_nonZero_b2pred_nonZero_b            | (number of b-ions with non-zero predicted intensity) / (number of b-ions with non-zero predicted intensity)               |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_nonZero_y2pred_nonZero_y            | (number of y-ions with non-zero predicted intensity) / (number of y-ions with non-zero predicted intensity)               |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_not_seen_b2pred_nonZero_b           | (number of b-ions with non-zero predicted intensity not observed) / (number of b-ions with non-zero predicted intensity)  |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_not_seen_y2pred_nonZero_y           | (number of y-ions with non-zero predicted intensity not observed) / (number of y-ions with non-zero predicted intensity)  |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_not_seen2pred_nonZero_fragments     | (number of ions with non-zero predicted intensity not observed) / (number of ions with non-zero predicted intensity)      |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_seen_nonzero_b2pred_nonZero_b       | (number of observed b-ions with non-zero predicted intensity) / (number of b-ions with non-zero predicted intensity)      |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_seen_nonzero_y2pred_nonZero_y       | (number of observed y-ions with non-zero predicted intensity) / (number of y-ions with non-zero predicted intensity)      |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_seen_nonzero2pred_nonZero_fragments | (number of observed non-zero ions with non-zero predicted intensity) / (number of ions with non-zero predicted intensity) |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_seen_zero_b2pred_nonZero_b          | (number of predicted and observed zero-intensity b-ions) / (number of b-ions with non-zero predicted intensity)           |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_seen_zero_y2pred_nonZero_y          | (number of predicted and observed zero-intensity y-ions) / (number of y-ions with non-zero predicted intensity)           |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | relpred_pred_seen_zero2pred_nonZero_fragments    | (number of predicted and observed zero-intensity ions) / (number of ions with non-zero predicted intensity)               |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+
    | abs_rt_diff                                      | absolute difference between retention time and aligned predicted retention time                                           |
    +--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------+

Features adopted from MS2PIP/MS2Rescore
---------------------------------------

.. table::
    :class: fixed-table

    +-------------------------------------+--------------------------------------------------+
    | Feature                             | Description                                      |
    +=====================================+==================================================+
    | pearson_corr                        | Pearson correlation on all potential y- and      |
    |                                     | b-ions                                           |
    +-------------------------------------+--------------------------------------------------+
    | spearman_corr                       | Spearman correlation on all potential y- and     |
    |                                     | b-ions                                           |
    +-------------------------------------+--------------------------------------------------+
    | mse                                 | Mean square error                                |
    +-------------------------------------+--------------------------------------------------+
    | cos                                 | Cosine similarity                                |
    +-------------------------------------+--------------------------------------------------+
    | std_abs_diff                        | Standard deviation of the absolute differences   |
    +-------------------------------------+--------------------------------------------------+
    | abs_diff_Q3                         | Quantile 3 of the absolute differences           |
    +-------------------------------------+--------------------------------------------------+
    | abs_diff_Q2                         | Quantile 2 of the absolute differences           |
    +-------------------------------------+--------------------------------------------------+
    | abs_diff_Q1                         | Quantile 1 of the absolute differences           |
    +-------------------------------------+--------------------------------------------------+
    | min_abs_diff                        | Minimum absolute difference                      |
    +-------------------------------------+--------------------------------------------------+
    | max_abs_diff                        | Maximum absolute difference                      |
    +-------------------------------------+--------------------------------------------------+
    | spectral_angle_single_charge        | Normalized spectral contrast angle (SA) on       |
    |                                     | singly charged ions                              |
    +-------------------------------------+--------------------------------------------------+
    | spectral_angle_double_charge        | Normalized spectral contrast angle (SA) on       |
    |                                     | doubly charged ions                              |
    +-------------------------------------+--------------------------------------------------+
    | spectral_angle_triple_charge        | Normalized spectral contrast angle (SA) on       |
    |                                     | triply charged ions                              |
    +-------------------------------------+--------------------------------------------------+
    | spectral_angle_b_ions               | Normalized spectral contrast angle (SA) on       |
    |                                     | b-ions                                           |
    +-------------------------------------+--------------------------------------------------+
    | spectral_angle_y_ions               | Normalized spectral contrast angle (SA) on       |
    |                                     | y-ions                                           |
    +-------------------------------------+--------------------------------------------------+
    | pearson_corr_single_charge          | Pearson correlation on singly charged ions       |
    +-------------------------------------+--------------------------------------------------+
    | pearson_corr_double_charge          | Pearson correlation on doubly charged ions       |
    +-------------------------------------+--------------------------------------------------+
    | pearson_corr_triple_charge          | Pearson correlation on triply charged ions       |
    +-------------------------------------+--------------------------------------------------+
    | pearson_corr_b_ions                 | Pearson correlation on b-ions                    |
    +-------------------------------------+--------------------------------------------------+
    | pearson_corr_y_ions                 | Pearson correlation on y-ions                    |
    +-------------------------------------+--------------------------------------------------+
    | spearman_corr_single_charge         | Spearman correlation on singly charged ions      |
    +-------------------------------------+--------------------------------------------------+
    | spearman_corr_double_charge         | Spearman correlation on doubly charged ions      |
    +-------------------------------------+--------------------------------------------------+
    | spearman_corr_triple_charge         | Spearman correlation on triply charged ions      |
    +-------------------------------------+--------------------------------------------------+
    | spearman_corr_b_ions                | Spearman correlation on b-ions                   |
    +-------------------------------------+--------------------------------------------------+
    | spearman_corr_y_ions                | Spearman correlation on y-ions                   |
    +-------------------------------------+--------------------------------------------------+

Single-cell features (``sc_features``)
--------------------------------------

These features are **off by default**. They are added when the ``sc_features`` config option is set,
and implied by ``allFeatures``. They were developed on low-input (single-cell) TMT data, where the
standard Prosit features leave a lot of the target / decoy separation on the table, but nothing about
them is specific to that setting.

The reporter-ion features are only defined for TMT data and are NaN otherwise; the channels are
matched using the run's own ``massTolerance`` when that is given in ppm, and a 20 ppm window
otherwise. Note that they let quantification influence identification; that is defensible as a
spectrum-quality measure, but it has to be stated whenever they are used.

.. table::
    :class: fixed-table

    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Feature                 | Description                                                                                                                                                        |
    +=========================+====================================================================================================================================================================+
    | mean_ppm_error          | Mean mass deviation of the matched fragments, in ppm                                                                                                               |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | max_ppm_error           | Largest mass deviation of the matched fragments, in ppm                                                                                                            |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | std_ppm_error           | Standard deviation of the matched fragments' mass deviation, in ppm                                                                                                |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | intensity_coverage      | (summed intensity of matched fragments) / (summed intensity of all observed peaks)                                                                                 |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | annotated_frac_avg20    | (number of matched peaks) / (matched + competing peaks), where a competing peak is an unmatched observed peak with at least 20% of the mean matched-peak intensity |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | annotated_frac_all      | (number of matched peaks) / (number of observed peaks), with no significance test on the unmatched ones                                                            |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | longest_b_series        | Longest run of consecutive b-ion positions with at least one matched charge state                                                                                  |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | longest_y_series        | Longest run of consecutive y-ion positions with at least one matched charge state                                                                                  |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | longest_series_frac     | max(longest_b_series, longest_y_series) / (peptide length - 1); percolator is a linear model and cannot form this ratio from the parts                             |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | n_reporter_channels     | Number of the 11 TMT reporter channels carrying a peak; NaN for non-TMT data                                                                                       |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | reporter_intensity_frac | (summed reporter-ion intensity) / (summed intensity of all observed peaks); NaN for non-TMT data                                                                   |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | reporter_max_frac       | (largest single reporter channel) / (summed reporter-ion intensity), i.e. carrier dominance; NaN for non-TMT data                                                  |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | spectral_angle_no_b1    | Normalized spectral contrast angle with the b1 ions masked out; b1 is rarely observed under HCD, so its absence should not penalise a correct PSM                  |
    +-------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+
