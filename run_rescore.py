from oktoberfest.ce_calibration import CeCalibration
from oktoberfest.re_score import ReScore
import datetime
re_score_raw = ReScore(search_path="/media/kusterlab/internal_projects/active/ProteomeTools/TMT/datasets/iTRAQ_8_021401/combined/txt/msms.txt",
                       raw_path="/media/kusterlab/internal_projects/active/ProteomeTools/TMT/datasets/iTRAQ_8_021401/",
                       out_path="/media/kusterlab/internal_projects/active/ProteomeTools/TMT/datasets/iTRAQ_8_021401/labeled_all_feat/",
                       mzml_reader_package='pyteomics')

re_score_raw.get_raw_files()
re_score_raw.split_msms()
re_score_raw.calculate_features()
re_score_raw.merge_input('prosit')
re_score_raw.merge_input('andromeda')

re_score_raw.rescore_with_perc('prosit')
re_score_raw.rescore_with_perc('andromeda')

