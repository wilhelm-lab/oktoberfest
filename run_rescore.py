from oktoberfest.re_score import ReScore

re_score_raw = ReScore(
    search_path="/media/kusterlab/internal_projects/active/microFlAIMS/06_Fractionation_LF/txt_100per/msms.txt",
    raw_path="/media/kusterlab/internal_projects/active/microFlAIMS/06_Fractionation_LF/",
    out_path="/media/kusterlab/internal_projects/active/microFlAIMS/06_Fractionation_LF/test_oktoberfest/",
    mzml_reader_package="pyteomics",
)

re_score_raw.get_raw_files()
re_score_raw.split_msms()
re_score_raw.calculate_features()
re_score_raw.merge_input("prosit")
# re_score_raw.merge_input('andromeda')

re_score_raw.rescore_with_perc("prosit")
# re_score_raw.rescore_with_perc('andromeda')
