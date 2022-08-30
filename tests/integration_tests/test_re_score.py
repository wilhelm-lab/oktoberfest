import os

from oktoberfest.re_score import ReScore


def integration_test_1():
    ce_cal = ReScore(search_path="D:/Compmass/workDir/HCD_OT/msms.txt", raw_path="D:/Compmass/workDir/HCD_OT/")
    df_search = ce_cal._load_search()
    grouped_search = df_search.groupby("RAW_FILE")
    raw_files = grouped_search.groups.keys()
    re_score_raw = {}
    for raw_file in raw_files:
        re_score_raw[raw_file] = ReScore(
            search_path="D:/Compmass/workDir/HCD_OT/msms.txt",
            raw_path="D:/Compmass/workDir/HCD_OT/" + raw_file + ".mzml",
        )
        msms_raw = grouped_search.get_group(raw_file)
        re_score_raw[raw_file].predict_with_aligned_ce(msms_raw)
        re_score_raw[raw_file].gen_perc_metrics("prosit_" + raw_file + ".tab")

    raw_files_list = list(raw_files)
    with open("merged_prosit.tab", "wb") as fout:
        # first file:
        with open("prosit_" + raw_files_list[0] + ".tab", "rb") as f:
            fout.write(f.read())
        # now the rest:
        for raw_file in list(raw_files)[1:]:
            with open("prosit_" + rawfile + ".tab", "rb") as f:
                next(f)  # skip the header
                fout.write(f.read())


def integration_test_2():
    search_dir = "/media/processing_results/bierdimpfl/workDir/627"
    msms_path = os.path.join(search_dir, "msms.txt")
    config_path = os.path.join(search_dir, "config_new.json")
    re_score = ReScore(search_path=msms_path, raw_path=search_dir, config_path=config_path)
    re_score.get_raw_files()
    re_score.split_msms()
    re_score.calculate_features()
    re_score.merge_input()
    re_score.rescore_with_perc()


if __name__ == "__main__":
    # integration_test_1()
    integration_test_2()
