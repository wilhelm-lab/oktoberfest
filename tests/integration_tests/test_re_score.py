import os

from oktoberfest.re_score import ReScore


def integration_test():
    """Integration test for the rescoring."""
    search_dir = "data/plasma"
    msms_path = os.path.join(search_dir, "msms.txt")
    config_path = os.path.join(search_dir, "config.json")
    re_score = ReScore(
        search_path=msms_path, raw_path=search_dir, out_path=os.path.join(search_dir, "out"), config_path=config_path
    )
    re_score.get_raw_files()
    re_score.split_msms()
    re_score.calculate_features()
    re_score.merge_input()
    re_score.merge_input('andromeda')
    re_score.rescore_with_perc()
    re_score.rescore_with_perc('andromeda')


if __name__ == "__main__":
    integration_test()
