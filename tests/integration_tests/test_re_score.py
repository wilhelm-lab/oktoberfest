from pathlib import Path

from oktoberfest.re_score import ReScore


def integration_test():
    """Integration test for the rescoring."""
    search_dir = Path("data/plasma")
    msms_path = search_dir / "msms.txt"
    config_path = search_dir / "config.json"
    out_path = search_dir / "out"
    out_path.mkdir(exist_ok=True)
    re_score = ReScore(search_path=msms_path, raw_path=search_dir, out_path=out_path, config_path=config_path)
    re_score.get_raw_files()
    re_score.split_msms()
    re_score.calculate_features()
    re_score.merge_input()
    re_score.rescore()


if __name__ == "__main__":
    integration_test()
