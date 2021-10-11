import os
import sys

from .re_score import ReScore

__version__ = "0.1.0"
__copyright__ = '''Copyright (c) 2020-2021 Oktoberfest dev-team. All rights reserved.
Written by 
- Wassim Gabriel (wassim.gabriel@tum.de),
- Ludwig Lautenbacher (ludwig.lautenbacher@tum.de),
- Matthew The (matthew.the@tum.de),
- Firas Hamood (firas.hamood@tum.de),
- Cecilia Jensen (cecilia.jensen@tum.de)
at the Technical University of Munich.'''


def main():
    print('Oktoberfest version %s\n%s' % (__version__, __copyright__))
    print('Issued command:', os.path.basename(__file__) + " " + " ".join(map(str, sys.argv[1:])))
    
    args = parse_args()  

    run_oktoberfest(args.search_dir, args.config_path)


def parse_args():
    import argparse
    apars = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('search_dir', default=None, metavar = "S",
                       help='''Directory containing the msms.txt and raw files
                            ''')
    
    apars.add_argument('--config_path', default=None, metavar = "C",
                       help='''Path to config file in json format. If this argument is not specified, we try to find and use a file called config.json in <search_dir>.
                            ''')
                            
    args = apars.parse_args()
    
    return args


def run_oktoberfest(search_dir, config_path):
    msms_path = os.path.join(search_dir, "msms.txt")
    if not config_path:
        config_path = os.path.join(search_dir, "config.json")
    
    re_score = ReScore(search_path=msms_path,
                       raw_path=search_dir)
    re_score.get_raw_files()
    re_score.split_msms()
    re_score.calculate_features()
    re_score.merge_input()
    re_score.rescore_with_perc()
    
