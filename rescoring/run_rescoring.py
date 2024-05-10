from oktoberfest.runner import run_job
import argparse
import json
import os

"""
Code to run rescoring on server
"""

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i')
parser.add_argument('--output', '-o')
parser.add_argument('--type', '-t', choices=['sqrt', 'sum', 'basic', 'both'])

args = parser.parse_args()

task_config_rescoring = {}

if args.type == 'sqrt':
    task_config_rescoring = {
        "type": "Rescoring",
        "tag": "",
        "inputs": {
            "search_results": args.input + "/txt/msms.txt",
            "search_results_type": "Maxquant",
            "spectra": args.input,
            "spectra_type": "raw"
        },
        "output": args.output,
        "models": {
            "intensity": "Prosit_2023_intensity_sqrt",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "10.157.98.62:9500",
        "ssl": False,
        "thermoExe": "/cmnfs/home/m.khanh/oktoberfest/rescoring/ThermoRawFileParser1.4.3/ThermoRawFileParser.exe",
        "numThreads": 1,
        "fdr_estimation_method": "mokapot",
        "regressionMethod": "spline",
        "allFeatures": False,
        "massTolerance": 20,
        "unitMassTolerance": "ppm"
    }
elif args.type == 'sum':
    task_config_rescoring = {
        "type": "Rescoring",
        "tag": "",
        "inputs": {
            "search_results": args.input + "/txt/msms.txt",
            "search_results_type": "Maxquant",
            "spectra": args.input,
            "spectra_type": "raw"
        },
        "output": args.output,
        "models": {
            "intensity": "Prosit_2024_intensity_sum", #"Prosit_2023_intensity_sum",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "10.157.98.62:9500",
        "ssl": False,
        "thermoExe": "/cmnfs/home/m.khanh/oktoberfest/rescoring/ThermoRawFileParser1.4.3/ThermoRawFileParser.exe",
        "numThreads": 1,
        "fdr_estimation_method": "mokapot",
        "regressionMethod": "spline",
        "allFeatures": False,
        "massTolerance": 20,
        "unitMassTolerance": "ppm"
    }
elif args.type == 'basic':
    task_config_rescoring = {
        "type": "Rescoring",
        "tag": "",
        "inputs": {
            "search_results": args.input + "/txt/msms.txt",
            "search_results_type": "Maxquant",
            "spectra": args.input,
            "spectra_type": "raw"
        },
        "output": args.output,
        "models": {
            "intensity": "Prosit_2020_intensity_HCD",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "koina.proteomicsdb.org:443",
        "ssl": True,
        "thermoExe": "/cmnfs/home/m.khanh/oktoberfest/rescoring/ThermoRawFileParser1.4.3/ThermoRawFileParser.exe",
        "numThreads": 1,
        "fdr_estimation_method": "mokapot",
        "regressionMethod": "spline",
        "allFeatures": False,
        "massTolerance": 20,
        "unitMassTolerance": "ppm"
    }
elif args.type == 'both':
    task_config_rescoring = {
        "type": "Rescoring",
        "tag": "",
        "inputs": {
            "search_results": args.input + "/txt/msms.txt",
            "search_results_type": "Maxquant",
            "spectra": args.input,
            "spectra_type": "raw"
        },
        "output": args.output,
        "models": {
            "intensity": "Prosit_2024_intensity_single_cell",
            "irt": "Prosit_2019_irt"
        },
        "prediction_server": "10.157.98.62:9500",
        "ssl": False,
        "thermoExe": "/cmnfs/home/m.khanh/oktoberfest/rescoring/ThermoRawFileParser1.4.3/ThermoRawFileParser.exe",
        "numThreads": 1,
        "fdr_estimation_method": "mokapot",
        "regressionMethod": "spline",
        "allFeatures": False,
        "massTolerance": 20,
        "unitMassTolerance": "ppm"
    }

if not os.path.isdir(args.output):
    os.mkdir(args.output)

with open(args.output + "/rescoring_config.json", 'w') as fp:
    json.dump(task_config_rescoring, fp)

run_job(args.output + "/rescoring_config.json")
