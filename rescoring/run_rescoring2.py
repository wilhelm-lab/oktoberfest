from oktoberfest.runner import run_job
import argparse
import json
import os

"""
Code to re-run rescoring on all files in /cmnfs/home/m.khanh/rescoring_output on server
"""

input_dir = '/cmnfs/data/proteomics/single_cell/PXD043355/'
output_dir = '/cmnfs/home/m.khanh/rescoring_output'

dir = [f for f in os.listdir(output_dir) if os.path.isdir(f)]

for d in dir:
    run_type = d.split("_")[-1].strip('\n')
    inp = d[1:-1]
    out = output_dir + d

    if run_type == 'sqrt':
        task_config_rescoring = {
            "type": "Rescoring",
            "tag": "",
            "inputs": {
                "search_results": inp + "/txt/msms.txt",
                "search_results_type": "Maxquant",
                "spectra": inp,
                "spectra_type": "raw"
            },
            "output": out,
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
    elif run_type == 'sum':
        task_config_rescoring = {
            "type": "Rescoring",
            "tag": "",
            "inputs": {
                "search_results": inp + "/txt/msms.txt",
                "search_results_type": "Maxquant",
                "spectra": inp,
                "spectra_type": "raw"
            },
            "output": out,
            "models": {
                "intensity": "Prosit_2023_intensity_sum",
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
    elif run_type == 'basic':
        task_config_rescoring = {
            "type": "Rescoring",
            "tag": "",
            "inputs": {
                "search_results": inp + "/txt/msms.txt",
                "search_results_type": "Maxquant",
                "spectra": inp,
                "spectra_type": "raw"
            },
            "output": out,
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
    with open(out + "/rescoring_config.json", 'w') as fp:
        json.dump(task_config_rescoring, fp)

    run_job(out + "/rescoring_config.json")
