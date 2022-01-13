from typing import Dict
import json
import logging

logger = logging.getLogger(__name__)


class Config:
    data: Dict
    
    def __init__(self):
        self.data = {}
    
    def read(self, config_path: str):
        logger.info(f"Reading configuration from {config_path}")
        with open(config_path) as f:
            self.data = json.load(f)
    
    def get_num_threads(self):
        if "numThreads" in self.data:
            return self.data["numThreads"]
        else:
            return 1
    
    def get_fasta(self):
        if "fileUploads" in self.data:
            return self.data['fileUploads']['fasta']
        elif "uploads" in self.data:
            return self.data['uploads']['FASTA']
        else:
            raise ValueError("No fasta file specified in config file")
    
    def get_prosit_server(self):
        return self.data['prosit_server']
    
    def get_models(self):
        if "models" in self.data:
            return self.data['models']
        elif "selectedIntensityModel" in self.data:
            return { "selectedIntensityModel" : self.data["selectedIntensityModel"],
                     "selectedIRTModel": self.data["selectedIRTModel"],
                     "selectedProteotypicityModel": self.data["selectedProteotypicityModel"] }
    def get_tag(self):
        if "tag" in self.data:
            return self.data['tag']
        else:
            return "tmt"

    def get_all_features(self):
        if "all_features" in self.data:
            return self.data['all_features']
        else:
            return False

    def get_job_type(self):
        if "jobType" in self.data:
            return self.data['jobType']
        elif "type" in self.data:
            return self.data['type']
        else:
            raise ValueError("No job type specified in config file")
    
    def get_raw_type(self):
        if "fileUploads" in self.data: # thermo or mzml
            return self.data["fileUploads"]["raw_type"]
        else:
            return "thermo"
    
    def get_search_type(self):
        if "fileUploads" in self.data: # maxquant or internal
            return self.data["fileUploads"]["search_type"]
        else:
            return "maxquant"
    
