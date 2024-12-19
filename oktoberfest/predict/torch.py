import pandas as pd
import torch as th
import os
import yaml
import re
from tqdm import tqdm
import numpy as np
from oktoberfest.predict.peptide_encoder import PeptideEncoder

device = th.device('cuda' if th.cuda.is_available() else 'cpu')

def tokenize_modified_sequence(modseq):
    tokenized = []
    modseq = re.sub('-|(\[])', '', modseq) # remove - or []
    #modseq = re.sub('(\[]-)|(\-\[])','',modseq)

    pos = 0
    while pos < len(modseq):
        character = modseq[pos]
        hx = ord(character)
        if character == '[':
            ahead = 1
            mod = []
            while character != ']':
                mod.append(character)
                character = modseq[pos+ahead]
                ahead += 1
            token = "".join(mod) + ']'
            if pos != 0:
                tokenized[-1] += token
            else:
                tokenized.append(token)
            pos += ahead - 1
        else:
            tokenized.append(character)
        pos += 1

    return tokenized

def integerize_sequence(seq, dictionary):
    return [dictionary[aa] for aa in seq]

class TorchModel:
    def __init__(
        self,
        model_path: str,
        ion_dict_path: str,
        token_dict_path: str,
        yaml_dir_path: str,
    ):
        ion_dict = pd.read_csv(ion_dict_path, index_col='full')
        with open(os.path.join(yaml_dir_path, "model.yaml")) as f: model_config = yaml.safe_load(f)
        with open(os.path.join(yaml_dir_path, "loader.yaml")) as f: load_config = yaml.safe_load(f)
        num_tokens = len(open(token_dict_path).read().strip().split("\n")) + 1
        self.model = PeptideEncoder(
            tokens = num_tokens,
            final_units = len(ion_dict),
            max_charge = load_config['charge'][-1],
            **model_config
        )
        self.model.to(device)
        self.model.load_state_dict(th.load(model_path, map_location=device, weights_only=True))
        self.model.eval()

        self.token_dict = self.create_dictionary(token_dict_path)

        method_list = load_config['method_list']
        self.method_dic = {method: m for m, method in enumerate(method_list)}
        self.method_dicr = {n:m for m,n in self.method_dic.items()}

    def create_dictionary(self, dictionary_path):
        amod_dic = {
            line.split()[0]:m for m, line in enumerate(open(dictionary_path))
        }
        amod_dic['X'] = len(amod_dic)
        amod_dic[''] = amod_dic['X']
        return amod_dic

    def predict(self, spectra_dataframe):
        data = spectra_dataframe.obs
        
        tokenized_sequence = data.apply(
            lambda x: (
                tokenize_modified_sequence(x['MODIFIED_SEQUENCE']) + 
                ['X']*(30-x['PEPTIDE_LENGTH'])
            ), axis=1
        )
        intseq = tokenized_sequence.map(lambda x: integerize_sequence(x, self.token_dict))
        intseq = th.tensor(intseq.to_list(), dtype=th.int32)
        charge = th.tensor(data['PRECURSOR_CHARGE'].to_list(), dtype=th.int32)
        method = th.tensor(
            data['FRAGMENTATION'].map(lambda x: self.method_dic[x]).to_list(), dtype=th.int32
        )
        
        chunked_intseq = intseq.split(100,0)

        iterable = tqdm(zip(chunked_intseq, charge.split(100,0), method.split(100, 0)), total=len(chunked_intseq))
        outputs = []
        for intseq_, charge_, method_ in iterable:
            inp = {
                'intseq': intseq_.to(device),
                'charge': charge_.to(device),
                #'energy': batchdev['ce'],
                'method': method_.to(device),
            }
            with th.no_grad():
                prediction = self.model(**inp)
            prediction /= prediction.max(1, keepdim=True)[0]
            outputs.append(prediction.cpu().numpy())
        full_output = np.concatenate(outputs, axis=0)

        return {
            'intensities': full_output,
            'annotation': np.tile(spectra_dataframe.var.index.to_numpy()[None], [full_output.shape[0],1])
        }
