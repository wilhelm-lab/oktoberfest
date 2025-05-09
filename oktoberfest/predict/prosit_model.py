#imports
import torch
import torch.nn as nn
import torch.nn.functional as F
#from model_scripts.prosit.feature_extractor import PTM_LOSS_LOOKUP, PTM_MOD_DELTA_MASS_LOOKUP, PTM_GAIN_LOOKUP, PTM_ATOM_COUNT_LOOKUP, PTM_RED_SMILES_LOOKUP
import numpy as np

def masked_spectral_distance(y_pred, y_true, reduction='mean'):
    """
    Calculates the masked spectral distance between true and predicted intensity vectors.
    The masked spectral distance is a metric for comparing the similarity between two intensity vectors.

    Masked. normalized spectral angles between true and predicted vectors
    > arccos(1*1 + 0*0) = 0 -> SL = 0 -> high correlation
    > arccos(0*1 + 1*0) = pi/2 -> SL = 1 -> low correlation
    :param y_true: tensor containing the true values with shape (batch_size, num_values)
    :param y_pred: tensor containing the predicted values with the same shape as y_true
    :param reduction: specifies the reduction method used on the output; Pytorch can only implicitly create a gradient from a scalar; options: 'mean', 'sum' 'none'
    """

    # define epsilon to avoid numerical instability
    epsilon = 1e-7

    pred_masked = ((y_true + 1) * y_pred) / (y_true + 1 + epsilon)
    true_masked = ((y_true + 1) * y_true) / (y_true + 1 + epsilon)

    pred_norm = torch.nn.functional.normalize(pred_masked, p=2, dim=-1)
    true_norm = torch.nn.functional.normalize(true_masked, p=2, dim=-1)

    product = torch.sum(pred_norm * true_norm, dim=1)
    arccos = torch.acos(torch.clamp(product, -1 + epsilon, 1 - epsilon))

    batch_values = 2 * arccos / np.pi

    if reduction == "mean":
        return torch.mean(batch_values)
    elif reduction == "sum":
        return torch.sum(batch_values)
    elif reduction == "none":
        return batch_values
    else:
        raise NotImplementedError(f"Reduction Method {reduction} is not implemented. Available reductions are 'mean', 'sum' and 'none'.")

# Global constants
_ALPHABET_UNMOD_ORDERED = "ACDEFGHIKLMNPQRSTVWY"
ALPHABET_UNMOD = {k: v for v, k in enumerate(_ALPHABET_UNMOD_ORDERED, start=1)}
ALPHABET_UNMOD.update({"[]-": 21, "-[]": 22})
ALPHABET_UNMOD.update({"[UNIMOD:737]-": 21})

"""
# define the feature extractor
FEATURE_EXTRACTORS_PARAMETERS = {
    "mod_loss": {
        # Default value is 1 for all SIX atoms
        "feature_default_value": [1] * 6,
        "lookup_table": PTM_LOSS_LOOKUP,
        "description": "Loss of atoms due to PTM.",
    },
    "delta_mass": {
        # Default value is 0 for the whole sequence
        "feature_default_value": 0,
        "lookup_table": PTM_MOD_DELTA_MASS_LOOKUP,
        "description": "Delta mass of PTM.",
    },
    "mod_gain": {
        # Default value is 1 for all SIX atoms
        "feature_default_value": [1] * 6,
        "lookup_table": PTM_GAIN_LOOKUP,
        "description": "Gain of atoms due to PTM.",
    },
    "atom_count": {
        # Default value is 1 for all SIX atoms
        "feature_default_value": [1] * 6,
        "lookup_table": PTM_ATOM_COUNT_LOOKUP,
        "description": "Atom count of PTM.",
    },
    "red_smiles": {
        # Default value is 0 for the whole PTM smiles representation (currently 60)
        "feature_default_value": [0] * 60,
        "lookup_table": PTM_RED_SMILES_LOOKUP,
        "description": "Reduced SMILES representation of PTM.",
    },
}
"""

# Custom Attention Layer
class TorchAttentionLayer(nn.Module):
    def __init__(self, input_shape, context=False, bias=True):
        super(TorchAttentionLayer, self).__init__()
        self.context = context
        self.bias = bias
        self.W = nn.Parameter(torch.Tensor(input_shape[-1]))
        self.u = nn.Parameter(torch.Tensor(input_shape[-1])) if context else None
        self.b = nn.Parameter(torch.zeros(input_shape[1])) if bias else None
        nn.init.normal_(self.W)
        if context:
            nn.init.normal_(self.u)

    def forward(self, x, mask=None):
        a = torch.tanh(torch.matmul(x, self.W))
        if self.bias:
            a = a + self.b
        if self.context:
            a = torch.matmul(x, self.u)
        a = torch.exp(a)
        if mask is not None:
            a = a * mask.float()
        tmp = torch.sum(a, dim=1, keepdim=True) + 1e-10
        a = a / tmp
        weighted_input = x * a.unsqueeze(-1)
        return weighted_input.sum(dim=1)

# custom Decoder Attention Layer
class TorchDecoderAttentionLayer(nn.Module):
    def __init__(self, time_steps):
        super(TorchDecoderAttentionLayer, self).__init__()
        self.time_steps = time_steps
        self.dense = nn.Linear(time_steps, time_steps, bias=True)
        self.softmax = nn.Softmax(dim=-1)

    def forward(self, input):
        x = input.transpose(1,2)
        x = self.dense(x)
        x = self.softmax(x)
        x = x.transpose(1,2)
        out = input * x
        return out

# Prosit Retention Time
class TorchPrositRetentionTimePredictor(nn.Module):
    def __init__(self, embedding_output_dim=16, seq_length=32, alphabet=ALPHABET_UNMOD, dropout_rate=0.5, latent_dropout_rate=0.1, recurrent_layers_sizes=(256, 512), regressor_layer_size=512):
        super(TorchPrositRetentionTimePredictor, self).__init__()
        self.embedding = nn.Embedding(len(alphabet), embedding_output_dim)

        self.sequence_encoder = SequenceEncoder(embedding_output_dim, recurrent_layers_sizes[0], recurrent_layers_sizes[1], dropout_rate)

        self.attention = TorchAttentionLayer(input_shape=(1, seq_length, recurrent_layers_sizes[1]))

        self.regressor = nn.Sequential(
            nn.Linear(regressor_layer_size, regressor_layer_size),
            nn.ReLU(),
            nn.Dropout(latent_dropout_rate)
        )
        self.output_layer = nn.Linear(regressor_layer_size, 1)

    def forward(self, inputs, mode):
        # the mode parameter is not needed and only exists to keep the forward call signature of other models
        x = self.embedding(inputs)
        x = self.sequence_encoder(x)
        x = self.attention(x)
        x = self.regressor(x)
        x = self.output_layer(x)
        return x

class MetaEncoder(nn.Module):
    def __init__(self, input_size, recurrent_layer_size, dropout_rate):
        super(MetaEncoder, self).__init__()
        self.dense = nn.Linear(input_size, recurrent_layer_size)
        self.dropout = nn.Dropout(dropout_rate)

    def forward(self, inputs):
        x = torch.cat(inputs, dim=-1)
        x = self.dense(x)
        x = self.dropout(x)
        return x

class MetaEncoder2(nn.Module):
    def __init__(self, 
        recurrent_layer_size,
        dropout_rate,
        use_charge=True, 
        use_energy=False, 
        use_method=False, 
        max_charge=None, 
        num_methods=None,
    ):
        super(MetaEncoder2, self).__init__()
        self.use_charge = use_charge
        self.use_energy = use_energy
        self.use_method = use_method
        self.max_charge = max_charge
        self.num_methods = num_methods
        input_size = max_charge if self.use_charge else 0
        input_size += 1 if self.use_energy else 0
        input_size += num_methods if self.use_method else 0
        self.meta_encoder = MetaEncoder(input_size, recurrent_layer_size, dropout_rate)

    def forward(self, charge=None, energy=None, method=None):
        inputs = []
        if self.use_charge:
            charge_vector = F.one_hot(charge.long(), self.max_charge).float() # bs, max_charge
            inputs.append(charge_vector)
        if self.use_energy:
            energy_vector = energy[:,None]
            inputs.append(energy_vector) # bs, 1
        if self.use_method:
            method_vector = F.one_hot(method.long(), self.num_methods).float() # bs, num_methods
            inputs.append(method_vector)
        return self.meta_encoder(inputs)


class MetaDataFusionLayer(nn.Module):
    def __init__(self, max_ion):
        super(MetaDataFusionLayer, self).__init__()
        self.max_ion = max_ion

    def forward(self, x, encoded_meta):
        x = x * encoded_meta
        x = x.unsqueeze(1).repeat(1, self.max_ion, 1)
        return x

class MetaDataFusionLayer2(nn.Module):
    def __init__(self, *args):
        super(MetaDataFusionLayer2, self).__init__()

    def forward(self, x, encoded_meta):
        # x: bs, sl, units
        # encoded_meta: bs, units
        x = x * encoded_meta[:, None]
        return x

class SequenceEncoder(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, dropout_rate):
        super(SequenceEncoder, self).__init__()
        self.encoder1 = nn.GRU(input_size, hidden_size, batch_first=True, bidirectional=True)
        self.encoder2 = nn.GRU(hidden_size*2, output_size, batch_first=True)
        self.dropout = nn.Dropout(dropout_rate)

    def forward(self, x):
        x, _ = self.encoder1(x)
        x = self.dropout(x)
        x, _ = self.encoder2(x)
        x = self.dropout(x)
        return x

class Decoder(nn.Module):
    def __init__(self, regressor_layer_size, seq_length, dropout_rate):
        super(Decoder, self).__init__()
        self.gru = nn.GRU(regressor_layer_size, regressor_layer_size, batch_first=True)
        self.dropout = nn.Dropout(dropout_rate)
        self.decoder_attention = TorchDecoderAttentionLayer(seq_length)

    def forward(self, x):
        x, _ = self.gru(x)
        x = self.dropout(x)
        x = self.decoder_attention(x)
        return x

class Regressor(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(Regressor, self).__init__()
        self.dense = nn.Linear(input_size, hidden_size)
        self.activation = nn.LeakyReLU()

    def forward(self, x):
        x = self.dense(x)
        x = self.activation(x)
        x = x.flatten(start_dim=-2)
        return x

class Regressor2(nn.Module):
    def __init__(self, input_size, final_units):
        super(Regressor2, self).__init__()
        self.dense = nn.Linear(input_size, final_units)
        self.activation = nn.LeakyReLU()
    
    def forward(self, x):
        x = self.dense(x)
        x = self.activation(x)
        return x.mean(1)

# Prosit Intensity Predictor
class TorchPrositIntensityPredictor(nn.Module):
    DEFAULT_INPUT_KEYS = {
        "SEQUENCE_KEY": "sequence",
        "COLLISION_ENERGY_KEY": "collision_energy_aligned_normed",
        "PRECURSOR_CHARGE_KEY": "precursor_charge_onehot",
    }

    # can be extended to include all possible meta data
    META_DATA_KEYS = {
        "COLLISION_ENERGY_KEY": "collision_energy_aligned_normed",
        "PRECURSOR_CHARGE_KEY": "precursor_charge_onehot",
    }

    # retrieve the Lookup PTM feature keys
    #PTM_INPUT_KEYS = [*FEATURE_EXTRACTORS_PARAMETERS.keys()]

    def __init__(
        self, 
        embedding_output_dim=16, 
        seq_length=30, 
        final_units=174, # len_fion=6
        tokens=20, # alphabet=ALPHABET_UNMOD
        use_charge=True,
        use_energy=False,
        use_method=True,
        max_charge=6,
        num_methods=4,
        dropout_rate=0.2, 
        latent_dropout_rate=0.1, 
        recurrent_layers_sizes=(256, 512),
        regressor_layer_size=512, 
        use_prosit_ptm_features=False, 
        with_termini=True,
        **kwargs
    ):
        super(TorchPrositIntensityPredictor, self).__init__()
        self.dropout_rate = dropout_rate
        self.latent_dropout_rate = latent_dropout_rate
        self.embedding_output_dim = embedding_output_dim
        self.embeddings_count = tokens

        self.use_prosit_ptm_features = use_prosit_ptm_features

        self.max_ion = seq_length - 3 if with_termini else seq_length - 1
        self.atleast1 = use_charge or use_energy or use_method

        self.embedding = nn.Embedding(self.embeddings_count, self.embedding_output_dim)

        self.sequence_encoder = SequenceEncoder(embedding_output_dim, recurrent_layers_sizes[0], recurrent_layers_sizes[1], dropout_rate)

        #TODO ptm encoders that are not relavant to this imidiate project
        self.meta_encoder = MetaEncoder2(
            recurrent_layers_sizes[1], 
            self.dropout_rate,
            use_charge=use_charge,
            use_energy=use_energy,
            use_method=use_method,
            max_charge=max_charge+1,
            num_methods=num_methods,
        ) if self.atleast1 else None

        #self.attention = TorchAttentionLayer(input_shape=(1, seq_length, recurrent_layers_sizes[1]))

        self.meta_data_fusion_layer = MetaDataFusionLayer2(self.max_ion)

        self.decoder = Decoder(regressor_layer_size, seq_length, dropout_rate)

        self.regressor = Regressor2(regressor_layer_size, final_units)

        self.global_step = nn.Parameter(torch.tensor(0), requires_grad=False)

    def forward(self, intseq, charge=None, energy=None, method=None):
        
        x = self.embedding(intseq)
        x = self.sequence_encoder(x)
        #x = self.attention(x)
        
        encoded_meta = self.meta_encoder(charge, energy, method) if self.atleast1 else None
        if self.meta_data_fusion_layer and encoded_meta is not None:
            x = self.meta_data_fusion_layer(x, encoded_meta)
        x = self.decoder(x)
        x = self.regressor(x)
        return x

    def _collect_values_from_inputs_if_exists(self, inputs, keys_mapping):
        collected_values = []

        keys = []
        if isinstance(keys_mapping, dict):
            keys = keys_mapping.values()

        elif isinstance(keys_mapping, list):
            keys = keys_mapping

        for key_in_inputs in keys:
            # get the input under the specified key if exists
            single_input = inputs.get(key_in_inputs, None)
            if single_input is not None:
                if single_input.ndim == 1:
                    single_input = single_input.unsqueeze(1)
                collected_values.append(single_input)
        return collected_values


def PrositModel(
    final_units,
    tokens,
    max_charge,
    kwargs
):
    return TorchPrositIntensityPredictor(
        final_units=final_units,
        tokens=tokens,
        max_charge=max_charge,
        **kwargs
    )

if __name__ == "__main__":
    # prepare dummy input data and labels
    sequences = torch.randint(low=0, high=len(ALPHABET_UNMOD), size=(5,30))

    random_charges = torch.randint(low=0, high=6, size=(5,))
    onehot_charges = F.one_hot(random_charges, num_classes=6).float()
    random_ce = torch.rand(5)

    random_labels = torch.rand((5, 174))
    mask = torch.rand((5, 174)) < 0.25
    random_labels[mask] = -1

    input_dict = {
        'intseq': sequences,
        'energy': random_ce,
        'charge': random_charges
    }

    # define model, lossfunction and optimizer
    model = PrositModel(
        final_units=174,
        tokens=len(ALPHABET_UNMOD),
        max_charge=6,
        kwargs={
            'use_charge': True,
            'use_energy': True,
            'use_method': False
        },
    )
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4)
    output = model(**input_dict)
    loss = masked_spectral_distance(output, random_labels)
    optimizer.step()
