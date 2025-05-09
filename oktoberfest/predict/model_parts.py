import torch as th
from torch import nn
import numpy as np
twopi = 2*np.pi

class BatchTorch1d(nn.Module):
    def __init__(self, units):
        super(BatchTorch1d, self).__init__()
        self.norm = nn.BatchNorm1d(units)
    def forward(self, x):
        return self.norm(x.transpose(-1,-2)).transpose(-1,-2)

def get_norm_type(string):
    if string.lower() == 'layer':
        return nn.LayerNorm
    elif string.lower() == 'batch':
        return BatchTorch1d

class QKVAttention(nn.Module):
    def __init__(self,
                 heads,
                 dim,
                 sl=None,
                 is_relpos=False, 
                 max_rel_dist=None,
    ):
        super(QKVAttention, self).__init__()
        self.heads = heads
        self.dim = dim
        self.sl = sl
        self.is_relpos = is_relpos
        self.maxd = max_rel_dist
        
        self.scale = dim**-0.5
        if is_relpos:
            assert sl is not None
            self.maxd = sl if max_rel_dist==None else max_rel_dist
            self.ak = self.build_relpos_tensor(sl, self.maxd) # sl, sl, dim
            self.av = self.build_relpos_tensor(sl, self.maxd) # same
    
    def build_relpos_tensor(self, seq_len, maxd=None):
        # This is the equivalent of an input dependent ij bias, rather than one
        # that is a directly learned ij tensor
        maxd = seq_len-1 if maxd==None else (maxd-1 if maxd==seq_len else maxd)
        a = th.arange(seq_len, dtype=th.int32)
        b = th.arange(seq_len, dtype=th.int32)
        relpos = a[:, None] - b[None]
        tsr = th.zeros(
            2*seq_len-1, self.dim
        ).normal_(0, seq_len**-0.5).type(th.float32)
        # set maximum distance
        relpos = relpos.clamp(-maxd, maxd)
        relpos += maxd
        relpos_tsr = tsr[relpos]
        relpos_tsr.requires_grad = True
        
        return relpos_tsr
    
    def forward(self, Q, K, V, mask=None, bias=None, gate=None, return_full=False):#, before_lambda=nn.Identity(), after_lambda=nn.Identity()):
        bsp, sl, dim = Q.shape
        # shape: batch/heads/etc, sequence_length, dim
        QK = self.scale * th.einsum('abc,adc->abd', Q, K)
        if self.is_relpos:
            QK += th.einsum('abc,bec->abe', Q, self.ak)
        QK = QK.reshape(-1, self.heads, sl, K.shape[1])
        if bias is not None:
            QK += bias

        #QK = before_lambda(QK)

        # Make mask fit 4 dimensional QK matrix
        if mask == None:
            mask = th.zeros_like(QK)
        elif len(mask.shape) == 2:
            mask = mask[:, None, None, :] # for every head and query
        elif len(mask.shape) == 3:
            mask = mask[:, None]
        
        weights = th.softmax(QK-mask, dim=-1)
        weights = weights.reshape(-1, sl, V.shape[1])
        
        #weights = after_lambda(weights)

        att = th.einsum('abc,acd->abd', weights, V)
        if self.is_relpos:
            att += th.einsum('abc,bcd->abd', weights, self.av)

        if gate is not None:
            att = att * gate

        other = [QK, weights, att] if return_full else None
        
        return att, other

class BaseAttentionLayer(nn.Module):
    def __init__(
        self, 
        indim, 
        d, 
        h, 
        out_units=None, 
        gate=False, 
        dropout=0,
        alphabet=False
    ):
        super(BaseAttentionLayer, self).__init__()
        self.indim = indim
        self.d = d
        self.h = h
        self.out_units = indim if out_units==None else out_units
        self.drop = nn.Dropout(dropout) if dropout > 0 else nn.Identity()
        self.alphabet = alphabet

        self.attention_layer = QKVAttention(h, d)
        
        shape = (d*h, self.out_units)
        self.Wo = nn.Linear(*shape, bias=True)
        #self.Wo.weight = nn.Parameter( 
        #    nn.init.normal_(th.empty(self.Wo.weight.shape), 0.0, 0.3*(h*d)**-0.5)
        #)

        self.shortcut = (
            nn.Identity()
            if self.out_units == indim else
            nn.Linear(indim, out_units)
        )
        
        self.gate = gate
        if gate:
            self.Wg = nn.Linear(indim, d*h)
        
        if alphabet:
            self.alpha = nn.Parameter(th.tensor(1.), requires_grad=True)
            self.beta = nn.Parameter(th.tensor(1.), requires_grad=True)
        
class SelfAttention(BaseAttentionLayer):
    def __init__(self, 
                 indim, 
                 d, 
                 h, 
                 out_units=None, 
                 gate=False,
                 bias=False,
                 bias_in_units=None,
                 modulator=False,
                 dropout=0,
                 alphabet=False,
                 rotation_matrix=False,
                 rotation_values=None,
    ):
        super().__init__(
            indim=indim, d=d, h=h, out_units=out_units, 
            gate=gate, dropout=dropout, alphabet=alphabet
        )

        self.qkv = nn.Linear(indim, 3*d*h, bias=True)
        self.rotation_matrix = rotation_matrix

        """
        self.before = nn.Linear(h, h, bias=False)
        self.before.weight = nn.init.normal_(self.before.weight, 0, 0.01)
        self.after = nn.Linear(h, h, bias=False)
        self.after.weight = nn.init.normal_(self.after.weight, 0, 0.01)
        self.before_lambda = lambda x: (
            self.before(x.permute([0,2,3,1])).permute([0,3,1,2])
        )
        self.after_lambda = lambda x: (
            self.after(x.reshape(-1,h,x.shape[1],x.shape[2]).permute([0,2,3,1])).
            permute([0,3,1,2]).reshape(-1,x.shape[1],x.shape[2])
        )
        """

        self.bias = bias
        if bias == 'pairwise':
            self.Wpw = nn.Linear(bias_in_units, h)
        elif bias == 'regular':
            self.Wb = nn.Linear(indim, h)
        
        self.modulator = modulator
        if modulator:
            self.alphaq = nn.Parameter(th.tensor(0.))
            self.alphak = nn.Parameter(th.tensor(0.))
            self.alphav = nn.Parameter(th.tensor(0.))

        if rotation_matrix:
            if type(rotation_values)==int:
                values = th.arange(rotation_values, dtype=th.float32)[None]
                min_lambda = 1
                max_lambda = 1000
            else:
                values = rotation_values
                min_lambda = 0.001
                max_lambda = 10000
            rcos, rsin = RotationMatrix(
                values=values,
                units=d,
                min_lam=min_lambda,
                max_lam=max_lambda,
            )
            self.rcos = nn.Parameter(rcos, requires_grad=False)
            self.rsin = nn.Parameter(rsin, requires_grad=False)

        
    def get_qkv(self, qkv):
        bs, sl, units = qkv.shape
        Q, K, V = qkv.split(units//3, -1)
        Q = Q.reshape(-1, sl, self.d, self.h)
        Q = Q.permute([0,3,1,2]).reshape(-1, sl, self.d)
        K = K.reshape(-1, sl, self.d, self.h)
        K = K.permute([0,3,1,2]).reshape(-1, sl, self.d)
        V = V.reshape(-1, sl, self.d, self.h)
        V = V.permute([0,3,1,2]).reshape(-1, sl, self.d)
        if self.modulator:
            Q *= th.sigmoid(self.alphaq)
            K *= th.sigmoid(self.alphak)
            V *= th.sigmoid(self.alphav)
        
        if self.rotation_matrix:
            cos = self.rcos[:,:sl].tile([bs*self.h, 1, 1])
            sin = self.rsin[:,:sl].tile([bs*self.h, 1, 1])
            Q_ = Q*cos
            Q_[..., ::2] += -sin[..., ::2]*Q_[..., 1::2]
            Q_[..., 1::2] += sin[..., ::2]*Q_[..., ::2]
            Q = Q_
            K_ = K*cos
            K_[..., ::2] += -sin[..., ::2]*K_[..., 1::2]
            K_[..., 1::2] += sin[..., ::2]*K_[..., ::2]
            K = K_
        
        return Q, K, V # bs*h, sl, d
    
    def forward(self, x, mask=None, biastsr=None, return_full=False):
        bs, sl, units = x.shape
        qkv = self.qkv(x) # bs, sl, 3*d*h
        Q, K, V = self.get_qkv(qkv) # bs*h, sl, d
        if self.bias == 'regular':
            B = self.Wb(x)[:,None]
            B = B.permute([0,3,1,2])
        elif self.bias == 'pairwise':
            B = self.Wpw(biastsr) # bs, sl, sl, h
            B = B.permute([0,3,1,2]) # bs, h, sl, sl
        else:
            B = None
        if self.gate:
            G = th.sigmoid(self.Wg(x)) # bs, sl, d*h
            G = (
                G.reshape(bs, sl, self.d, self.h)
                .permute([0,3,1,2])
                .reshape(bs*self.h, sl, self.d)
            )
        else:
            G = None
        att, other = self.attention_layer(Q, K, V, mask, bias=B, gate=G, return_full=return_full)#, before_lambda=self.before_lambda, after_lambda=self.after_lambda) # bs*h, sl, d
        att = att.reshape(-1, self.h, sl, self.d)
        att = att.permute([0,2,3,1])
        att = att.reshape(-1, sl, self.d*self.h) # bs, sl, d*h
        resid = self.Wo(att)
        
        if self.alphabet:
            output = self.alpha*self.shortcut(x) + self.beta*self.drop(resid)
        else:
            output = self.shortcut(x) + self.drop(resid)
        
        other = [Q, K, V] + other + [resid] if return_full else None
        
        return {'out': output, 'other': other}

class CrossAttention(BaseAttentionLayer):
    def __init__(self, indim, kvindim, d, h, out_units=None, dropout=0, alphabet=False):
        super().__init__(
            indim=indim, d=d, h=h, out_units=out_units, 
            dropout=dropout, alphabet=alphabet
        )
        
        self.Wq = nn.Linear(indim, d*h, bias=False)
        self.Wkv = nn.Linear(kvindim, 2*d*h, bias=False)
        
    def get_qkv(self, q, kv):
        bs, sl, units = q.shape
        bs, sl2, kvunits = kv.shape
        Q = q.reshape(bs, sl, self.d, self.h)
        Q = Q.permute([0,3,1,2]).reshape(-1, sl, self.d)
        K, V = kv.split(kvunits//2, -1)
        K = K.reshape(bs, sl2, self.d, self.h)
        K = K.permute([0,3,1,2]).reshape(-1, sl2, self.d)
        V = V.reshape(bs, sl2, self.d, self.h)
        V = V.permute([0,3,1,2]).reshape(-1, sl2, self.d)
        
        return Q, K, V
    
    def forward(self, q_feats, kv_feats, mask=None):
        _, slq, _ = q_feats.shape
        Q = self.Wq(q_feats)
        KV = self.Wkv(kv_feats)
        Q, K, V = self.get_qkv(Q, KV)
        att, other = self.attention_layer(Q, K, V, mask)
        att = att.reshape(-1, self.h, slq, self.d)
        att = att.permute([0,2,1,3])
        att = att.reshape(-1, slq, self.h*self.d)
        resid = self.Wo(att)
        
        return self.shortcut(q_feats) + self.drop(resid)

class FFN(nn.Module):
    def __init__(self, indim, unit_multiplier=1, out_units=None, dropout=0, alphabet=False):
        super(FFN, self).__init__()
        self.indim = indim
        self.mult = unit_multiplier
        self.out_units = indim if out_units==None else out_units
        self.alphabet = alphabet

        self.W1 = nn.Linear(indim, indim*self.mult)
        self.W2 = nn.Linear(indim*self.mult, self.out_units, bias=False)
        shape = self.W2.weight.shape
        self.W2.weight = nn.Parameter( 
            nn.init.normal_(th.empty(shape), 0.0, 0.3*(indim*self.mult)**-0.5)
        )

        self.drop = nn.Dropout(dropout) if dropout>0 else nn.Identity()
        
        if alphabet:
            self.alpha = nn.Parameter(th.tensor(1.), requires_grad=True)
            self.beta = nn.Parameter(th.tensor(1.), requires_grad=True)
    
    def forward(self, x, embed=None, return_full=False):
        out1 = self.W1(x)
        out2 = th.relu(out1 + (0 if embed==None else embed))
        out3 = self.W2(out2)
        
        other = [out1, out3] if return_full else None
        
        if self.alphabet:
            out = self.alpha*x + self.beta*self.drop(out3)
        else:
            out = x + self.drop(out3)

        return {'out': out, 'other': other}

class TransBlock(nn.Module):
    def __init__(self,
                 attention_dict,
                 ffn_dict,
                 norm_type='layer',
                 prenorm=True,
                 embed_type=None, # preembed | ffnembed | normembed | None
                 embed_indim=256,
                 is_cross=False,
                 kvindim=256,
                 channel_alpha=False,
    ):
        super(TransBlock, self).__init__()
        self.norm_type = norm_type
        self.mult = ffn_dict['unit_multiplier']
        self.prenorm = prenorm
        self.embed_type = embed_type
        self.is_cross = is_cross

        norm = get_norm_type(norm_type)

        # How to embed, if at all, precursor level information (charge, energy, etc.)
        elementwise_affine = True
        if embed_type is not None:
            units = attention_dict['indim']
            if embed_type == 'preembed':
                num = attention_dict['indim'] if channel_alpha else 1
                self.alpha = nn.Parameter(0.1*th.ones(num), requires_grad=True)
                units = units
            elif embed_type == 'ffnembed':
                units = units * self.mult
            elif embed_type in ['normembed', 'preembed_wb']:
                units = 2 * units
                if embed_type == 'normembed':
                    elementwise_affine = False
            else:
                raise NotImplementedError("Choose a real embedding option")
        
            assert type(embed_indim) == int
            self.embed = nn.Linear(embed_indim, units)
            
        indim = attention_dict['indim']
        self.norm1 = norm(indim, elementwise_affine=elementwise_affine)
        self.norm2 = norm(ffn_dict['indim'])
        self.selfattention = SelfAttention(**attention_dict)
        if is_cross:
            cross_dict = attention_dict.copy()
            if 'pairwise_bias' in cross_dict.keys(): cross_dict.pop('pairwise_bias')
            if 'bias_in_units' in cross_dict.keys(): cross_dict.pop('bias_in_units')
            cross_dict['kvindim'] = kvindim
            self.crossnorm = norm(indim)
            self.crossattention = CrossAttention(**cross_dict)
        self.ffn = FFN(**ffn_dict)
        
        
    def forward(self, 
                x, 
                kv_feats=None, 
                embed_feats=None, 
                spec_mask=None, 
                seq_mask=None,
                biastsr=None,
                return_full=False
    ):
        bs, sl, units = x.shape
        selfmask = seq_mask if self.is_cross else spec_mask
        
        out = x
        # Embed precursor level information (?)
        if self.embed_type is not None:
            Emb = self.embed(embed_feats)[:,None,:]
    
            if self.embed_type == 'preembed':
                out = out + self.alpha[None, None] * Emb
                if self.prenorm: out = self.norm1(out)
            elif self.embed_type == 'ffnembed':
                out = self.norm1(out) if self.prenorm else out
            elif self.embed_type in ['normembed', 'preembed_wb']:
                weight, bias = Emb.split(units, -1)
                weight = 1 + weight
                bias = bias
                if self.embed_type == 'preembed_wb':
                    out = weight * out + bias
                elif self.prenorm & (self.embed_type == 'normembed'): 
                    out = weight * self.norm1(out) + bias
        else:
            if self.prenorm:
                out = self.norm1(out)
            Emb = None
        
        # Self attention
        outsa = self.selfattention(out, selfmask, biastsr, return_full=return_full)
        out = outsa['out']
        if not self.prenorm:
            out = self.norm1(out)
            if self.embed_type == 'normembed':
                out = weight * out + bias
        
        # Cross attention
        if self.is_cross:
            out = self.crossnorm(out) if self.prenorm else out
            out = self.crossattention(out, kv_feats, spec_mask)
            out = out if self.prenorm else self.crossnorm(out)
        
        # FFN
        if self.prenorm: out = self.norm2(out)
        if self.embed_type != 'ffnembed': Emb = None
        outffn = self.ffn(out, Emb, return_full=return_full)
        out = outffn['out']
        out = out if self.prenorm else self.norm2(out)
        
        other = outsa['other'] + outffn['other'] + [out] if return_full else None
        
        return {'out': out, 'other': other}

class ActModule(nn.Module):
    def __init__(self, activation):
        super(ActModule, self).__init__()
        self.act = activation
    def forward(self, x):
        return self.act(x)

def FourierFeatures(t, min_lam, max_lam, embedsz):
    x = th.arange(embedsz//2).type(th.float32).to(t.device)
    x /= embedsz//2 - 1
    """embed = ( 
        t[..., None] * 
        th.exp(-x*log(max_lam / min_lam) )[None] /
        (min_lam / 6.2831853)
    )"""
    denom = (min_lam / twopi) * (max_lam / min_lam) ** x
    embed = t[...,None] / denom[None]

    return th.cat([embed.sin(), embed.cos()], dim=-1)

def subdivide_float(x):
    mul = -1*(x < 0).type(th.float32) + (x >= 0).type(th.float32)
    X = abs(x)
    a = X.floor_divide(100)
    b = (X-a*100).floor_divide(1)
    X_ = ((X - X.floor_divide(1))*10000).round()
    c = X_.floor_divide(100)
    d = (X_ - c*100).floor_divide(1)
    
    return mul[...,None] * th.cat([a[...,None], b[...,None], c[...,None], d[...,None]], -1)

def delta_tensor(mz, shift=0):
    return mz[..., None] - mz[:, None] + shift # bs, sl, sl

def RotationMatrix(
    values: th.tensor,
    units: int,
    min_lam: float=0.001,
    max_lam: float=10000,
    dense=False,
):
    bs, sl = values.shape

    halfway = units // 2
    d = th.arange(halfway)
    Theta = (min_lam/twopi)*(max_lam/min_lam)**(2*d/units)
    pos_by_theta = values[...,None] / Theta[None,None]
    cos = np.cos(pos_by_theta)
    sin = np.sin(pos_by_theta)

    full_cos = th.tile(cos[...,None], [1,1,1,2]).reshape(cos.shape[0], cos.shape[1], -1)
    full_sin = th.tile(sin[...,None], [1,1,1,2]).reshape(sin.shape[0], sin.shape[1], -1)
    if dense:
        out = th.zeros(bs, sl, units, units)
        out[:,:,th.arange(units),th.arange(units)] = full_cos
        out[:,:,th.arange(0,units,2),th.arange(1,units,2)] = -sin
        out[:,:,th.arange(1,units,2),th.arange(0,units,2)] = sin
    else:
        out = [full_cos, full_sin]

    return out

