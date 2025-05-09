import numpy as np
import pandas as pd
import re

def theoretical_ions(
    ions: list, 
    max_length: int=30, 
    max_charge: int=5, 
    neutral_loss: list=[''],
    include_p=True
):
    if '' not in neutral_loss:
        neutral_loss = [''] + neutral_loss
    full_ions = [
        'p', 'p+i', 'p+2i', 'p+3i', 'p+4i',
        'p^1', 'p^1+i', 'p^1+2i', 'p^1+3i',
    ] if include_p else []
    lengths = len(full_ions)*[0]
    charges = len(full_ions)*[0]
    ion_types = len(full_ions)*['p'] if include_p else []
    NeutralLosses = len(full_ions)*['']
    for i in range(1, max_length, 1):
        for j in range(1, max_charge+1, 1):
            for k in ions:
                for l in neutral_loss:
                    charge = '' if j==1 else '^%d'%int(j)
                    nl = '' if l=='' else f'-{l}'
                    ion_types.append(k)
                    lengths.append(i)
                    charges.append(j)
                    NeutralLosses.append(l)
                    full_ions.append(f"{k}{str(i)}{nl}{charge}")
    df = pd.DataFrame({
        'ion': ion_types, 
        'length': lengths, 
        'charge': charges, 
        'neutral_loss': NeutralLosses,
    })
    df.index = pd.Index(full_ions)
    df.index.name = 'full'
    return df

def select_ion_dictionary(
    ion_types: list,
    max_length: int,
    max_charge: int,
    criteria: list,
    counts_path: str=None,
):
    passed_ions = []
    all_ions = theoretical_ions(ion_types, max_length=max_length, max_charge=max_charge, include_p=True)
    if counts_path is None:
        passed_ions = all_ions.index.tolist()
    else:
        with open(counts_path) as f:
            for line in f:
                ion, counts = line.strip().split('\t')
                counts = int(counts)
                if eval("&".join(['(%s)'%s for s in criteria])):
                    passed_ions.append(ion)
    select_ions = all_ions.loc[passed_ions]
    select_ions['index'] = np.arange(len(select_ions))
    return select_ions

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

class Scale:
    def __init__(self):
        self.ppm = lambda theor, exp: 1e6*(exp[None] - theor[:,None])/theor[:,None]
        self.diff = lambda theor, exp: exp[None] - theor[:,None]
        self.sepmod = lambda token: re.sub('\[|\]', '', token).split("UNIMOD:")
        self.mass = {
                # Amino acids
                'A': 71.037113805,'R':156.101111050,'N':114.042927470,'D':115.026943065,
                'C':103.009184505,'Q':128.058577540,'E':129.042593135,'G': 57.021463735,
                'H':137.058911875,'I':113.084064015,'L':113.084064015,'K':128.094963050,
                'M':131.040484645,'F':147.068413945,'P': 97.052763875,'S': 87.032028435,
                'T':101.047678505,'W':186.079312980,'Y':163.063328575,'V': 99.068413945,
                'O':237.147726925,'U':150.953633405,
                # Neutral losses
                'NH2':16.0187,'NH3':17.026549105,'H2O':18.010565,'CO':27.994915,
                'C2H5NOS':91.009184,'CH2SH':46.9955+0.0458,'CH3SOH':63.99828544,
                'HPO3':79.966335,'H3PO4':97.9769,'H5PO5':115.987465,'H7PO6':133.99803,
                'TMT':229.17,'RP126':154.1221,'RP127N':155.1192,'RP127C':155.1254,
                'RP128N':156.1225,'RP128C':156.1287,'RP129N':157.1258,'RP129C':157.1322,
                'RP130N':158.1291,'RP130C':158.1356,'RP131':159.1325,
                # Modifications
                'Acetyl':42.010565,'Carbamidomethyl':57.021464,'Oxidation':15.994915,'Deamidation': 0.984016,
                'Gln->pyro-Glu':-17.026549, 'Glu->pyro-Glu':-18.010565,'Phospho':79.966331,
                'Pyro-carbamidomethyl':39.994915,'CAM':57.021464,'TMT6plex':231.17747,'nAcetyl':203.079373,
                # Single atoms
                'hydrogen':1.007825035,'oxygen':15.9949146,'nitrogen':14.003074,
                'sulphur':31.9720707,'carbon':12,'phosphorus':30.973762,
                # Isotopes
                'i':1.00727646688,'iso1':1.003,'iso2':1.002,
                # Ions, immoniums, et al.
                'a':-26.9871,'b':1.007276, 'p':20.02656, 'y':19.0184,
                'ICA':76.021545, 'IDA':88.039304, 'IDB':70.028793, 'IEA':102.054954,
                'IFA':120.080775,'IFB':91.054226, 'IHA':110.071273,'IHB':82.05255,
                'IHC':121.039639,'IHD':123.055289,'IHE':138.066188,'IHF':156.076753,
                'IIA':86.096425, 'IIC':72.080776, 'IKA':101.107324,'IKB':112.07569,
                'IKC':84.080776, 'IKD':129.102239,'IKE':56.049476, 'IKF':175.118952,
                'ILA':86.096426, 'ILC':72.080776, 'IMA':104.052846,'IMB':61.010647,
                'IMC':120.047761,'INA':87.055289, 'INB':70.02874,  'IPA':70.065126, 
                'IQA':101.070939,'IQB':56.049476, 'IQC':84.04439,  'IQD':129.065854,
                'IRA':129.113473,'IRB':59.060375, 'IRC':70.065126, 'IRD':73.076025, 
                'IRE':87.091675, 'IRF':100.086924,'IRG':112.086924,'IRH':60.055624,
                'IRI':116.070605,'IRJ':175.118952,'ISA':60.04439,  'ITA':74.06004, 
                'ITB':119.0814,  'IVA':72.080776, 'IVC':55.054227, 'IVD':69.033491,
                'IWA':159.091675,'IWB':77.038577, 'IWC':117.057301,'IWD':130.065126,
                'IWE':132.080776,'IWF':170.06004, 'IWH':142.065126,'IYA':136.07569,
                'IYB':91.054227, 'IYC':107.049141,
                'ICCAM':133.04301,
                'TMTpH':230.17,'TMT126':126.1277,'TMT127N':127.1248,'TMT127C':127.1311,
                'TMT128N':128.1281,'TMT128C':128.1344,'TMT129N':129.1315,
                'TMT129C':129.1378,'TMT130N':130.1348,'TMT130C':130.1411,'TMT131':131.1382
    	}
        self.mass['1'] = self.mass['Acetyl']
        self.mass['4'] = self.mass['Carbamidomethyl']
        self.mass['7'] = self.mass['Deamidation']
        self.mass['35'] = self.mass['Oxidation']

    def calcmass(self, modseq, precursor_charge, ion, delta=0.0):
        """
        Calculating the mass of fragments

        Parameters
        ----------
        seq : Modified peptide sequence, with [UNIMOD:#] notation (str)
        precursor_charge : Precursor charge (int)
        ion : Ion type (str)

        Returns
        -------
        mass as a float

        """

        ## modification
        #Mstart = mods.find('(') if mods!='0' else 1
        #modamt = int(mods[0:Mstart])
        #modlst = []
        #if modamt>0:
        #    Mods = [re.sub("[()]",'',m).split(',') for m in 
        #             mods[Mstart:].split(')(')]
        #    for mod in Mods:
        #        [pos,aa,typ] = mod # mod position, amino acid, and type
        #        modlst.append([int(pos), self.mass[typ]])
        tokenized = [self.sepmod(m) if ":" in m else m for m in tokenize_modified_sequence(modseq)]
        seq = "".join([m[0] for m in tokenized])
        modlst = [[i, self.mass[mod[-1]]] for i, mod in enumerate(tokenized) if type(mod)==list]
        
        # isotope
        isomass = self.mass['iso1'] if '+i' in ion else self.mass['iso2']
        if ion[-1]=='i': # evaluate isotope and set variable iso
            hold = ion.split("+")
            iso = 1 if hold[-1]=='i' else int(hold[-1][:-1])
            ion = "+".join(hold[:-1]) # join back if +CO was split
        else:
            iso = 0
        
        # If internal calculate here and return
        if ion[:3]=='Int':
            ion = ion.split('+')[0] if iso!=0 else ion
            ion = ion.split('-')
            nl = self.mass[ion[1]] if len(ion)>1 else 0
            [start,extent] = [int(z) for z in ion[0][3:].split('>')]
            modmass = sum([ms[1] for ms in modlst 
                           if ((ms[0]>=start)&(ms[0]<(start+extent)))
                          ])
            return (sum([self.mass[aa] for aa in seq[start:start+extent]]) - nl
                    + iso*isomass + self.mass['i'] + modmass)
        # if TMT, calculate here and return
        if ion[:3]=='TMT':
            ion = ion.split('+')[0] if iso!=0 else ion
            return self.mass[ion] + iso*isomass
        
        # product charge
        hold = ion.split("^") # separate off the charge at the end, if at all
        charge = 1 if len(hold)==1 else int(hold[-1]) # no ^ means charge 1
        # extent
        letnum = hold[0].split('-')[0];let = letnum[0] # ion type and extent is always first string separated by -
        num = int(letnum[1:]) if ((let!='p')&(let!='I')) else 0 # p type ions never have number
        
        # neutral loss
        nl=0
        hold = hold[0].split('-')[1:] # most are minus, separated by -
        """If NH2-CO-CH2SH, make the switch to C2H5NOS. Get rid of CO and CH2SH.""" 
        #if len(hold)>0 and ('NH2' in hold[0]):
        #    mult = (int(hold[0][0]) 
        #            if ((ord(hold[0][0])>=48) & (ord(hold[0][0])<=57)) else '')
        #    hold[0] = str(mult)+'C2H5NOS';del(hold[1],hold[1])
        for item in hold:
            if '+' in item: # only CO can be a +
                items = item.split('+') # split it e.g. H2O+CO
                nl-=self.mass[items[0]] # always minus the first
                nl+=self.mass[items[1]] # always plus the second
            else:
                if (ord(item[0])>=48) & (ord(item[0])<=57): # if there are e.g. 2 waters -> 2H2O 
                    mult = int(item[0])
                    item = item[1:]
                else:
                    mult = 1
                nl-=mult*self.mass[item]

        if let == 'A':
            sm = sum([self.mass[aa] for aa in seq[:num]])
            modmass = sum([mod[1] for mod in modlst if num>mod[0]]) # if modification is before extent from n terminus
            return self.mass['i'] + (sm + modmass - self.mass['CO'] + self.mass['hydrogen'] + nl + iso*isomass + delta) / charge
        elif let == 'a':
            sm = sum([self.mass[aa] for aa in seq[:num]])
            modmass = sum([mod[1] for mod in modlst if num>mod[0]]) # if modification is before extent from n terminus
            return self.mass['i'] + (sm + modmass - self.mass['CO'] + nl + iso*isomass + delta) / charge
        elif let=='b':
            sm = sum([self.mass[aa] for aa in seq[:num]])
            modmass = sum([mod[1] for mod in modlst if num>mod[0]]) # if modification is before extent from n terminus
            return self.mass['i'] + (sm + modmass + nl + iso*isomass + delta) / charge
        elif let=='C': # c-prime (PRIMARY PRODUCT)
            sm = sum([self.mass[aa] for aa in seq[:num]])
            modmass = sum([mod[1] for mod in modlst if num>mod[0]])
            return self.mass['i'] + (sm + self.mass['NH3'] + modmass + nl + iso*isomass + delta) / charge
        elif let=='c': # c-dot
            sm = sum([self.mass[aa] for aa in seq[:num]])
            modmass = sum([mod[1] for mod in modlst if num>mod[0]])
            return self.mass['i'] + (sm + self.mass['NH2'] + modmass + nl + iso*isomass + delta) / charge 
        elif let=='p':
            # p eq.: proton + (aa+H2O+mods+i)/charge
            sm = sum([self.mass[aa] for aa in seq])
            charge = int(precursor_charge) if ((charge==1)&('^1' not in ion)) else charge
            modmass = sum([mod[1] for mod in modlst]) # add all modifications
            return self.mass['i'] + (sm + self.mass['H2O'] + modmass + nl + iso*isomass + delta) / charge
        elif let=='v':
            sm = sum([self.mass[aa] for aa in seq[-(num-1):]])
            modmass = sum([mod[1] for mod in modlst if ((len(seq)-num)<=mod[0])])
            return self.mass['i'] + (sm + 2*(self.mass['carbon']+self.mass['oxygen']+self.mass['hydrogen']) + modmass + nl + iso*isomass + delta) / charge       
        elif let=='X':
            sm = sum([self.mass[aa] for aa in seq[-num:]])
            modmass = sum([mod[1] for mod in modlst if ((len(seq)-num)<=mod[0])])
            return self.mass['i'] + (sm + 2*self.mass['oxygen'] + self.mass['carbon'] + self.mass['hydrogen'] + modmass + nl + iso*isomass + delta) / charge
        elif let=='x':
            sm = sum([self.mass[aa] for aa in seq[-num:]])
            modmass = sum([mod[1] for mod in modlst if ((len(seq)-num)<=mod[0])])
            return self.mass['i'] + (sm + 2*self.mass['oxygen'] + self.mass['carbon'] + modmass + nl + iso*isomass + delta) / charge
        elif let=='y':
            sm = sum([self.mass[aa] for aa in seq[-num:]])
            modmass = sum([mod[1] for mod in modlst if ((len(seq)-num)<=mod[0])]) # if modification is before extent from c terminus
            return self.mass['i'] + (sm + self.mass['H2O'] + modmass + nl + iso*isomass + delta) / charge #0.9936 1.002
        elif let=='Z': # z-prime
            sm = sum([self.mass[aa] for aa in seq[-num:]])
            modmass = sum([mod[1] for mod in modlst if ((len(seq)-num)<=mod[0])])
            return self.mass['i'] + (sm + self.mass['oxygen'] - self.mass['nitrogen'] + self.mass['hydrogen'] + modmass + nl + iso*isomass + delta) / charge
        elif let=='z': # z_r in spectrum fundamentals, z-dot (PRIMARY PRODUCT)
            sm = sum([self.mass[aa] for aa in seq[-num:]])
            modmass = sum([mod[1] for mod in modlst if ((len(seq)-num)<=mod[0])])
            return self.mass['i'] + (sm + self.mass['oxygen'] - self.mass['nitrogen'] + modmass + nl + iso*isomass + delta) / charge # + 0.0027545
        elif let=='I':
            sm = self.mass[letnum]
            return (sm + iso*isomass + delta) / charge
        else:
            return False

    def LowMz(self, mz1, mz2, maxmz=300, thr=3e-3, thr2=4e-3, ions=[]):
        """
        Special criterion for very low mz, using absolute difference

        Parameters
        ----------
        mz1 : mzs for predicted peaks
        mz2 : mzs for experimental peaks
        maxmz : Maximum m/z, below which criterion will be applied. 
                The default is 300.
        thr : Threshold for maximum absolute value difference. The default 
              is 3e-3.
        thr2 : Threshold for monoisotopic peaks. The default is 4e-3.

        Returns
        -------
        Indices of true positives to remove

        """
        
        lt = mz1<maxmz
        eye = np.array([True if 'i' in ion else False for ion in ions])
        diff = abs(self.diff(mz1, mz2))
        S = np.where(lt&(eye==False))[0]
        S2 = np.where(lt&eye)[0]
        tp1 = np.sort(np.append(
            S[np.where(diff[S].min(1)<thr)[0]], S2[np.where(diff[S2].min(1)<thr2)[0]]
        ))
        tp2 = diff[tp1].argmin(1)
        
        return (tp1,tp2)

    def match(self, mz1, mz2, thr=10, spl=[0], 
              lowmz=False, pions=None, typ='ppm'):
        """
        Find matches between 2 peaks lists' m/z values
        - All non-matches will be classifiedfalse positives for mz1, and false 
           negatives for mz2

        Parameters
        ----------
        mz1 : First vector of peak m/z's', preferably the predicted peaks
        mz2 : Second vector of peak m/z's, preferably the experimental peaks'
        thr : Thresholds, under which peaks are considered matched. There
               should be 1 more value than length of 'spl'. The default is 
               [15,20,25].
        spl : m/z values on which to split dataset for different thresholds.
               Splits will be under first value, between middle values, and 
               over the top value. The default is '[800,1200]'.
        lowmz : Boolean whether to apply LowMz criterion or not after matching.
                 Note that between the thresholds and LowMz, this is an OR
                 criterion, so only 1 needs to be satisfied. Default False.
        pions : Array of predicted ion types. This is necessary for LowMz=True.
        typ : Distance metric between peaks to be matched. Use either absolute 
               difference ('abs'), or ppm. The default is 'ppm'.

        Returns
        -------
        TP : Tuple of true positive indices for [0] mz1 peaks and [1] mz2 peaks
        FP1 : False positives indices for mz1 peaks
        FN2 : False negatives indices for mz2 peaks

        """
        if type(thr) in [float, int]:
            thr = [thr]
        
        delta = self.ppm(mz1,mz2) if typ=='ppm' else self.diff(mz1,mz2)
        TP1=[];TP2=[];FP1=[];FN2=[]
        for i,s in enumerate(thr):
            # smaller than first, between middle, over the top
            split = (
                mz1>spl[i-1] if i==(len(thr)-1) else 
                (mz1<spl[i] if i==0 else 
                ( (mz1>spl[i-1]) & (mz1<spl[i]) ) 
                ) 
            )
            S = np.where(split)[0]
            tp1 = S[np.where(abs(delta[split]).min(1) < s)[0]]
            tp2 = abs(delta[tp1]).argmin(1)
            # assert(sum((diff<thr).sum(1)>1)==0)
            fp1 = S[np.where(abs(delta[split]).min(1) > s)[0]]
            split = (mz2>spl[i-1] if i==(len(thr)-1) else 
                     (mz2<spl[i] if i==0 else ((mz2>spl[i-1])&(mz2<spl[i])) ) )
            fn2 = np.where(split)[0][np.where(abs(delta[:,split]).min(0) > s)[0]]
            TP1.append(tp1);TP2.append(tp2);FP1.append(fp1);FN2.append(fn2)
        TP1 = np.concatenate(TP1);TP2 = np.concatenate(TP2)
        FP1 = np.concatenate(FP1);FN2 = np.concatenate(FN2)
        
        # Put any rejected matches, that satisfy the LowMz criterion, into the
        # True positive lists, and remove from the False lists.
        if lowmz:
            # Get indices of matches satisfying LowMz criterion
            ltp = self.LowMz(mz1, mz2, ions=pions)
            # Iterate through that list
            for tp1,tp2 in zip(*ltp):
                # See if they were rejected by previous matching thresholds
                if tp1 in FP1:
                    FP1 = np.delete(FP1, np.where(FP1==tp1))
                    TP1 = np.append(TP1, tp1)
                    # exp list can have duplicates (before tiebreaker)
                    TP2 = np.append(TP2, tp2)
                if tp2 in FN2:
                    FN2 = np.delete(FN2, np.where(FN2==tp2))
            TP1 = np.sort(TP1)
            TP2 = np.sort(TP2)
        return (TP1,TP2), FP1, FN2

def tiebreak(TP, theoretical_df, real_mzs, thr=2e-2):
    """
    Must get rid of multiple predicted peaks (peaks1) that match a
     single experimental peak (peaks2).

	Parameters
	----------
	theoretical_df : dataframe with columns [[ion, mz]]
	real_mzs : full raw mz vector
	TP : Original true positive indices (TP1, TP2)
    thr : Tiebreaker threshold (daltons)

	Returns
	-------
	peaks1 : peaks1 without the updated mz and ab values, and without
			  the duplicate peaks.
	TP1 : New, smaller, true positive indices for predicted peaks
	TP2 : New, smaller, true positive indices for experimental peaks

	"""
    TP1, TP2 = TP
    theor_mzs_ = theoretical_df['mz'].to_numpy()[TP1]
    theor_dict_ = theoretical_df['ion'].to_numpy()[TP1]
    real_mzs_  = real_mzs[TP2]

    # unique experimental TP indices, and the counts for each
    uniqs, i, cnts = np.unique(TP2, return_index=True, return_counts=True)
    if len(TP2)==len(uniqs): 
        return TP1, TP2
    multi_inds = uniqs[np.where(cnts>1)[0]] # all duplicate TP2 values
    dels=[]
    for ind in multi_inds:
        # indices of TP2 array that match its duplicate value (which is in turn an index for the peaks2 array)
        nest_inds = np.where(TP2==ind)[0]
        expmz = real_mzs_[nest_inds]
        predmz = theor_mzs_[nest_inds]
        
        # Make sure duplicates are below tiebreaker threshold
        #inds2 = np.where(abs(predmz-expmz)<thr)[0]
        #if len(inds2)==0:
        #    for m in nest_inds[np.argsort(abs(predmz-expmz))[1:]]: dels.append(m)
        #else:
        
        I = TP1[nest_inds] # Values of TP1 for duplicates
        # use m/z of the closer duplicate predicted peak
        winner = abs(predmz-expmz).argmin() # index of nest_inds
        closer_tp1_index = nest_inds[winner] # index of TP1 vector
        closer_mz_index = I[winner] # index of thoer_mzs
        # collect indices of others to get rid of
        for m in nest_inds:
            if m != closer_tp1_index: dels.append(m)
			
	#global_dels = TP1[dels]
    TP1 = np.delete(TP1, dels)
    TP2 = np.delete(TP2, dels)
    assert(len(TP2)==len(uniqs))
    return TP1, TP2#, global_dels

def my_annotation_function(psms, theor_dict, threshold_ppm=20, p_window=0):
    ion_counts = {}
    scale = Scale()
    df = {
        'mz': [], 'int': [],  
        'matched_inds': [], 'matched_ions': [], 'matched_ppm': [],
    }
    for seq, modseq, charge, raw_mz, raw_int, mass in zip(
        psms['SEQUENCE'],
        psms['MODIFIED_SEQUENCE'], 
        psms['PRECURSOR_CHARGE'],
        psms['MZ'],
        psms['INTENSITIES'],
        psms['MASS'],
    ):
        # Calculate theoretical spectrum
        peptide_length = len(seq)
        ions = theor_dict.query(f"length < {peptide_length} and charge <= {charge}")
        theoretical_mz = np.array([scale.calcmass(modseq, charge, ion) for ion in ions.index])
        
        # Exclude all peaks with m/z between [p-p_window, p+p_window]
        # - Create a mask to filter variables 'theoretical_mz' and 'ions'
        if p_window > 0:
            p = scale.calcmass(modseq, charge, 'p')
            exclusion_mask = (theoretical_mz > (p + p_window)) | (theoretical_mz < (p - p_window))
        else:
            exclusion_mask = np.array(len(theoretical_mz)*[True])
        theoretical_mz = theoretical_mz[exclusion_mask]
        theor_df = pd.DataFrame({'ion': ions.index.to_numpy()[exclusion_mask], 'mz': theoretical_mz})
        
        # Match peaks
        TP, FP, FN = scale.match(theoretical_mz, raw_mz, thr=threshold_ppm)
        TP = tiebreak(TP, theor_df, raw_mz)
        raw_mz_ = raw_mz[TP[1]]
        raw_int_ = raw_int[TP[1]]
        theoretical_mz_ = theoretical_mz[TP[0]]
        ppm = 1e6 * (raw_mz_ - theoretical_mz_) / theoretical_mz_
        matched_ions = ions.index.to_numpy()[exclusion_mask][TP[0]]
        
        # Collect results
        #df['ions'].append(matched_ions)
        df['mz'].append(raw_mz_)
        df['int'].append(raw_int_)
        df['matched_inds'].append(TP[1])
        df['matched_ions'].append(matched_ions)
        df['matched_ppm'].append(ppm)
        
        for I in matched_ions:
            if I not in ion_counts:
                ion_counts[I] = 0
            ion_counts[I] += 1

    return {
        'dataframe': pd.DataFrame(df),
        'statistics': ion_counts,
    }

# if __name__ == '__main__':
#     scale = Scale()
    
#     import pandas as pd
#     df = pd.read_parquet("/cmnfs/data/proteomics/shabaz_exotic/processed/merged_search/ECD/all_psms_Trypsin.parquet")
    
#     ions = theoretical_ions(['C','Z','z'])

#     for sequence, charge in zip(df['MODIFIED_SEQUENCE'], df['PRECURSOR_CHARGE']):
#         theoretical_masses = [scale.calcmass(sequence, charge, ion) for ion in ions]

#     (TP1, TP2), FP, FN = scale.match(mz_array, theoretical_masses)
        

