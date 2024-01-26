import numpy as np

ppm = lambda theor, exp: 1e6*(exp[None] - theor[:,None])/theor[:,None]

def match(
    mz1, 
    mz2, 
    thr=[15,20,25], 
    spl=[800,1200], 
    lowmz=False, 
    pions=None
):
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
    
    delta = ppm(mz1,mz2)
    TP1=[];TP2=[];FP1=[];FN2=[]
    for i,s in enumerate(thr):
        # smaller than first, between middle, over the top
        split = (mz1>spl[i-1] if i==(len(thr)-1) else 
                 (mz1<spl[i] if i==0 else ((mz1>spl[i-1])&(mz1<spl[i])) ) )
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
    """if lowmz:
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
        TP2 = np.sort(TP2)"""
    return (TP1,TP2), FP1, FN2