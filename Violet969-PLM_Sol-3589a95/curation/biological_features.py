"""Label-blind sequence annotation primitives, protocol version 1.0.

Coordinates returned by regions() are 1-based inclusive. Unknown amino acids
are never assigned zero propensity. See BIOLOGICAL_ANNOTATION_METHODS.md.
"""
from collections import Counter
import numpy as np
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint

AA = 'ACDEFGHIKLMNPQRSTVWY'
KD = dict(zip(AA, [1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]))
# Conchillo-Sole et al. 2007, Additional file 1 (not fitted to UESolDS).
A3V = dict(zip('IFVLYWMCATSPGKHQRNED', [1.822,1.754,1.594,1.380,1.159,1.037,.910,.604,-.036,-.159,-.294,-.334,-.535,-.931,-1.033,-1.231,-1.240,-1.302,-1.412,-1.836]))
GROUPS = dict(positive='KRH', negative='DE', charged='DEKRH',
              polar_uncharged='STNQCY', nonpolar='AGILMFPVW', aromatic='FWY',
              hydrophobic='AVILMFWY', small='AGST', proline='P', cysteine='C')
HST = -.02
WINDOWS = (9, 19)

def canonical(s):
    return bool(s) and not (set(s) - set(AA))

def regions(mask, minimum=1):
    mask = np.asarray(mask, dtype=bool)
    edges = np.diff(np.r_[False, mask, False].astype(int))
    return [(int(a+1), int(b)) for a,b in zip(np.where(edges==1)[0], np.where(edges==-1)[0]) if b-a>=minimum]

def clip_regions(intervals, length):
    return [(a, min(b,length)) for a,b in intervals if a<=length]

def window_values(values, width):
    """Only complete windows; NaN propagates to any overlapping window."""
    v=np.asarray(values, dtype=float)
    return np.convolve(v, np.ones(width)/width, mode='valid') if len(v)>=width else np.array([])

def sidechain_charge(s, ph):
    pos={'K':10., 'R':12., 'H':5.98}
    neg={'D':4.05, 'E':4.45, 'C':9., 'Y':10.}
    return np.array([1/(1+10**(ph-pos[a])) if a in pos else
                     -1/(1+10**(neg[a]-ph)) if a in neg else
                     0. if a in AA else np.nan for a in s])

def charge_profile(s, ph):
    """Fractional charge assigned to residues, including actual free termini."""
    q=sidechain_charge(s, ph)
    if not canonical(s): return q*np.nan
    ip=IsoelectricPoint(s)
    q[0]+=1/(1+10**(ph-ip.pos_pKs['Nterm']))
    q[-1]-=1/(1+10**(ip.neg_pKs['Cterm']-ph))
    return q

def biochemistry(s):
    n=len(s); c=Counter(s)
    result={'length':n,'nonstandard_count':sum(v for a,v in c.items() if a not in AA),
            'nonstandard_residues':''.join(sorted(set(s)-set(AA)))}
    # Canonical fractions use original length, not number of recognized residues.
    result.update({'frac_'+a:c[a]/n if n else np.nan for a in AA})
    result.update({'frac_'+k:sum(c[a] for a in v)/n if n else np.nan for k,v in GROUPS.items()})
    result['entropy_bits']=-sum((v/n)*np.log2(v/n) for v in c.values()) if canonical(s) else np.nan
    result['charge_proxy']=(c['K']+c['R']-c['D']-c['E'])/n if canonical(s) else np.nan
    result['hydropathy_mean']=np.mean([KD[a] for a in s]) if canonical(s) else np.nan
    if canonical(s):
        ip=IsoelectricPoint(s)
        # Wider bracket than Biopython default; otherwise acidic pI clips at 4.05.
        result['pi_estimate']=ip.pi(pH=7., min_=0., max_=14.)
        for ph in (7.,7.4):
            key=str(ph).replace('.','p')
            result['net_charge_ph'+key]=ip.charge_at_pH(ph)
            result['charge_density_ph'+key]=ip.charge_at_pH(ph)/n
    else:
        result['pi_estimate']=np.nan
        for key in ('7p0','7p4'):
            result['net_charge_ph'+key]=result['charge_density_ph'+key]=np.nan
    for name, values in [('hydropathy', [KD.get(a,np.nan) for a in s]),
                         ('sidechain_charge_ph7p0',sidechain_charge(s,7.)),
                         ('sidechain_charge_ph7p4',sidechain_charge(s,7.4)),
                         ('positive_fraction',[float(a in 'KRH') if a in AA else np.nan for a in s]),
                         ('negative_fraction',[float(a in 'DE') if a in AA else np.nan for a in s])]:
        for width in WINDOWS:
            v=window_values(values,width)
            for op in ('min','max'):
                key=f'{name}_w{width}_{op}'
                if len(v) and np.isfinite(v).any():
                    i=int(np.nanargmin(v) if op=='min' else np.nanargmax(v))
                    result[key]=v[i];result[key+'_start']=i+1;result[key+'_end']=i+width
                else:
                    result[key]=result[key+'_start']=result[key+'_end']=np.nan
    return result

def aggregation(s, terminal_policy='paper'):
    """Published AGGRESCAN sequence algorithm, independently reimplemented.

    This is not the authors' C executable. Virtual terminal residues and copied
    edge-window values follow the 2007 paper. Does not compute area descriptors.
    """
    if not canonical(s) or len(s)<5:
        return np.full(len(s),np.nan), [], 'unsupported_input'
    w=5 if len(s)<=75 else 7 if len(s)<=175 else 9 if len(s)<=300 else 11
    if terminal_policy not in ('paper','published_example'):raise ValueError(terminal_policy)
    nt,ct=(A3V['K']+A3V['R'])/2,(A3V['D']+A3V['E'])/2
    # Supplementary example numerically uses the opposite terminal order.
    # Keep the prose-defined rule primary; expose discrepancy as sensitivity.
    if terminal_policy=='published_example':nt,ct=ct,nt
    v=np.array([nt, *[A3V[a] for a in s], ct])
    # Scale values are specified to 0.001; virtual averages require 0.0005.
    # Integer summation avoids classifying an exactly-threshold window positive
    # merely because floating-point convolution perturbs -0.02.
    sums=np.convolve(np.rint(v*2000).astype(np.int64),np.ones(w,dtype=np.int64),mode='valid')
    profile=np.pad(sums/(2000*w), (w//2-1,w//2-1), mode='edge')
    assert len(profile)==len(s)
    intervals=regions((profile>HST)&np.array([a!='P' for a in s]),5)
    return profile,intervals,'ok' if intervals else 'zero_features'

def regional_summary(scores, intervals):
    v=np.asarray(scores)
    if not len(v) or not np.isfinite(v).all():
        return dict(mean=np.nan,maximum=np.nan,n_regions=np.nan,fraction=np.nan,longest=np.nan)
    return dict(mean=float(v.mean()), maximum=float(v.max()), n_regions=len(intervals),
                fraction=sum(b-a+1 for a,b in intervals)/len(v), longest=max([b-a+1 for a,b in intervals],default=0))

def membrane_summary(segments,length):
    """DeepTMHMM named regions: crossing segments versus membrane association."""
    crossing=clip_regions([(a,b) for name,a,b in segments if name in ('TMhelix','Beta sheet')],length)
    associated=clip_regions([(a,b) for name,a,b in segments if name in ('TMhelix','Beta sheet','reentrant','interfacial')],length)
    return dict(tm_segment_count=len(crossing),tm_fraction=sum(b-a+1 for a,b in crossing)/length,
                membrane_associated_fraction=sum(b-a+1 for a,b in associated)/length)
