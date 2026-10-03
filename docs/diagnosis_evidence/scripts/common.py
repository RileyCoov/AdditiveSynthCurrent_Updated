import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # sibling modules (load.py etc.)
from load import *
import numpy as np
from scipy.signal import butter, sosfiltfilt, stft
R={'original':1.0,'unity':1.0,'up_fifth':2**(7/12),'down_fifth':2**(-7/12)}
def bp(x,fs,lo,hi,order=4):
    hi=min(hi,fs/2*0.98)
    sos=butter(order,[lo/(fs/2),hi/(fs/2)],'band',output='sos'); return sosfiltfilt(sos,x)
def rms_env(x,win):
    k=np.ones(win)/win; return np.sqrt(np.convolve(x**2,k,'same')+1e-20)
def onsets(x,fs,lo,hi,nfft=1024,hop=128,thr_k=1.5,min_gap_s=0.09):
    f,t,X=stft(x,fs,'hann',nperseg=nfft,noverlap=nfft-hop)
    sel=(f>=lo)&(f<hi); M=np.log1p(1000*np.abs(X[sel]))
    flux=np.maximum(np.diff(M,axis=1),0).sum(0); flux=np.concatenate([[0],flux])
    # adaptive threshold
    from scipy.ndimage import median_filter
    med=median_filter(flux,size=int(0.5*fs/hop)); 
    th=med+thr_k*np.std(flux)
    pk=[i for i in range(1,len(flux)-1) if flux[i]>th[i] and flux[i]>=flux[i-1] and flux[i]>=flux[i+1]]
    out=[]
    for i in pk:
        ts=t[i]
        if not out or ts-out[-1]>min_gap_s: out.append(ts)
        elif flux[i]>flux[np.searchsorted(t,out[-1])]: out[-1]=ts
    # STFT frame centre time -> onset approx at centre; refine with band envelope rise
    return np.array(out), (t,flux)
