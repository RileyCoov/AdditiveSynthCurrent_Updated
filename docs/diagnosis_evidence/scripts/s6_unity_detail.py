import sys; sys.path.insert(0,'an')
from load import *
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.signal import stft, butter, sosfiltfilt
LSB=1/32768
def runs(mask, minlen=64):
    d=np.diff(np.concatenate([[0],mask.astype(int),[0]])); s=np.where(d==1)[0]; e=np.where(d==-1)[0]
    keep=(e-s)>=minlen; return s[keep], e[keep]
def analyze(n, plot=False):
    xo,fs=load(n,'original'); xu,_=load(n,'unity'); e=xu-xo
    ex=np.abs(e)<=1.01*LSB; s,t=runs(ex,64)
    pt=np.zeros(len(xo),bool)
    for a,b in zip(s,t): pt[a:b]=True
    # dilate by 1024 each side to exclude crossfades
    from scipy.ndimage import binary_dilation
    ptd=binary_dilation(pt, np.ones(2*2048+1,bool))
    k=512; pad=np.pad(xo**2,(k//2,k//2),mode='edge'); loc=np.sqrt(np.convolve(pad,np.ones(k)/k,'valid')[:len(xo)])
    act=loc>10**(-50/20)
    m=act&~ptd
    ser_model=10*np.log10(np.sum(xo[m]**2)/np.sum(e[m]**2)) if m.sum() else np.nan
    # crossfade zone only
    xf=ptd&~pt&act
    ser_xf=10*np.log10(np.sum(xo[xf]**2)/np.sum(e[xf]**2)) if xf.sum() else np.nan
    # band SER in model-only region via STFT masks
    f,tt,X=stft(xo,fs,'hann',nperseg=2048,noverlap=1536); _,_,E=stft(e,fs,'hann',nperseg=2048,noverlap=1536)
    fr_idx=np.clip((tt*fs).astype(int),0,len(xo)-1); mf=m[fr_idx]
    bands=[(20,60),(60,100),(100,200),(200,500),(500,2000),(2000,5000),(5000,12000),(12000,20000)]
    bs=[]
    for lo,hi in bands:
        sel=(f>=lo)&(f<hi)
        bs.append(10*np.log10(np.sum(np.abs(X[sel][:,mf])**2)/np.sum(np.abs(E[sel][:,mf])**2)))
    print(f"{n:30s} model-only frac {m.sum()/act.sum():.2f}  SER_model {ser_model:6.2f}  SER_xfade {ser_xf:6.2f} | band SER(model) "+" ".join(f"{lo}-{hi}:{b:5.1f}" for (lo,hi),b in zip(bands,bs)))
    return xo,xu,e,pt,ptd,fs
for n in NAMES:
    analyze(n)
