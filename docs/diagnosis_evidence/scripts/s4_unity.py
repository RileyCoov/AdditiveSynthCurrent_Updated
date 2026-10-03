import sys; sys.path.insert(0,'an')
from load import *
import numpy as np
from scipy.signal import butter, sosfiltfilt, stft
LSB=1/32768
bands=[(20,100),(100,200),(200,500),(500,1000),(1000,2500),(2500,5000),(5000,10000),(10000,20000)]
def bp(x,fs,lo,hi):
    sos=butter(4,[lo/(fs/2),min(hi,fs/2*0.99)/(fs/2)],'band',output='sos'); return sosfiltfilt(sos,x)
print("UNITY vs ORIGINAL")
hdr=f"{'name':30s} {'SER_dB':>7s} {'exact%':>7s} {'exactActive%':>12s} "+" ".join(f"{lo}-{hi}"[:9].rjust(9) for lo,hi in bands)
print(hdr)
res={}
for n in NAMES:
    xo,fs=load(n,'original'); xu,_=load(n,'unity'); e=xu-xo
    ser=10*np.log10(np.sum(xo**2)/np.sum(e**2))
    exact=np.abs(e)<=1.01*LSB
    # active samples: local rms > -50 dBFS
    k=512; pad=np.pad(xo**2,(k//2,k//2),mode='edge'); loc=np.sqrt(np.convolve(pad,np.ones(k)/k,'valid')[:len(xo)])
    act=loc>10**(-50/20)
    bser=[]
    for lo,hi in bands:
        xb=bp(xo,fs,lo,hi); eb=bp(e,fs,lo,hi)
        bser.append(10*np.log10((np.sum(xb**2)+1e-20)/(np.sum(eb**2)+1e-20)))
    res[n]=(ser,exact.mean()*100,(exact&act).sum()/act.sum()*100,bser)
    print(f"{n:30s} {ser:7.2f} {exact.mean()*100:7.1f} {(exact&act).sum()/act.sum()*100:12.1f} "+" ".join(f"{b:9.1f}" for b in bser))
np.save('an/unity_res.npy',res,allow_pickle=True)
