import sys; sys.path.insert(0,'an')
from load import *
import numpy as np
LSB=1/32768
def runs(mask, minlen=64):
    d=np.diff(np.concatenate([[0],mask.astype(int),[0]])); s=np.where(d==1)[0]; e=np.where(d==-1)[0]
    keep=(e-s)>=minlen; return s[keep], e[keep]
for n in ['1985','DrumLoopShort','HappyMono','take-me-out','PianoSampleMono','FairlightAdditiveWaveTA1_C3','300hzSine_5sec']:
    xo,fs=load(n,'original'); xu,_=load(n,'unity'); e=xu-xo
    ex=np.abs(e)<=1.01*LSB
    s,t=runs(ex,64)
    L=t-s
    print(f"\n{n}: {len(s)} exact runs >=64 smp; len median {np.median(L) if len(L) else 0:.0f}, min {L.min() if len(L) else 0}, max {L.max() if len(L) else 0}")
    print("  start mod 1024:", np.bincount((s%1024)//128, minlength=8), " start mod 256:", np.bincount((s%256)//32,minlength=8))
    print("  len hist (<=512,<=1024,<=2048,<=4096,>4096):", [(L<=512).sum(),((L>512)&(L<=1024)).sum(),((L>1024)&(L<=2048)).sum(),((L>2048)&(L<=4096)).sum(),(L>4096).sum()])
    if n=='1985':
        for a,b in zip(s,t):
            print(f"   {a:7d}-{b:7d} ({a/fs:6.3f}-{b/fs:6.3f}s, {b-a:5d} smp, start%1024={a%1024}, start%256={a%256})")
