from common import *
from scipy.signal import hilbert
fs=48000
for n,v,f0 in [('440sawtooth','down_fifth',319.34),('440sawtooth','up_fifth',717.77),('300hzSaw_5sec','down_fifth',199.95)]:
    x,_=load(n,v)
    y=bp(x,fs,f0-60,f0+60,order=4); a=np.abs(hilbert(y))
    blk=1200
    g=[];pk=[]
    for i in range(blk*2,len(x)-blk*2,blk):
        g.append(a[i:i+blk].mean()); pk.append(np.abs(x[i:i+blk]).max())
    g=np.array(g); g=20*np.log10(g/np.median(g)); pk=np.array(pk)
    atceil=pk>=0.9504
    print(f"{n} {v}: blocks {len(g)}; blocks with peak at ceiling: {atceil.sum()}; mean gain(dB) when at ceiling {g[atceil].mean() if atceil.any() else np.nan:+.2f} vs not {g[~atceil].mean():+.2f}; corr(gain, peak) {np.corrcoef(g,pk)[0,1]:+.2f}")
    print("   peak per 25ms (first 40):",np.round(pk[:40],3).tolist())
