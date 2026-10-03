from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
n=sys.argv[1] if len(sys.argv)>1 else '1985'
lo,hi=float(sys.argv[2]) if len(sys.argv)>2 else 40, float(sys.argv[3]) if len(sys.argv)>3 else 160
fs=48000; X={v:load(n,v)[0] for v in VERS}
win=int(0.010*fs)
def bandpow(v):
    r=R[v]; y=bp(X[v],fs,lo*r,hi*r,order=6)
    return np.convolve(y**2,np.ones(win)/win,'same')
P={v:bandpow(v) for v in VERS}
ons=np.load(f'an/{n}_onsets.npy')
D={v:10*np.log10((P[v]+1e-12)/(P['original']+1e-12)) for v in VERS[1:]}
t=np.arange(len(X['original']))/fs
fig,ax=plt.subplots(3,1,figsize=(18,9),sharex=True)
for i,v in enumerate(VERS[1:]):
    ax[i].plot(t[::48],D[v][::48],lw=0.6); ax[i].set_ylim(-25,15); ax[i].axhline(0,c='k',lw=.5); ax[i].axhline(-6,c='r',lw=.5,ls='--')
    for o in ons: ax[i].axvline(o,c='g',lw=0.4)
    for o in [1.350,2.977,4.174,4.739]: ax[i].axvline(o,c='m',lw=1.0,ls=':')
    ax[i].set_title(f'{n} {v}: band power ratio vs original, src band {lo}-{hi} Hz (mapped), 10 ms; green=onsets, magenta=README times')
plt.tight_layout(); plt.savefig(f'an/{n}_lfratio_{int(lo)}_{int(hi)}.png',dpi=60)
# list intervals below -6 dB lasting >=10 ms where original band power is within 25 dB of its 95th pct (i.e., bass present)
act=P['original']>np.percentile(P['original'],95)*10**(-25/10)
for v in VERS[1:]:
    m=(D[v]<-6)&act
    d=np.diff(np.concatenate([[0],m.astype(int),[0]])); s=np.where(d==1)[0]; e=np.where(d==-1)[0]
    keep=(e-s)>=int(0.010*fs)
    s,e=s[keep],e[keep]
    near=[np.min(np.abs(ons-(a/fs))) for a in s]
    print(f"{v}: {len(s)} intervals <-6 dB >=10 ms (bass active); total {np.sum(e-s)/fs*1000:.0f} ms; frac within 40 ms of onset: {np.mean(np.array(near)<0.04) if len(s) else 0:.2f}")
    if v!='unity':
        for a,b,nr in zip(s,e,near):
            print(f"   {a/fs:6.3f}-{b/fs:6.3f}s ({b-a:5d} smp, {1000*(b-a)/fs:4.1f} ms) min {D[v][a:b].min():6.1f} dB  nearest onset {nr*1000:5.1f} ms")
