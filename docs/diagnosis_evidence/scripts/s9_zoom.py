from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
n='1985'
X={v:load(n,v)[0] for v in VERS}; fs=48000
times=[float(a) for a in sys.argv[1].split(',')] if len(sys.argv)>1 else [1.350,2.977,4.174,4.739]
fig,ax=plt.subplots(len(times),1,figsize=(14,3.2*len(times)))
for i,t0 in enumerate(times):
    a=int((t0-0.12)*fs); b=int((t0+0.12)*fs); tt=np.arange(a,b)/fs
    for v,c in zip(VERS,['k','tab:orange','tab:green','tab:red']):
        r=R[v]; y=bp(X[v],fs,25*r,180*r)
        ax[i].plot(tt,y[a:b]+{'original':0.3,'unity':0.1,'up_fifth':-0.1,'down_fifth':-0.3}[v],c=c,lw=0.8,label=v)
    for k in range(a//1024,b//1024+1): ax[i].axvline(k*1024/fs,c='gray',lw=0.3)
    ax[i].axvline(t0,c='b',lw=1); ax[i].set_title(f'LF (25-180 Hz src, mapped) around {t0}s; grid=1024-hop'); ax[i].legend(fontsize=6,loc='upper right')
plt.tight_layout(); plt.savefig('an/1985_zoomLF.png',dpi=65)
