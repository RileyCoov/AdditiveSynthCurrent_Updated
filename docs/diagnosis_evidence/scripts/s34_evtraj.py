from tf import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fs=48000; n='1985'
cases=[(1.0,1.6,530,'530'),(1.3,1.9,1205,'1205'),(4.2,4.9,3441,'3441'),(1.85,2.35,2777,'2777'),(5.3,5.8,3500,'3500')]
X={v:load(n,v)[0] for v in VERS}
S={}
for v in VERS:
    f,t,Y=spec(X[v],fs); S[v]=(f,t,20*np.log10(Y+1e-9))
fig,ax=plt.subplots(len(cases),2,figsize=(16,3.4*len(cases)))
for i,(t0,t1,fsrc,lab) in enumerate(cases):
    for v,c in zip(VERS,['k','tab:orange','tab:green','tab:red']):
        f,t,L=S[v]; r=R[v]
        sel=(t>=t0)&(t<=t1)
        # level at mapped freq: max within ±1.5% of r*fsrc
        band=(f>=r*fsrc*0.985)&(f<=r*fsrc*1.015)
        ax[i,0].plot(t[sel],L[band][:,sel].max(0),c=c,label=v,lw=.9)
        # local spectrum slice at time of max output excess
    ax[i,0].set_title(f'level at src {fsrc} Hz (±1.5%), mapped; {t0}-{t1}s'); ax[i,0].legend(fontsize=6); ax[i,0].grid(alpha=.3)
    # spectral slice at event centre in source coords
    tc=0.5*(t0+t1)
    for v,c in zip(VERS,['k','tab:orange','tab:green','tab:red']):
        f,t,L=S[v]; r=R[v]; j=np.argmin(np.abs(t-tc))
        # find time of max at band in up
        ax[i,1].plot(f/r,L[:,j],c=c,lw=.7,label=v)
    ax[i,1].set_xlim(fsrc*0.8,fsrc*1.2); ax[i,1].set_title(f'spectrum slice at {tc:.2f}s, x = output f / r'); ax[i,1].legend(fontsize=6); ax[i,1].grid(alpha=.3)
plt.tight_layout(); plt.savefig('an/1985_evtraj.png',dpi=60)
