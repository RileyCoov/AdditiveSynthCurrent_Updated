from segs import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
n='1985'; fs=48000; X={v:load(n,v)[0] for v in VERS}
segs=[s for s in transient_segments(n) if s[1]-s[0]==4096]
pick=[s for s in segs if abs(s[0]/fs-2.944)<0.01 or abs(s[0]/fs-4.715)<0.01 or abs(s[0]/fs-1.387)<0.01 or abs(s[0]/fs-6.0)<0.1]
fig,ax=plt.subplots(len(pick),2,figsize=(18,3.2*len(pick)))
for i,(S,E) in enumerate(pick):
    A=S-1536; B=S+4096+1024; tt=np.arange(A,B)-S
    for v,c,off in zip(VERS,['k','tab:orange','tab:green','tab:red'],[0.9,0.3,-0.3,-0.9]):
        ax[i,0].plot(tt,X[v][A:B]+off,c=c,lw=.5,label=v)
        r=R[v]; y=bp(X[v],fs,35*r,160*r,order=6)
        ax[i,1].plot(tt,y[A:B]*2+off*0.5,c=c,lw=.7,label=v)
    for a in ax[i]:
        for k in range(-1,6): a.axvline(k*1024,c='gray',lw=.5)
        a.set_title(f'segment S={S} ({S/fs:.3f}s); x = samples from S',fontsize=8); a.legend(fontsize=5,loc='lower left')
    ax[i,0].set_ylabel('full band'); ax[i,1].set_ylabel('LF 35-160 src (x2)')
plt.tight_layout(); plt.savefig('an/1985_segwave.png',dpi=60)
