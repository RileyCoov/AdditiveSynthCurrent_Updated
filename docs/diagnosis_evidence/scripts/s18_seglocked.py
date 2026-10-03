from segs import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
n=sys.argv[1] if len(sys.argv)>1 else '1985'
fs=48000; X={v:load(n,v)[0] for v in VERS}
segs=transient_segments(n)
L=[E-S for S,E in segs]
print(n,"transient segments:",len(segs)," lengths (hops):",np.bincount(np.array(L)//1024))
ons=np.load(f'an/{n}_onsets.npy') if os.path.exists(f'an/{n}_onsets.npy') else None
if ons is not None:
    rel=[]
    for S,E in segs:
        inside=[o*fs-S for o in ons if S-2048<=o*fs<=E+2048]
        rel+=inside
    print(" onset position relative to segment start (samples): median %.0f, IQR [%.0f, %.0f]"%(np.median(rel),np.percentile(rel,25),np.percentile(rel,75)))
bands={'40-160':(40,160),'160-500':(160,500),'500-2k':(500,2000),'2k-6k':(2000,6000),'6k-16k':(6000,16000)}
pre=3072; post=4096+3072
tt=np.arange(-pre,post)
fig,ax=plt.subplots(len(bands),1,figsize=(12,14),sharex=True)
win=240
for bi,(bn,(lo,hi)) in enumerate(bands.items()):
    curves={}
    for v in VERS:
        r=R[v]; y=bp(X[v],fs,lo*r,hi*r,order=6); p=np.convolve(y**2,np.ones(win)/win,'same')
        acc=[]
        for S,E in segs:
            if E-S!=4096 or S-pre<0 or S+post>len(p): continue
            acc.append(p[S-pre:S+post])
        curves[v]=np.mean(acc,0)   # average POWER (not dB) across segments
        ax[bi].plot(tt,10*np.log10(curves[v]+1e-12),label=v,lw=.8)
    for k in range(-3,8): ax[bi].axvline(k*1024,c='k',lw=.4)
    for k in range(0,17): ax[bi].axvline(k*256,c='b',lw=.15)
    ax[bi].set_title(f'{n} band {bn} Hz (src, mapped), mean power over {len(acc)} 4-hop transient segments; x=samples from segment start')
    ax[bi].legend(fontsize=6)
    # numeric: ratio to original per 256 block from -1024 to 5120
    line=f"{bn:8s}"
    for v in VERS[1:]:
        rr=[10*np.log10(curves[v][pre+j*256:pre+(j+1)*256].mean()/curves['original'][pre+j*256:pre+(j+1)*256].mean()) for j in range(-4,20)]
        line+=f"\n   {v:10s} "+" ".join(f"{x:+5.1f}" for x in rr)
    print(line)
print("columns: 256-sample blocks from -1024 to +5120 relative to segment start (segment = blocks 0..15)")
plt.tight_layout(); plt.savefig(f'an/{n}_seglocked.png',dpi=60)
