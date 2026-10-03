from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
n=sys.argv[1] if len(sys.argv)>1 else '1985'
xo,fs=load(n,'original')
X={v:load(n,v)[0] for v in VERS}
o,_=onsets(xo,fs,30,16000)
# refine onsets with 1 ms envelope of broadband (HP 100) signal: steepest log-rise
hp=bp(xo,fs,100,20000); e1=rms_env(hp,48); le=20*np.log10(e1)
ref=[]
for t in o:
    a=int((t-0.03)*fs); b=int((t+0.01)*fs); a=max(a,48); b=min(b,len(xo)-48)
    d=le[a+48:b+48]-le[a:b]   # 1ms-difference
    ref.append((a+np.argmax(d))/fs)
ref=np.array(ref)
np.save(f'an/{n}_onsets.npy',ref)
print("refined onsets:",np.round(ref,3).tolist())
bands={'LF 40-200':(40,200),'LM 200-800':(200,800),'M 800-3k':(800,3000),'HF 3k-10k':(3000,10000)}
pre,post=0.08,0.15
tt=np.arange(-int(pre*fs),int(post*fs))/fs
fig,ax=plt.subplots(len(bands),1,figsize=(10,12))
summary={}
for bi,(bn,(lo,hi)) in enumerate(bands.items()):
    for v in VERS:
        r=R[v]; lo2,hi2=lo*r,hi*r
        y=bp(X[v],fs,lo2,hi2); env=rms_env(y,240)  # 5 ms
        segs=[]
        for t in ref:
            c=int(t*fs); a=c-int(pre*fs); b=c+int(post*fs)
            if a<0 or b>len(env): continue
            segs.append(20*np.log10(env[a:b]))
        segs=np.array(segs); m=segs.mean(0)
        summary[(bn,v)]=segs
        ax[bi].plot(tt*1000,m,label=v)
    ax[bi].set_title(f'{n} {bn} (mapped band for shifted), mean over {len(segs)} onsets'); ax[bi].axvline(0,c='k',lw=0.5); ax[bi].legend(fontsize=7); ax[bi].grid(alpha=.3)
plt.tight_layout(); plt.savefig(f'an/{n}_eventbands.png',dpi=70)
# numeric: per band, mean dB of (version - original) in windows
wins=[(-60,-20),(-20,0),(0,10),(10,30),(30,60),(60,120)]
for bn in bands:
    o_=summary[(bn,'original')]
    line=f"{bn:12s}"
    for v in VERS[1:]:
        d=summary[(bn,v)]-o_
        line+=f" | {v[:5]}: "+" ".join(f"{np.mean(d[:, (tt*1000>=a)&(tt*1000<b)]):+5.1f}" for a,b in wins)
    print(line)
print("windows(ms):",wins)
