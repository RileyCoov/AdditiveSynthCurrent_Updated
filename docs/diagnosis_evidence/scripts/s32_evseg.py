from segs import *
import pickle
fs=48000
for n in ['1985','HappyMono','DrumLoopShort']:
    ev=pickle.load(open(f'an/{n}_events.pkl','rb'))
    segs=transient_segments(n); L=len(load(n,'original')[0])
    inseg=np.zeros(L,bool)
    for S,E in segs: inseg[S:E]=True
    cov=inseg.mean()
    for v in ['up_fifth','down_fifth']:
        for key in [v, v+'_sal']:
            E_=ev[key]
            # event spans window centres t0..t1; the synthesized artifact extends ±2048 samples (4096 analysis window)
            c=np.array([int(0.5*(e['t0']+e['t1'])*fs) for e in E_]); c=c[(c>=0)&(c<L)]
            frac=inseg[c].mean() if len(c) else np.nan
            print(f"{n:14s} {key:16s} N={len(c):4d} frac centred inside transient segments {frac:.2f} (time coverage {cov:.2f})")
    # pre-echo: HF 6-16k ratio in [S+512,S+1536]
    X={v:load(n,v)[0] for v in VERS}
    for v in ['up_fifth','down_fifth']:
        r=R[v]; po=bp(X['original'],fs,6000,16000)**2; pv=bp(X[v],fs,6000*r,min(16000*r,23500))**2
        d=[10*np.log10(pv[S+512:S+1536].mean()/po[S+512:S+1536].mean()) for S,E in segs if E-S==4096]
        d2=[10*np.log10(pv[S-4096:S-1024].mean()/po[S-4096:S-1024].mean()) for S,E in segs if E-S==4096 and S>4096]
        print(f"   {v}: HF(6-16k) level vs original in [S+512,S+1536] (pre-onset): median {np.median(d):+.1f} dB; baseline [S-4096,S-1024]: {np.median(d2):+.1f} dB")
