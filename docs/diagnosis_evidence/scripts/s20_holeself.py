from segs import *
fs=48000
def lfhole(n, lo=40, hi=160):
    X={v:load(n,v)[0] for v in VERS}; segs=[s for s in transient_segments(n) if s[1]-s[0]==4096]
    out={v:[] for v in VERS}
    P={}
    for v in VERS:
        r=R[v]; y=bp(X[v],fs,lo*r,hi*r,order=6); P[v]=np.convolve(y**2,np.ones(128)/128,'same')
    for S,E in segs:
        if S-2048<0 or E+2048>len(X['original']): continue
        for v in VERS:
            base=P[v][S-2048:S].mean()
            hole=P[v][S+1536:S+2304].mean()     # 768 smp window just before engine onset
            out[v].append(10*np.log10(hole/base+1e-12))
    return {v:np.array(out[v]) for v in VERS}, len(segs)
for n in ['1985','DrumLoopShort','HappyMono','PianoSampleMono','FairlightAdditiveWaveTA1_C3','FairlightAdditiveWaveTA1_C2','300hzSine_5sec']:
    for (lo,hi) in [(40,160),(500,2000)]:
        d,ns=lfhole(n,lo,hi)
        if ns==0: print(n,'no segments'); break
        o=d['original']
        s=f"{n:28s} {lo}-{hi}Hz segs={ns:2d} | self-ref level in [S+1536,S+2304] vs [S-2048,S]: "
        s+=" ".join(f"{v[:5]} {np.median(d[v]):+5.1f}" for v in VERS)
        s+=f" | (shift-orig) median: up {np.median(d['up_fifth']-o):+5.1f} dn {np.median(d['down_fifth']-o):+5.1f}; frac(dn-orig<-6dB) {np.mean(d['down_fifth']-o<-6):.2f}"
        print(s)
