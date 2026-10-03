from common import *
fs=48000
bands=[(30,80),(80,200),(200,500),(500,1500),(1500,4000),(4000,8000),(8000,15000)]
print("Long-term band level difference vs original (shifted: source band mapped by r), dB")
print(f"{'name':30s} {'ver':6s} "+" ".join(f"{lo}-{hi}".rjust(10) for lo,hi in bands)+"   total")
for n in NAMES:
    X={v:load(n,v)[0] for v in VERS}
    po=[np.mean(bp(X['original'],fs,lo,hi)**2) for lo,hi in bands]
    for v in VERS[1:]:
        r=R[v]
        pv=[np.mean(bp(X[v],fs,lo*r,min(hi*r,23500))**2) for lo,hi in bands]
        d=[10*np.log10(a/b) if b>1e-12 else np.nan for a,b in zip(pv,po)]
        tot=10*np.log10(np.mean(X[v]**2)/np.mean(X['original']**2))
        print(f"{n:30s} {v[:6]:6s} "+" ".join(f"{x:10.1f}" for x in d)+f"   {tot:+5.1f}")
