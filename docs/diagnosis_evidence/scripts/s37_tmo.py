from common import *
fs=48000
for n in ['take-me-out','1985']:
    X={v:load(n,v)[0] for v in VERS}
    o,_=onsets(X['original'],fs,30,16000)
    hp=bp(X['original'],fs,100,20000); le=20*np.log10(rms_env(hp,48))
    ref=[]
    for t in o:
        a=max(int((t-0.03)*fs),48); b=min(int((t+0.01)*fs),len(le)-48); d=le[a+48:b+48]-le[a:b]; ref.append(a+np.argmax(d))
    ref=[c for c in ref if c>4096 and c<len(hp)-4096]
    line=f"{n:12s} onsets {len(ref)}:"
    for lo,hi in [(40,160),(500,2000),(6000,16000)]:
        Pw={v:np.convolve(bp(X[v],fs,lo*R[v],min(hi*R[v],23500),order=6)**2,np.ones(128)/128,'same') for v in VERS}
        res={}
        for v in VERS:
            res[v]=np.median([10*np.log10(Pw[v][c-600:c-50].mean()/Pw[v][c-2600:c-800].mean()) for c in ref])
        line+=f" | {lo}-{hi}: pre-onset[-600,-50] vs [-2600,-800]: orig {res['original']:+.1f} up {res['up_fifth']:+.1f} dn {res['down_fifth']:+.1f}"
    print(line)
