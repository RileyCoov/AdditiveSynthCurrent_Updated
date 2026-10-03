from common import *
def steepest(le, c, w0, w1, d=48):
    a=c+w0; b=c+w1
    dd=le[a+d:b+d]-le[a:b]
    return a+np.argmax(dd), dd.max()
for n in ['1985','DrumLoopShort','HappyMono','take-me-out']:
    X={v:load(n,v)[0] for v in VERS}; fs=48000
    xo=X['original']
    o,_=onsets(xo,fs,30,16000)
    print(n, len(o),"onsets")
    for bn,(lo,hi) in {'40-200':(40,200),'200-800':(200,800),'800-3k':(800,3000),'3k-10k':(3000,10000)}.items():
        les={}
        for v in VERS:
            r=R[v]; y=bp(X[v],fs,lo*r,hi*r); les[v]=20*np.log10(rms_env(y,96))
        # reference onset in original band
        line=f"  {bn:8s}"
        offs={v:[] for v in VERS[1:]}; rises={v:[] for v in VERS}
        for t in o:
            c=int(t*fs)
            if c<4000 or c>len(xo)-4000: continue
            i0,r0=steepest(les['original'],c,-1440,480)
            if r0<3: continue  # need a clear onset in this band
            rises['original'].append(r0)
            for v in VERS[1:]:
                i1,r1=steepest(les[v],i0,-1440,1440); offs[v].append((i1-i0)/fs*1000); rises[v].append(r1)
        for v in VERS[1:]:
            a=np.array(offs[v]); 
            line+=f" | {v[:5]} N={len(a)} med {np.median(a) if len(a) else np.nan:+5.1f} ms IQR[{np.percentile(a,25) if len(a) else 0:+5.1f},{np.percentile(a,75) if len(a) else 0:+5.1f}] rise {np.median(rises[v]) if len(a) else 0:4.1f}dB"
        line+=f" | orig rise {np.median(rises['original']) if rises['original'] else 0:4.1f}dB"
        print(line)
