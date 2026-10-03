# Re-derives every headline number from the WAVs in one pass. Output -> KEY_NUMBERS.txt
from segs import *
from tf import *
import pickle
fs=48000; LSB=1/32768
out=[]
def P(s): print(s); out.append(s)
# 1. passthrough
for n in ['1985','DrumLoopShort','HappyMono','take-me-out']:
    xo,_=load(n,'original'); xu,_=load(n,'unity'); e=xu-xo
    k=512; pad=np.pad(xo**2,(k//2,k//2),mode='edge'); loc=np.sqrt(np.convolve(pad,np.ones(k)/k,'valid')[:len(xo)]); act=loc>10**(-50/20)
    ex=np.abs(e)<=1.01*LSB; sg=transient_segments(n)
    P(f"[1] {n}: unity bit-exact (|e|<=1 LSB) share of active samples {100*(ex&act).sum()/act.sum():.1f}%; 4-hop grid segments {sum(1 for s in sg if s[1]-s[0]==4096)}/{len(sg)}; unity SER {10*np.log10(np.sum(xo**2)/np.sum(e**2)):.2f} dB")
# 2. segment-locked hole and pre-echo
for n in ['1985','HappyMono','DrumLoopShort']:
    X={v:load(n,v)[0] for v in VERS}; sg=[s for s in transient_segments(n) if s[1]-s[0]==4096 and s[0]>=4096 and s[1]+2048<len(X['original'])]
    for lo,hi in [(40,160),(500,2000)]:
        Pw={v:np.convolve(bp(X[v],fs,lo*R[v],hi*R[v],order=6)**2,np.ones(128)/128,'same') for v in VERS}
        d={v:np.array([10*np.log10(Pw[v][S+1536:S+2304].mean()/Pw[v][S-2048:S].mean()) for S,E in sg]) for v in VERS}
        P(f"[2] {n} {lo}-{hi} Hz: self-ref level [S+1536,S+2304] vs [S-2048,S], median (shift - original): up {np.median(d['up_fifth']-d['original']):+.1f} dB, down {np.median(d['down_fifth']-d['original']):+.1f} dB; unity {np.median(d['unity']-d['original']):+.1f}  (N={len(sg)})")
    po=bp(X['original'],fs,6000,16000)**2
    for v in ['up_fifth','down_fifth']:
        pv=bp(X[v],fs,6000*R[v],min(16000*R[v],23500))**2
        a=np.median([10*np.log10(pv[S+512:S+1536].mean()/po[S+512:S+1536].mean()) for S,E in sg])
        b=np.median([10*np.log10(pv[S-4096:S-1024].mean()/po[S-4096:S-1024].mean()) for S,E in sg])
        P(f"[3] {n} {v}: HF 6-16k pre-onset [S+512,S+1536] {a:+.1f} dB vs original; between-hit baseline {b:+.1f} dB")
# 4. symmetry
ev=pickle.load(open('an/1985_events.pkl','rb'))
def top(E,N): return sorted(E,key=lambda e:-(e['exc']+2*np.log2(e['nfr'])+0.3*e['rel_to_frame']))[:N]
def match(A,ra,B,rb):
    return sum(1 for a in A if any(abs(a['f']/ra-b['f']/rb)/(a['f']/ra)<0.03 and a['t0']<=b['t1']+0.02 and b['t0']<=a['t1']+0.02 for b in B))
if not ev['up_fifth']: P("[4] 1985: no up-shift tonal events -> symmetry test n/a")
else:
  A=top(ev['up_fifth'],25); h=match(A,R['up_fifth'],ev['down_fifth'],R['down_fifth'])
  rng=np.random.default_rng(0); T=max(e['t1'] for e in ev['up_fifth']); nl=[]
  for k in range(500):
      s=rng.uniform(0.3,T-0.3); A2=[dict(e,t0=(e['t0']+s)%T,t1=(e['t0']+s)%T+(e['t1']-e['t0'])) for e in A]; nl.append(match(A2,R['up_fifth'],ev['down_fifth'],R['down_fifth']))
  P(f"[4] 1985 top-25 up-shift tonal events with a down-shift twin at same source freq (±3%) and overlapping time: {h}/25; circular-shift null mean {np.mean(nl):.1f}, 99th pct {np.percentile(nl,99):.0f}, p≈{(np.sum(np.array(nl)>=h)+1)/(len(nl)+1):.3f}")
for n in ['1985','HappyMono','take-me-out']:
    ev=pickle.load(open(f'an/{n}_events.pkl','rb')); dur=len(load(n,'original')[0])/fs
    for v in ['up_fifth','down_fifth']:
        E=ev[v]; r=R[v]; rel=np.array([e['rel_to_frame'] for e in E]+[np.nan])[:len(E)]; f=np.array([e['f']/r for e in E]+[np.nan])[:len(E)]; pr=rel>-12
        P(f"[5] {n} {v}: all events {len(E)} ({len(E)/dur:.0f}/s), {100*np.mean(f>4000):.0f}% src>4 kHz, median rel-to-frame-max {np.median(rel):.1f} dB | prominent (within 12 dB of frame max) {pr.sum()} ({pr.sum()/dur:.1f}/s), median src f {np.median(f[pr]):.0f} Hz; unity events {len(ev['unity'])}")
open('KEY_NUMBERS.txt','w').write("\n".join(out)+"\n")
