from common import *
from scipy.signal import hilbert
n=sys.argv[1] if len(sys.argv)>1 else '1985'
fs=48000
X={v:load(n,v)[0] for v in VERS}
ons=np.load(f'an/{n}_onsets.npy') if os.path.exists(f'an/{n}_onsets.npy') else onsets(X['original'],fs,30,16000)[0]
extra=[1.350,2.977,4.174,4.739] if n=='1985' else []
E={}
for v in VERS:
    r=R[v]; y=bp(X[v],fs,35*r,160*r,order=6)
    env=np.abs(hilbert(y)); k=int(0.004*fs); env=np.convolve(env,np.ones(k)/k,'same')
    E[v]=20*np.log10(env+1e-9)
def dip(v,t):
    c=int(t*fs); base=np.median(E[v][c-int(0.08*fs):c-int(0.02*fs)])
    seg=E[v][c-int(0.01*fs):c+int(0.05*fs)]
    i=np.argmin(seg); depth=seg[i]-base
    # duration below base-6dB
    below=seg<base-6; dur=below.sum()/fs*1000
    return depth, (i/fs*1000-10), dur
rows=[]
for t in list(ons)+extra:
    if t<0.1 or t>len(X['original'])/fs-0.06: continue
    rows.append((t,)+tuple(dip(v,t) for v in VERS))
print(f"{n}: LF(35-160 Hz src, mapped) envelope dip in [-10,+50] ms vs median of [-80,-20] ms")
print(f"{'t':>6s} | {'orig depth/at/dur<-6':>22s} | {'unity':>22s} | {'up':>22s} | {'down':>22s}")
for r_ in rows:
    print(f"{r_[0]:6.3f} | "+" | ".join(f"{d:+6.1f}dB @{a:+5.1f}ms {du:4.1f}ms" for d,a,du in r_[1:]))
arr=np.array([[x[0] for x in r_[1:]] for r_ in rows]); dur=np.array([[x[2] for x in r_[1:]] for r_ in rows])
print("median depth:",np.round(np.median(arr,0),1)," mean dur<-6dB (ms):",np.round(dur.mean(0),1), " count depth<-10dB:",(arr<-10).sum(0))
