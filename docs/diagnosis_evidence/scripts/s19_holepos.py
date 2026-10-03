from segs import *
n=sys.argv[1] if len(sys.argv)>1 else '1985'
fs=48000; X={v:load(n,v)[0] for v in VERS}
segs=transient_segments(n)
ons=np.load(f'an/{n}_onsets.npy')
win=256
def bbpow(v,lo=40,hi=16000):
    r=R[v]; y=bp(X[v],fs,lo*r,min(hi*r,23000),order=4); return np.convolve(y**2,np.ones(win)/win,'same')
P={v:bbpow(v) for v in VERS}
rows=[]
for S,E in segs:
    if E-S!=4096: continue
    o=[x*fs-S for x in ons if S<=x*fs<=E]
    if not o: continue
    op=o[0]
    for v in ['up_fifth','down_fifth']:
        D=10*np.log10((P[v][S:E+1024]+1e-12)/(P['original'][S:E+1024]+1e-12))
        Ds=np.convolve(D,np.ones(512)/512,'same')
        i=np.argmin(Ds[512:-512])+512
        rows.append((v,op,i,Ds[i], D[2560:3584].mean(), D[0:1536].mean()))
import collections
for v in ['up_fifth','down_fifth']:
    R_=np.array([(r[1],r[2],r[3],r[4],r[5]) for r in rows if r[0]==v])
    sl,ic=np.polyfit(R_[:,0],R_[:,1],1)
    cc=np.corrcoef(R_[:,0],R_[:,1])[0,1]
    print(f"{n} {v}: N={len(R_)}  hole-min position vs onset position: slope {sl:.2f}, r={cc:.2f}; hole pos median {np.median(R_[:,1]):.0f} IQR[{np.percentile(R_[:,1],25):.0f},{np.percentile(R_[:,1],75):.0f}]  onset pos IQR[{np.percentile(R_[:,0],25):.0f},{np.percentile(R_[:,0],75):.0f}]")
    print(f"    mean broadband ratio in blocks 2560-3584: median {np.median(R_[:,3]):+.1f} dB, frac < -3 dB: {(R_[:,3]<-3).mean():.2f}, frac < -6 dB: {(R_[:,3]<-6).mean():.2f};  min depth (512-smoothed) median {np.median(R_[:,2]):+.1f}")
