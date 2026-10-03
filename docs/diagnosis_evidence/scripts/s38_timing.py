from segs import *
fs=48000; n=sys.argv[1] if len(sys.argv)>1 else '1985'
X={v:load(n,v)[0] for v in VERS}
sg=[s for s in transient_segments(n) if s[1]-s[0]==4096 and s[0]>=3072 and s[0]+4096+3072<len(X['original'])]
pre=2048; post=4096+2048
for name,lo,hi in [('40-160',40,160),('500-2k',500,2000),('2k-6k',2000,6000),('6-16k',6000,16000)]:
    cur={}
    for v in VERS:
        r=R[v]; p=np.convolve(bp(X[v],fs,lo*r,min(hi*r,23500),order=6)**2,np.ones(240)/240,'same')
        cur[v]=10*np.log10(np.mean([p[S-pre:S+post] for S,E in sg],0))
    def rise(c):
        base=np.median(c[pre-1024:pre+512]); top=np.max(c[pre+1024:pre+3584]); half=base+0.5*(top-base)
        idx=np.where(c[pre+1024:pre+3584]>half)[0]; return idx[0]+1024 if len(idx) else -1
    ro=rise(cur['original'])
    line=f"{n} {name:7s}: original half-rise at S+{ro}"
    for v in ['up_fifth','down_fifth']:
        c=cur[v]; base=np.median(c[pre-1536:pre+0])
        dip=np.where(c[pre+1024:pre+3072]<base-2)[0]+1024
        rv=rise(c)
        line+=f" | {v[:5]}: half-rise S+{rv} ({(rv-ro)/fs*1000:+.1f} ms), dip(<pre-2dB) S+{dip.min() if len(dip) else -1}..S+{dip.max() if len(dip) else -1} ({len(dip)/fs*1000:.1f} ms), min {c[pre+1024:pre+3072].min()-base:+.1f} dB"
    print(line)
