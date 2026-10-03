from segs import *
fs=48000
def clicks(x, thr=8.0):
    # HP > 10 kHz, then |x| relative to local (10 ms) median envelope
    y=bp(x,fs,10000,23500,order=4); a=np.abs(y)
    from scipy.ndimage import median_filter
    k=481; med=median_filter(a,size=k,mode='nearest')+1e-6
    z=a/med
    idx=np.where(z>thr)[0]
    # merge within 64
    ev=[]
    for i in idx:
        if not ev or i-ev[-1][-1]>64: ev.append([i])
        else: ev[-1].append(i)
    return [(e[0],max(z[e[0]:e[-1]+1])) for e in ev]
for n in ['1985','HappyMono','take-me-out','DrumLoopShort','SaintSaensMonoChord','choir-burst_C_major_Mono_48k','Female_Sung_Line_3_48','FairlightAdditiveWaveTA1_C3','300hzSaw_5sec','300hzSine_5sec']:
    line=f"{n:28s}"
    for v in VERS:
        x,_=load(n,v); c=clicks(x)
        ph=np.array([i%1024 for i,_ in c])
        near=np.mean((ph<64)|(ph>960)) if len(c) else np.nan
        line+=f" | {v[:5]} {len(c)/(len(x)/fs):5.1f}/s gridfrac {near:.2f}"
    print(line)
print("(click = HP>10k sample magnitude > 8x local 10 ms median; gridfrac = fraction within ±64 smp of k*1024; chance = 0.125)")
