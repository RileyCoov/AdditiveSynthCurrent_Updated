from segs import *
fs=48000
def fold_ratio(n, lo, hi, excl_segs=True):
    X={v:load(n,v)[0] for v in VERS}
    segs=transient_segments(n)
    mask=np.ones(len(X['original']),bool)
    for S,E in segs: mask[max(0,S-2048):E+2048]=False
    out={}
    for v in VERS[1:]:
        r=R[v]; y=bp(X[v],fs,lo*r,min(hi*r,23000),order=4)**2; yo=bp(X['original'],fs,lo,hi,order=4)**2
        ph=np.arange(len(y))%1024
        m=mask.copy(); m[:4096]=False; m[-4096:]=False
        num=np.bincount(ph[m]//64,weights=y[m],minlength=16); den=np.bincount(ph[m]//64,weights=yo[m],minlength=16)
        rr=10*np.log10(num/den); out[v]=rr-rr.mean()
    return out
for n in ['1985','HappyMono','take-me-out','SaintSaensMonoChord','choir-burst_C_major_Mono_48k','Female_Sung_Line_3_48','300hzSaw_5sec','440sawtooth']:
    o=fold_ratio(n,40,16000)
    print(f"{n:28s} power ratio vs position mod 1024 (64-smp bins, dB, mean-removed): up pk-pk {np.ptp(o['up_fifth']):.2f} dn pk-pk {np.ptp(o['down_fifth']):.2f} unity pk-pk {np.ptp(o['unity']):.2f}")
    print("     up:",np.round(o['up_fifth'],2).tolist()); print("     dn:",np.round(o['down_fifth'],2).tolist())
