from segs import *
fs=48000
for n in ['1985','HappyMono','DrumLoopShort','FairlightAdditiveWaveTA1_C3','FairlightAdditiveWaveTA1_C2']:
    X={v:load(n,v)[0] for v in VERS}; dur=len(X['original'])/fs
    segs=[s for s in transient_segments(n) if s[1]-s[0]==4096 and s[0]>=2048 and s[1]+2048<len(X['original'])]
    tot=np.mean(X['original']**2)
    line=f"{n:28s} segs/s {len(segs)/dur:4.2f} |"
    for lo,hi in [(40,160),(160,500),(500,2000)]:
        share=10*np.log10(np.mean(bp(X['original'],fs,lo,hi)**2)/tot)
        P={}
        for v in VERS:
            r=R[v]; P[v]=np.convolve(bp(X[v],fs,lo*r,hi*r,order=6)**2,np.ones(128)/128,'same')
        steady=0; exp_=[]
        for S,E in segs:
            ch_o=10*np.log10(P['original'][S+1536:S+2304].mean()/P['original'][S-2048:S].mean())
            if abs(ch_o)<2:
                steady+=1
                d=np.mean([10*np.log10(P[v][S+1536:S+2304].mean()/P[v][S-2048:S].mean())-ch_o for v in ['up_fifth','down_fifth']])
                exp_.append(d)
        line+=f" {lo}-{hi}: share {share:5.1f} dB, steady-bg segs {steady:2d}/{len(segs)}, dip in steady bg {np.median(exp_) if exp_ else np.nan:+5.1f} dB, exposed holes/s {steady/dur:4.2f} |"
    print(line)
