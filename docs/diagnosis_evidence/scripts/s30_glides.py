from common import *
from scipy.ndimage import median_filter
fs=48000; N=4096; H=256
def peaks_tracks(x, fmin=150, fmax=8000):
    f,t,Z=stft(x,fs,'hann',nperseg=N,noverlap=N-H,boundary=None,padded=False)
    M=np.abs(Z); L=20*np.log10(M+1e-12)
    med=median_filter(L,size=(31,1),mode='nearest')
    df=fs/N
    frames=[]
    gmax=L.max()
    for j in range(L.shape[1]):
        l=L[:,j]; m=med[:,j]
        k=np.where((l[1:-1]>l[:-2])&(l[1:-1]>=l[2:])&(l[1:-1]-m[1:-1]>12)&(l[1:-1]>gmax-70))[0]+1
        k=k[(f[k]>=fmin)&(f[k]<=fmax)]
        a,b,c=l[k-1],l[k],l[k+1]; d=0.5*(a-c)/(a-2*b+c)
        frames.append((f[k]+d*df, 10**((b-0.25*(a-c)*d)/20)))
    # link
    tracks=[]; active={}
    for j,(fr,am) in enumerate(frames):
        new={}
        used=set()
        for tid,(lastf,lst) in active.items():
            if len(fr)==0: continue
            i=np.argmin(np.abs(fr-lastf))
            if abs(fr[i]-lastf)<=12 and i not in used:
                lst.append((j,fr[i],am[i])); used.add(i); new[tid]=(fr[i],lst)
            else: tracks.append(lst)
        for i in range(len(fr)):
            if i not in used:
                tid=len(tracks)+len(new)+j*100000+i; new[tid]=(fr[i],[(j,fr[i],am[i])])
        active=new
    tracks+= [lst for _,lst in active.values()]
    return tracks
def glide_stats(x):
    tr=peaks_tracks(x)
    E_tot=0; E_fast=0; E_vfast=0; slopes=[]
    for lst in tr:
        if len(lst)<5: continue
        j=np.array([p[0] for p in lst]); fr=np.array([p[1] for p in lst]); am=np.array([p[2] for p in lst])
        # local slope via 5-point regression windows
        for k in range(2,len(lst)-2):
            s=np.polyfit(j[k-2:k+3]*H/fs, fr[k-2:k+3],1)[0]   # Hz/s
            e=am[k]**2; E_tot+=e
            perhop=abs(s)*1024/fs
            if perhop>25: E_fast+=e
            if perhop>50: E_vfast+=e
            slopes.append((fr[k],perhop,e))
    return E_fast/E_tot, E_vfast/E_tot, slopes
for n in ['1985','HappyMono','take-me-out','Female_Sung_Line_3_48','choir-burst_C_major_Mono_48k','SaintSaensMonoChord','PianoSampleMono','DrumLoopShort','300hzSaw_5sec']:
    x,_=load(n,'original')
    a,b,sl=glide_stats(x)
    sl=np.array(sl)
    hi=sl[sl[:,0]>1500]
    fhi=(hi[:,2]*(hi[:,1]>25)).sum()/hi[:,2].sum() if len(hi) else np.nan
    print(f"{n:30s} tonal energy (150-8000 Hz) on partials moving >25 Hz/hop: {100*a:5.1f}%   >50 Hz/hop: {100*b:5.1f}%   | above 1.5 kHz: {100*fhi:5.1f}%")
