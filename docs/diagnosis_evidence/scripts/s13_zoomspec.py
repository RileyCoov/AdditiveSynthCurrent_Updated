from tf import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
n=sys.argv[1]; t0,t1=float(sys.argv[2]),float(sys.argv[3]); f0,f1=float(sys.argv[4]),float(sys.argv[5]); out=sys.argv[6]
nfft=int(sys.argv[7]) if len(sys.argv)>7 else 4096
fs=48000
fig,ax=plt.subplots(1,4,figsize=(20,6),sharey=True)
for i,v in enumerate(VERS):
    x,_=load(n,v); a=max(0,int((t0-0.1)*fs)); b=int((t1+0.1)*fs)
    f,t,Z=stft(x[a:b],fs,'hann',nperseg=nfft,noverlap=nfft-128)
    S=20*np.log10(np.abs(Z)+1e-10); r=R[v]
    # plot in SOURCE-frequency coordinates (divide output freq by r)
    ax[i].pcolormesh(t+a/fs,f/r,S,vmin=S.max()-70,vmax=S.max(),shading='auto',cmap='magma')
    ax[i].set_ylim(f0,f1); ax[i].set_xlim(t0,t1); ax[i].set_title(f'{v} (y = output f / {r:.3f})')
plt.tight_layout(); plt.savefig(out,dpi=60)
