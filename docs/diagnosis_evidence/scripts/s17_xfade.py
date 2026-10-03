from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fs=48000; xo,_=load('1985','original'); xu,_=load('1985','unity'); e=xu-xo
runs=[(53837,56663),(67150,70001),(141985,144750),(201227,204161),(226958,229741)]
fig,ax=plt.subplots(len(runs),1,figsize=(14,2.6*len(runs)))
for i,(a,b) in enumerate(runs):
    A=a-3000; B=b+3000; tt=np.arange(A,B)
    env=20*np.log10(np.convolve(np.abs(e[A:B]),np.ones(32)/32,'same')+1e-7)
    envx=20*np.log10(np.convolve(np.abs(xo[A:B]),np.ones(32)/32,'same')+1e-7)
    ax[i].plot(tt,envx,c='gray',lw=.6,label='|orig|'); ax[i].plot(tt,env,c='r',lw=.7,label='|unity-orig|')
    for k in range(A//256,B//256+1): ax[i].axvline(k*256,c='b',lw=0.2)
    for k in range(A//1024,B//1024+1): ax[i].axvline(k*1024,c='k',lw=0.6)
    ax[i].axvspan(a,b,color='g',alpha=.1); ax[i].set_ylim(-110,0); ax[i].legend(fontsize=6)
plt.tight_layout(); plt.savefig('an/1985_xfade.png',dpi=60)
# estimate where error first becomes small: find boundaries with 256 grid
for a,b in runs:
    pre=e[a-2048:a]; post=e[b:b+2048]
    # error energy per 256-block before and after
    pb=[20*np.log10(np.sqrt(np.mean(e[a-2048+j*256:a-2048+(j+1)*256]**2))+1e-9) for j in range(8)]
    qb=[20*np.log10(np.sqrt(np.mean(e[b+j*256:b+(j+1)*256]**2))+1e-9) for j in range(8)]
    print(a,b,'pre blocks dB',np.round(pb,1),' post',np.round(qb,1))
