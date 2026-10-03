from segs import *
from tf import spec
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.edgecolor':'#9a9993','axes.labelcolor':'#52514e','xtick.color':'#52514e','ytick.color':'#52514e','axes.titlecolor':'#0b0b0b','figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
C={'original':'#0b0b0b','unity':'#2a78d6','up_fifth':'#eb6834','down_fifth':'#1baf7a'}
LBL={'original':'Original','unity':'Unity','up_fifth':'Up a fifth','down_fifth':'Down a fifth'}
fs=48000; n='1985'
X={v:load(n,v)[0] for v in VERS}
sg=[s for s in transient_segments(n) if s[1]-s[0]==4096 and s[0]>=3072 and s[0]+4096+3072<len(X['original'])]
pre=2048; post=4096+2048; tt=(np.arange(-pre,post))/fs*1000
bands=[('Bass 40–160 Hz',40,160),('Mids 500–2000 Hz',500,2000),('Air 6–16 kHz',6000,16000)]
fig,ax=plt.subplots(3,1,figsize=(11,10),sharex=True)
for i,(name,lo,hi) in enumerate(bands):
    for v in ['original','unity','up_fifth','down_fifth']:
        r=R[v]; p=np.convolve(bp(X[v],fs,lo*r,min(hi*r,23500),order=6)**2,np.ones(240)/240,'same')
        m=np.mean([p[S-pre:S+post] for S,E in sg],0)
        ax[i].plot(tt,10*np.log10(m),c=C[v],lw=2 if v!='unity' else 1.2,ls='-' if v!='unity' else (0,(4,2)),label=LBL[v])
    ax[i].axvspan(1536/fs*1000,2304/fs*1000,color='#e34948',alpha=0.08,lw=0)
    ax[i].axvspan(512/fs*1000,1536/fs*1000,color='#eda100',alpha=0.08,lw=0)
    ax[i].axvline(2048/fs*1000,c='#9a9993',lw=1,ls=':')
    for k in range(-2,7): ax[i].axvline(k*1024/fs*1000,c='#e6e5e0',lw=0.6,zorder=0)
    ax[i].set_ylabel('Mean band power (dB)'); ax[i].set_title(name,loc='left',fontsize=11)
    ax[i].spines[['top','right']].set_visible(False)
ax[0].legend(loc='upper left',frameon=False,ncol=4,fontsize=9)
ax[0].text(1920/fs*1000,ax[0].get_ylim()[1],'dropout\n(~16 ms)',ha='center',va='top',fontsize=8,color='#a33')
ax[2].text(1024/fs*1000,ax[2].get_ylim()[1],'pre-echo',ha='center',va='top',fontsize=8,color='#8a6400')
ax[2].set_xlabel('ms from start of the engine’s 4-hop transient segment (gridlines = 1024-sample hops; dotted = segment midpoint ≈ hit)')
fig.suptitle(f'1985: average of all {len(sg)} transient segments, locked to the engine frame grid\nUnity tracks the original exactly (passthrough); both shifts show HF pre-echo, then a broadband dip right before the hit',x=0.01,ha='left',fontsize=12)
plt.tight_layout(); plt.savefig('deliver/fig1_transient_segment_dropout.png',dpi=110)

cases=[(1.05,1.55,530,'≈530 Hz, 1.24–1.30 s'),(1.30,1.85,1205,'≈1205 Hz, 1.50–1.58 s'),(4.20,4.90,3441,'≈3441 Hz, 4.37–4.69 s')]
S={}
for v in VERS:
    f,t,Y=spec(X[v],fs); S[v]=(f,t,20*np.log10(Y+1e-9))
fig,ax=plt.subplots(1,3,figsize=(15,4.6))
for i,(t0,t1,fsrc,lab) in enumerate(cases):
    for v in ['original','unity','up_fifth','down_fifth']:
        f,t,L=S[v]; r=R[v]; sel=(t>=t0)&(t<=t1); band=(f>=r*fsrc*0.985)&(f<=r*fsrc*1.015)
        y=L[band][:,sel].max(0); ref=S['original'][2][:, (S['original'][1]>=t0)&(S['original'][1]<=t1)].max()
        ax[i].plot(t[sel],y-ref,c=C[v],lw=2 if v!='unity' else 1.2,ls='-' if v!='unity' else (0,(4,2)),label=LBL[v])
    ax[i].set_title(f'Source frequency {lab}',loc='left',fontsize=10); ax[i].set_xlabel('time (s)')
    ax[i].spines[['top','right']].set_visible(False); ax[i].grid(alpha=0.25)
ax[0].set_ylabel('Level at that source frequency\n(dB re loudest component in window)')
ax[0].legend(frameon=False,fontsize=8,loc='lower left')
fig.suptitle('1985: the loudest “brief tones” are the same partial in both directions, about 15 dB above the source. Unity matches the original.',x=0.01,ha='left',fontsize=11.5)
plt.tight_layout(); plt.savefig('deliver/fig2_spurious_partials.png',dpi=110)
