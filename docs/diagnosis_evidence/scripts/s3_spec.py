import sys; sys.path.insert(0,'an')
from load import *
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.signal import stft
n='1985'
fig,ax=plt.subplots(4,1,figsize=(16,14),sharex=True)
for i,v in enumerate(VERS):
    x,fs=load(n,v)
    f,t,X=stft(x,fs,'hann',nperseg=4096,noverlap=4096-256)
    S=20*np.log10(np.abs(X)+1e-9)
    ax[i].pcolormesh(t,f,S,vmin=S.max()-90,vmax=S.max(),shading='auto',cmap='magma')
    ax[i].set_yscale('symlog',linthresh=200); ax[i].set_ylim(20,20000); ax[i].set_title(v)
plt.tight_layout(); plt.savefig('an/1985_spec.png',dpi=70)
# low band zoom
fig,ax=plt.subplots(4,1,figsize=(16,12),sharex=True)
for i,v in enumerate(VERS):
    x,fs=load(n,v)
    f,t,X=stft(x,fs,'hann',nperseg=8192,noverlap=8192-256)
    S=20*np.log10(np.abs(X)+1e-9)
    ax[i].pcolormesh(t,f,S,vmin=S.max()-70,vmax=S.max(),shading='auto',cmap='magma')
    ax[i].set_ylim(20,400); ax[i].set_title(v+' (0-400 Hz, 8192 win)')
plt.tight_layout(); plt.savefig('an/1985_lowspec.png',dpi=70)
