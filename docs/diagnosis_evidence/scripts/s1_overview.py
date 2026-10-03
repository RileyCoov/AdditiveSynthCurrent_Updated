import sys; sys.path.insert(0,'an')
from load import *
import numpy as np
from scipy.signal import correlate
print(f"{'name':38s} {'ver':10s} {'fs':>6s} {'ch':>2s} {'sub':>8s} {'N':>8s} {'dur':>6s} {'peak':>7s} {'rms_dB':>7s}")
for n in NAMES:
    for v in VERS:
        i=info(n,v); x,fs=load(n,v)
        print(f"{n:38s} {v:10s} {fs:6d} {i.channels:2d} {i.subtype:>8s} {len(x):8d} {len(x)/fs:6.2f} {np.abs(x).max():7.4f} {20*np.log10(np.sqrt(np.mean(x**2))+1e-12):7.2f}")
