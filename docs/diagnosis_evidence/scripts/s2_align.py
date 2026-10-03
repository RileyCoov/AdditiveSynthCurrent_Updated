import sys; sys.path.insert(0,'an')
from load import *
import numpy as np
from scipy.signal import fftconvolve, butter, sosfiltfilt
def env(x, fs, hop=48):
    # broadband log envelope via rms in 2ms
    n=len(x)//hop; e=np.sqrt(np.mean(x[:n*hop].reshape(n,hop)**2,1)+1e-12); return e
def lag_env(a,b,fs,maxlag_ms=60,hop=48):
    ea=env(a,fs,hop); eb=env(b,fs,hop)
    ea=ea-ea.mean(); eb=eb-eb.mean()
    L=int(maxlag_ms/1000*fs/hop)
    c=[np.dot(ea[max(0,-l):len(ea)-max(0,l)], eb[max(0,l):len(eb)-max(0,-l)]) for l in range(-L,L+1)]
    return (np.argmax(c)-L)*hop/fs*1000
def clipstats(x, thr=0.95):
    a=np.abs(x); m=a.max()
    at=(a>=m-1e-4)
    # runs of consecutive samples at max
    runs=0; longest=0; cur=0
    for v in at:
        if v: cur+=1
        else:
            if cur>=2: runs+=1
            longest=max(longest,cur); cur=0
    return m, int(at.sum()), runs, longest
print(f"{'name':30s} {'u_lag_smp':>9s} {'u_corr':>7s} {'up_lag_ms':>9s} {'dn_lag_ms':>9s} | clip(max,n_at_max,runs>=2,longest) up / down / unity")
for n in NAMES:
    xo,fs=load(n,'original'); xu,_=load(n,'unity'); xup,_=load(n,'up_fifth'); xdn,_=load(n,'down_fifth')
    c=fftconvolve(xu, xo[::-1], 'full'); L=np.argmax(np.abs(c))-(len(xo)-1)
    corr=np.dot(xo,xu)/np.sqrt(np.dot(xo,xo)*np.dot(xu,xu))
    print(f"{n:30s} {L:9d} {corr:7.4f} {lag_env(xo,xup,fs):9.2f} {lag_env(xo,xdn,fs):9.2f} | {clipstats(xup)} {clipstats(xdn)} {clipstats(xu)}")
