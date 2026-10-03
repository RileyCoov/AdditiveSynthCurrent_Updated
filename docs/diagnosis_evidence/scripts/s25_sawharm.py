from common import *
from scipy.signal import hilbert
fs=48000
for n,f0 in [('440sawtooth',478.9116),('300hzSaw_5sec',300.0)]:
    for v in VERS:
        x,_=load(n,v); r=R[v]
        f0o=f0*r
        if v in ('up_fifth','down_fifth'):
            # refine f0 from spectrum
            N=65536; seg=x[len(x)//2-N//2:len(x)//2+N//2]*np.hanning(N); S=np.abs(np.fft.rfft(seg)); f=np.fft.rfftfreq(N,1/fs)
            sel=(f>f0o*0.97)&(f<f0o*1.03); f0o=f[sel][np.argmax(S[sel])]
        line=f"{n:14s} {v:10s} f0={f0o:7.2f} AMpkpk(dB)/AM@hop(dB) h1..h12:"
        for h in range(1,13):
            fh=h*f0o
            if fh>20000: break
            bw=min(60,f0o*0.4)
            y=bp(x,fs,fh-bw,fh+bw,order=4); a=np.abs(hilbert(y))
            mid=a[len(a)//4:3*len(a)//4]
            A=np.abs(np.fft.rfft((mid-mid.mean())*np.hanning(len(mid)))); fm=np.fft.rfftfreq(len(mid),1/fs)
            j=np.argmin(np.abs(fm-46.875)); m46=20*np.log10(A[j-3:j+4].max()*4/len(mid)/mid.mean()+1e-12)
            line+=f" {20*np.log10(mid.max()/mid.min()):4.1f}/{m46:5.0f}"
        print(line)
