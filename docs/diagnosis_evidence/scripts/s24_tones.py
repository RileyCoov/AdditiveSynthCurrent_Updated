from common import *
from scipy.signal import hilbert, get_window
fs=48000
def tone_track(x, f0, fs, bw=40):
    y=bp(x,fs,f0-bw,f0+bw,order=4); a=hilbert(y); amp=np.abs(a); ph=np.unwrap(np.angle(a))
    inst=np.diff(ph)*fs/(2*np.pi)
    return amp, inst
def spectrum_spurs(x,fs,f0,nh,N=65536):
    seg=x[len(x)//2-N//2:len(x)//2+N//2]*get_window('blackmanharris',N)
    S=np.abs(np.fft.rfft(seg)); f=np.fft.rfftfreq(N,1/fs); S=20*np.log10(S/S.max()+1e-12)
    # mask harmonics ±15 Hz
    mask=np.ones_like(S,bool)
    for h in range(1,nh+1):
        mask&=np.abs(f-h*f0)>15
    mask&=f>20
    i=np.argmax(np.where(mask,S,-300)); 
    # top 5 spurs
    Sm=np.where(mask,S,-300); idx=np.argsort(Sm)[::-1]
    spurs=[];
    for j in idx:
        if all(abs(f[j]-f[k])>10 for k in spurs): spurs.append(j)
        if len(spurs)>=5: break
    return [(round(f[j],1),round(S[j],1)) for j in spurs]
for n,f0,nh in [('300hzSine_5sec',300,1),('300hzSaw_5sec',300,79),('440sawtooth',440,54)]:
    for v in VERS:
        x,_=load(n,v); r=R[v]; f=f0*r
        amp,inst=tone_track(x,f,fs)
        mid=slice(len(x)//4,3*len(x)//4)
        a=amp[mid]; fi=inst[mid]
        am_db=20*np.log10(a.max()/a.min())
        # AM spectrum: modulation at hop rate 46.875
        A=np.abs(np.fft.rfft((a-a.mean())*np.hanning(len(a)))); fm=np.fft.rfftfreq(len(a),1/fs)
        def pk(fr): j=np.argmin(np.abs(fm-fr)); return 20*np.log10(A[j-3:j+4].max()/a.mean()/len(a)*4+1e-12)
        nhh=int(min(nh, 20000/f))
        print(f"{n:16s} {v:10s} f_out={f:7.2f} | mean IF {np.mean(fi):8.3f} Hz (err {np.mean(fi)-f:+.3f}) IF std {np.std(fi):6.3f} Hz | AM pk-pk {am_db:5.2f} dB, AM@46.9Hz {pk(46.875):6.1f} dB, @93.75 {pk(93.75):6.1f} | spurs(rel dB): {spectrum_spurs(x,fs,f,nhh)}")
