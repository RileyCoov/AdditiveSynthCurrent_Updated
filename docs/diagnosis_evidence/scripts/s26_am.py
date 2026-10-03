from common import *
from scipy.signal import hilbert
fs=48000
for v,f0 in [('original',478.91),('unity',478.91),('up_fifth',717.77),('down_fifth',319.34)]:
    x,_=load('440sawtooth',v)
    for h in [1,4]:
        fh=h*f0; y=bp(x,fs,fh-60,fh+60,order=4); a=np.abs(hilbert(y)); mid=a[len(a)//6:5*len(a)//6]
        A=np.abs(np.fft.rfft((mid-mid.mean())*np.hanning(len(mid))))*4/len(mid)/mid.mean(); fm=np.fft.rfftfreq(len(mid),1/fs)
        sel=(fm>2)&(fm<300); idx=np.argsort(A[sel])[::-1][:6]
        print(v,'h',h,'AM components (Hz, dB):',[(round(fm[sel][i],1),round(20*np.log10(A[sel][i]),1)) for i in sorted(idx,key=lambda i:-A[sel][i])])
    # envelope over time of h1: print 50ms rms in dB
    fh=f0; y=bp(x,fs,fh-60,fh+60,order=4); a=np.abs(hilbert(y))
    blk=2400; e=[20*np.log10(a[i:i+blk].mean()) for i in range(0,len(a)-blk,blk)]
    print('   h1 env per 50ms:',np.round(np.array(e)-np.median(e),1).tolist())
