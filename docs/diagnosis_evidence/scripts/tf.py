from common import *
from scipy.ndimage import median_filter, maximum_filter, label, find_objects
NFFT=4096; HOP=256
def spec(x,fs):
    f,t,Z=stft(x,fs,'hann',nperseg=NFFT,noverlap=NFFT-HOP,boundary=None,padded=False)
    return f,t,np.abs(Z)*2/ (NFFT/2)   # ~ amplitude
def mapped_src(Xmag, f, r):
    # expected output magnitude at output bin f_out = |X(f_out/r)| (sinusoid amplitude preserved)
    out=np.empty_like(Xmag)
    fsrc=f/r
    for j in range(Xmag.shape[1]):
        out[:,j]=np.interp(fsrc,f,Xmag[:,j],right=0)
    return out
def added_tonal(Y, Xm, f, k_med=15, thr_peak=10, thr_excess=10, tol_bins=2):
    LY=20*np.log10(Y+1e-9); LX=20*np.log10(Xm+1e-9)
    medY=median_filter(LY,size=(2*k_med+1,1),mode='nearest')
    pkY=(LY>=maximum_filter(LY,size=(3,1)))&(LY-medY>thr_peak)
    medX=median_filter(LX,size=(2*k_med+1,1),mode='nearest')
    pkX=(LX>=maximum_filter(LX,size=(3,1)))&(LX-medX>thr_peak-4)  # lenient source peaks
    pkXd=maximum_filter(pkX.astype(np.uint8),size=(2*tol_bins+1,3))>0
    # source local max level near bin
    LXmax=maximum_filter(LX,size=(2*tol_bins+1,3))
    excess=LY-LXmax
    added=pkY&(~pkXd)&(excess>thr_excess)&(LY>LY.max()-80)
    return added, LY, LX, excess
def events(added, LY, excess, f, t, min_frames=2):
    # connect across time with ±1 bin freq tolerance
    from scipy.ndimage import binary_dilation
    st=np.ones((3,3),bool)
    lab,nl=label(binary_dilation(added,np.array([[0,1,0],[0,1,0],[0,1,0]],bool)),structure=st)
    lab=lab*added
    ev=[]
    for i,sl in enumerate(find_objects(lab),1):
        if sl is None: continue
        m=(lab[sl]==i)
        nfr=m.any(0).sum()
        if nfr<min_frames: continue
        fi=np.where(m)[0]+sl[0].start; ti=np.where(m)[1]+sl[1].start
        lev=LY[fi,ti]; exc=excess[fi,ti]
        ev.append(dict(t0=t[ti.min()],t1=t[ti.max()],f=np.median(f[fi]),fmin=f[fi].min(),fmax=f[fi].max(),
                       nfr=nfr,peak_db=lev.max(),exc=exc.mean(),score=np.sum(10**(exc/20))*1.0, rel=lev.max()))
    return ev
