from common import *
LSB=1/32768
def transient_segments(n):
    xo,fs=load(n,'original'); xu,_=load(n,'unity'); e=xu-xo
    ex=np.abs(e)<=1.01*LSB
    d=np.diff(np.concatenate([[0],ex.astype(int),[0]])); s=np.where(d==1)[0]; t=np.where(d==-1)[0]
    keep=(t-s)>=1000; s,t=s[keep],t[keep]
    segs=[]
    for a,b in zip(s,t):
        S=(a//1024)*1024; E=((b+1023)//1024)*1024
        segs.append((S,E))
    return segs
