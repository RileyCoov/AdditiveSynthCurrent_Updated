from common import *
def envlag(a,b,maxlag=0.06,hop=48):
    n=min(len(a),len(b))//hop
    ea=np.log(a[:n*hop].reshape(n,hop).mean(1)+1e-12); eb=np.log(b[:n*hop].reshape(n,hop).mean(1)+1e-12)
    ea=np.diff(ea); eb=np.diff(eb)   # onset-sensitive
    L=int(maxlag*48000/hop)
    c=[np.dot(ea[max(0,-l):len(ea)-max(0,l)], eb[max(0,l):len(eb)-max(0,-l)]) for l in range(-L,L+1)]
    c=np.array(c); i=np.argmax(c)
    # parabolic refine
    if 0<i<len(c)-1:
        d=0.5*(c[i-1]-c[i+1])/(c[i-1]-2*c[i]+c[i+1])
    else: d=0
    return (i-L+d)*hop/48
bands=[(40,200),(200,800),(800,3000),(3000,10000)]
for n in ['1985','DrumLoopShort','HappyMono','take-me-out','PianoSampleMono','48kCymbal']:
    X={v:load(n,v)[0] for v in VERS}; fs=48000
    line=f"{n:16s}"
    for lo,hi in bands:
        xo=bp(X['original'],fs,lo,hi)**2
        ls=[]
        for v in VERS[1:]:
            r=R[v]; y=bp(X[v],fs,lo*r,hi*r)**2
            ls.append(envlag(xo,y))
        line+=f" | {lo}-{hi}: u{ls[0]:+5.1f} up{ls[1]:+5.1f} dn{ls[2]:+5.1f}"
    print(line)
print("(ms, positive = output later than original; onset-flux xcorr, 1 ms hop)")
