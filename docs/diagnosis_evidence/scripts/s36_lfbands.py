from common import *
fs=48000
edges=[40,60,80,100,120,140,160,180,200,220,240,270,300,350,400]
for n in ['HappyMono','take-me-out','DrumLoopShort','1985','Female_Sung_Line_3_48']:
    X={v:load(n,v)[0] for v in VERS}
    line=f"{n:24s}"
    for v in ['up_fifth','down_fifth']:
        r=R[v]; d=[]
        for lo,hi in zip(edges[:-1],edges[1:]):
            po=np.mean(bp(X['original'],fs,lo,hi,order=6)**2); pv=np.mean(bp(X[v],fs,lo*r,hi*r,order=6)**2)
            d.append(10*np.log10(pv/po))
        line+=f"\n   {v:10s} "+" ".join(f"{lo}:{x:+4.1f}" for lo,x in zip(edges[:-1],d))
    print(line)
