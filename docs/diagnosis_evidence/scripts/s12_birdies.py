from tf import *
import json
n=sys.argv[1] if len(sys.argv)>1 else '1985'
fs=48000
xo,_=load(n,'original'); f,t,Xo=spec(xo,fs)
allev={}
for v in VERS[1:]:
    y,_=load(n,v); _,_,Y=spec(y,fs)
    Xm=mapped_src(Xo,f,R[v])
    added,LY,LX,exc=added_tonal(Y,Xm,f)
    ev=events(added,LY,exc,f,t,min_frames=4)   # >= 4 hops (~21 ms extent beyond window)
    # rank by loudness relative to frame's broadband level
    for e in ev:
        j0=np.searchsorted(t,e['t0']); j1=np.searchsorted(t,e['t1'])+1
        e['frame_max']=float(LY[:,j0:j1].max())
        e['rel_to_frame']=e['peak_db']-e['frame_max']
    ev.sort(key=lambda e:-(e['exc']+ e['rel_to_frame']*0.5 + 2*np.log2(e['nfr'])))
    allev[v]=ev
    per_s=added.sum()/ (len(xo)/fs)
    print(f"\n== {n} {v}: added-tonal TF cells/s {per_s:.1f}; events(>=4 hops) {len(ev)}; window center times")
    for e in ev[:25]:
        print(f"  {e['t0']:6.3f}-{e['t1']:6.3f}s  ({int(e['t0']*fs):7d}-{int(e['t1']*fs):7d})  f~{e['f']:7.1f} Hz [{e['fmin']:.0f}-{e['fmax']:.0f}]  hops={e['nfr']:3d}  peak {e['peak_db']:6.1f} dB  excess {e['exc']:5.1f} dB  rel-to-frame-max {e['rel_to_frame']:6.1f}")
import pickle; pickle.dump(allev,open(f'an/{n}_events.pkl','wb'))
