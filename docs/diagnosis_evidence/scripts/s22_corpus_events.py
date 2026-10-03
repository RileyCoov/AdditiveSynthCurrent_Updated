from tf import *
import pickle
fs=48000
rows=[]
for n in NAMES:
    xo,_=load(n,'original'); f,t,Xo=spec(xo,fs); dur=len(xo)/fs
    res={}
    for v in VERS[1:]:
        y,_=load(n,v); _,_,Y=spec(y,fs)
        Xm=mapped_src(Xo,f,R[v])
        added,LY,LX,exc=added_tonal(Y,Xm,f)
        ev=events(added,LY,exc,f,t,min_frames=4)
        for e in ev:
            j0=np.searchsorted(t,e['t0']); j1=np.searchsorted(t,e['t1'])+1
            e['rel_to_frame']=e['peak_db']-float(LY[:,j0:j1].max())
        res[v]=ev
        # salient = within 20 dB of frame max and excess>=15
        sal=[e for e in ev if e['rel_to_frame']>-20 and e['exc']>=15]
        res[v+'_sal']=sal
    pickle.dump(res,open(f'an/{n}_events.pkl','wb'))
    # symmetry among salient
    ru,rd=R['up_fifth'],R['down_fifth']
    def symfrac(A,ra,B,rb):
        if not A: return np.nan
        h=0
        for a in A:
            fa=a['f']/ra
            if any(abs(fa-b['f']/rb)/fa<0.03 and a['t0']<=b['t1']+0.02 and b['t0']<=a['t1']+0.02 for b in B): h+=1
        return h/len(A)
    s_up=symfrac(res['up_fifth_sal'],ru,res['down_fifth'],rd)
    rows.append((n,dur,len(res['unity']),len(res['up_fifth'])/dur,len(res['down_fifth'])/dur,len(res['up_fifth_sal'])/dur,len(res['down_fifth_sal'])/dur,s_up))
    print(f"{n:30s} dur {dur:5.1f}s | events/s (>=4 hops, exc>10dB): unity {len(res['unity'])/dur:5.1f} up {len(res['up_fifth'])/dur:6.1f} dn {len(res['down_fifth'])/dur:6.1f} | salient/s (exc>=15, within 20 dB of frame max): up {len(res['up_fifth_sal'])/dur:5.2f} dn {len(res['down_fifth_sal'])/dur:5.2f} | salient-up with down twin: {s_up:.2f}",flush=True)
