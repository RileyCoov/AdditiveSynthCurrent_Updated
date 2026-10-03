from segs import *
import pickle
fs=48000; n='1985'
ev=pickle.load(open(f'an/{n}_events.pkl','rb'))
segs=transient_segments(n)
print("1985 transient segments (s):", [(round(S/fs,3),round(E/fs,3)) for S,E in segs])
regions=[(1.20,1.31,'S'),(2.88,3.00,'S'),(4.08,4.29,'S'),(4.46,4.64,'S'),(6.31,6.38,'S'),(6.59,6.88,'S'),
         (2.10,2.14,'W'),(2.50,2.54,'W'),(5.46,5.50,'W'),(5.79,5.83,'W'),(7.37,7.41,'W'),(8.63,8.67,'W'),(8.91,8.95,'W')]
ru,rd=R['up_fifth'],R['down_fifth']
def top(E,k=150): return sorted(E,key=lambda e:-(e['exc']+2*np.log2(e['nfr'])+0.3*e['rel_to_frame']))[:k]
U=top(ev['up_fifth']); D=top(ev['down_fifth'])
for a,b,kind in regions:
    segov=[(round(S/fs,3),round(E/fs,3)) for S,E in segs if S/fs<b+0.02 and E/fs>a-0.02]
    def ovl(E,r):
        out=[]
        for e in E:
            # artifact audible span ≈ [t0-0.043, t1+0.043] (window half-length)
            if e['t0']-0.043<b and e['t1']+0.043>a:
                out.append(f"{e['t0']:.2f}-{e['t1']:.2f}s {e['f']:.0f}Hz(src {e['f']/r:.0f}) exc{e['exc']:.0f} rel{e['rel_to_frame']:.0f}")
        return out
    u=ovl(U,ru); d=ovl(D,rd)
    print(f"\n[{kind}] {a:.2f}-{b:.2f}s  ({int(a*fs)}-{int(b*fs)})  transient segs overlapping: {segov}")
    print("    UP  :", u[:6]); print("    DOWN:", d[:6])
