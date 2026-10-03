from tf import *
import pickle
n=sys.argv[1] if len(sys.argv)>1 else '1985'
ev=pickle.load(open(f'an/{n}_events.pkl','rb'))
up=ev['up_fifth']; dn=ev['down_fifth']
ru=R['up_fifth']; rd=R['down_fifth']
def match(A,ra,B,rb,tshift=0.0,ftol=0.03,ttol=0.02):
    hits=0; pairs=[]
    for a in A:
        fa=a['f']/ra
        for b in B:
            fb=b['f']/rb
            if abs(fa-fb)/fa<ftol and (a['t0']+tshift<=b['t1']+ttol) and (b['t0']<=a['t1']+tshift+ttol):
                hits+=1; pairs.append((a,b)); break
    return hits,pairs
# restrict to the more salient events: top N by excess*log(nfr)
def top(E,N): 
    return sorted(E,key=lambda e:-(e['exc']+2*np.log2(e['nfr'])+0.3*e['rel_to_frame']))[:N]
for N in [25,50,100,200]:
    A=top(up,N); B=top(dn,N)
    h,pairs=match(A,ru,dn,rd)
    # null: shift up event times by random offsets (circular) 
    rng=np.random.default_rng(0); T=max(e['t1'] for e in up)
    null=[]
    for k in range(200):
        s=rng.uniform(0.3,T-0.3)
        A2=[dict(e, t0=(e['t0']+s)%T, t1=(e['t0']+s)%T+(e['t1']-e['t0'])) for e in A]
        null.append(match(A2,ru,dn,rd)[0])
    print(f"top{N} up events: {h}/{N} have a down event at same SOURCE freq (±3%) & overlapping time; null (time-shuffled) mean {np.mean(null):.1f}, 95th pct {np.percentile(null,95):.0f}")
A=top(up,40)
h,pairs=match(A,ru,dn,rd)
print("\nSymmetric pairs among top-40 up events (source-coord freq):")
for a,b in sorted(pairs,key=lambda p:p[0]['t0']):
    print(f"  {a['t0']:.3f}-{a['t1']:.3f}s up {a['f']:.0f}Hz (src {a['f']/ru:.0f}) | down {b['t0']:.3f}-{b['t1']:.3f}s {b['f']:.0f}Hz (src {b['f']/rd:.0f})  exc up {a['exc']:.1f} dn {b['exc']:.1f}")
