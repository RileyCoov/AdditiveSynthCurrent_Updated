import pickle, numpy as np
for n in ['1985','HappyMono','take-me-out']:
    ev=pickle.load(open(f'an/{n}_events.pkl','rb'))
    for v,r in [('up_fifth',2**(7/12)),('down_fifth',2**(-7/12))]:
        E=ev[v]
        f=np.array([e['f']/r for e in E]); d=np.array([(e['t1']-e['t0'])*1000+85 for e in E]); rel=np.array([e['rel_to_frame'] for e in E]); ex=np.array([e['exc'] for e in E])
        prom=rel>-12
        print(f"{n:12s} {v:10s} all N={len(E)}: src f median {np.median(f):6.0f} Hz, frac src>4kHz {np.mean(f>4000):.2f}, dur median {np.median(d):4.0f} ms, rel-to-frame median {np.median(rel):5.1f} dB"
              f" || prominent (rel>-12 dB) N={prom.sum()}: src f median {np.median(f[prom]) if prom.any() else 0:6.0f}, frac>4k {np.mean(f[prom]>4000) if prom.any() else 0:.2f}, dur median {np.median(d[prom]) if prom.any() else 0:4.0f} ms, excess median {np.median(ex[prom]) if prom.any() else 0:4.1f}")
