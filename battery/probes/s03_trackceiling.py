"""Stage 0.2b - track-constrained ceiling (Ceiling B).

Ceiling A (s02) refits every frame independently with free frequencies and phases.
A real additive engine cannot do that: it must carry CONTINUOUS tracks with
continuous phase. This probe asks how much of the A-vs-engine gap is simply the
price of trackability.

Model (exactly the engine's MQ synthesis class):
    y[n] = sum_k a_k(n) * cos( phi0_k + 2*pi/sr * cumsum(f_k)[n] )
with per-frame breakpoints for f and a, linearly interpolated to sample rate, and
phase obtained by integrating frequency - so phase is continuous by construction
and every track is a single coherent partial for the whole excerpt.

Optimised by Adam on a waveform-domain loss. Initialised from the per-frame
matching-pursuit fit, so the reported value is a LOWER bound on Ceiling B; if even
this lower bound sits far above the engine, trackability is not the explanation.
"""
import sys, os, json, time
import numpy as np, torch
sys.path.insert(0, os.path.dirname(__file__))
from lib import read_wav, align, srr
from s02_ceiling import N, HOP, ZP

DUR = 1.0
K = 128
ITERS = 1500


def mp_select(x, w, k):
    """Matching pursuit: return (freqs cycles/sample, amps, phases) for one frame."""
    n = np.arange(len(x)); r = x.copy(); nz = len(x) * ZP
    fs, As, Bs = [], [], []
    for _ in range(k):
        S = np.abs(np.fft.rfft(w * r, nz)); S[0] = S[-1] = 0
        i = int(np.argmax(S))
        if 0 < i < len(S) - 1:
            a_, b_, c_ = S[i-1], S[i], S[i+1]; d = a_ - 2*b_ + c_
            p = 0.5 * (a_ - c_) / d if d != 0 else 0.0
        else:
            p = 0.0
        f = (i + p) / nz
        c, s = np.cos(2*np.pi*f*n), np.sin(2*np.pi*f*n)
        wc, ws = w*c, w*s
        g11, g12, g22 = c@wc, c@ws, s@ws
        r1, r2 = r@wc, r@ws
        det = g11*g22 - g12*g12
        if abs(det) < 1e-12: break
        A = (r1*g22 - r2*g12)/det; B = (r2*g11 - r1*g12)/det
        r -= A*c + B*s
        fs.append(f); As.append(A); Bs.append(B)
    fs, As, Bs = np.array(fs), np.array(As), np.array(Bs)
    amp = np.hypot(As, Bs)
    while len(fs) < k:   # pad if MP terminated early
        fs = np.append(fs, 0.25); amp = np.append(amp, 1e-8)
    o = np.argsort(fs)   # sort by frequency -> natural track continuity
    return fs[o], amp[o]


def track_ceiling(x, sr, k=K, iters=ITERS, dev='cpu'):
    L = len(x)
    w = np.hanning(N + 1)[:N]
    starts = list(range(0, max(1, L - N), HOP))
    T = len(starts)
    F0 = np.zeros((T, k)); A0 = np.zeros((T, k))
    for t, a in enumerate(starts):
        seg = x[a:a+N]
        if len(seg) < N or np.abs(seg).max() < 1e-7:
            F0[t] = np.linspace(0.001, 0.4, k); A0[t] = 1e-8; continue
        F0[t], A0[t] = mp_select(seg, w, k)

    # breakpoint times = frame centres, in samples
    tc = np.array([a + N // 2 for a in starts], dtype=np.float64)
    n = np.arange(L, dtype=np.float64)

    dev = torch.device(dev)
    f = torch.tensor(np.clip(F0, 1e-5, 0.49), dtype=torch.float32, device=dev, requires_grad=True)
    la = torch.tensor(np.log(np.clip(A0, 1e-8, None)), dtype=torch.float32, device=dev, requires_grad=True)
    ph = torch.tensor(np.random.RandomState(0).uniform(0, 2*np.pi, k), dtype=torch.float32,
                      device=dev, requires_grad=True)
    xt = torch.tensor(x, dtype=torch.float32, device=dev)
    nt = torch.tensor(n, dtype=torch.float32, device=dev)
    tct = torch.tensor(tc, dtype=torch.float32, device=dev)

    def interp(vals):                      # [T,k] breakpoints -> [L,k] per-sample
        idx = torch.clamp(torch.searchsorted(tct, nt.contiguous(), right=True) - 1, 0, len(tct) - 2)
        t0, t1 = tct[idx], tct[idx + 1]
        u = ((nt - t0) / (t1 - t0)).clamp(0, 1).unsqueeze(1)
        return vals[idx] * (1 - u) + vals[idx + 1] * u

    opt = torch.optim.Adam([f, la, ph], lr=3e-3)
    best = None
    for it in range(iters):
        opt.zero_grad()
        fn = interp(f).clamp(1e-6, 0.499)
        an = torch.exp(interp(la))
        phase = ph.unsqueeze(0) + 2 * np.pi * torch.cumsum(fn, dim=0)
        y = (an * torch.cos(phase)).sum(dim=1)
        loss = ((y - xt) ** 2).mean()
        loss.backward(); opt.step()
        if it % 100 == 0 or it == iters - 1:
            v = float(loss)
            if best is None or v < best[0]:
                best = (v, y.detach().cpu().numpy().copy())
    return best[1]


if __name__ == '__main__':
    IN = os.path.expanduser('~/Desktop/comparingOut/inputData')
    ENG = '/tmp/stage0/render'
    from s02_ceiling import oracle
    out = []
    for name in ['choir-burst_C_major_Mono_48k', 'take-me-out', '48kCymbal',
                 'SaintSaensMonoChord', 'Female_Sung_Line_3_48']:
        x, sr = read_wav(f'{IN}/{name}.wav')
        y, _ = read_wav(f'{ENG}/{name}_orig.wav')
        x, y = align(x, y)
        L = int(DUR * sr)
        env = np.convolve(x**2, np.ones(sr//10)/(sr//10), 'same')
        c = int(np.argmax(np.convolve(env, np.ones(L)/L, 'same')))
        a = max(0, min(len(x)-L, c - L//2))
        x, y = x[a:a+L], y[a:a+L]
        t0 = time.time()
        recA, valid = oracle(x, kmax=K, checkpoints=(K,))
        yb = track_ceiling(x, sr)
        eng = srr(x[valid], y[valid])[0]
        ca = srr(x[valid], recA[K][valid])[0]
        cb = srr(x[valid], yb[valid])[0]
        out.append(dict(name=name, engine=eng, ceilA=ca, ceilB=cb))
        print(f"{name:30s} engine={eng:5.1f}  CeilB(tracked,K{K})={cb:5.1f}  "
              f"CeilA(free,K{K})={ca:5.1f}   [{time.time()-t0:.0f}s]", flush=True)
    json.dump(out, open('/tmp/stage0/s03_trackceiling.json', 'w'), indent=1)
