"""Stage 0.2c - decompose the engine-to-ceiling gap into FREQUENCY error vs AMPLITUDE/PHASE error.

Well-posed and convex, unlike direct optimisation of a phase-continuous model.

For each frame we compare three reconstructions of the same audio:
  ENGINE   - what the engine actually renders.
  ENG-FREQ - the engine's OWN track frequencies for that frame, but with amplitude and
             phase solved in closed form by weighted least squares (the best any
             amplitude/phase estimator could do given those frequencies).
  ORACLE   - frequencies also chosen freely by matching pursuit (Ceiling A).

Reading:
  ENGINE -> ENG-FREQ  = headroom recoverable WITHOUT changing frequency estimation
                        (i.e. amplitude + phase estimation quality)
  ENG-FREQ -> ORACLE  = headroom that requires better FREQUENCY estimation / selection
"""
import sys, os, json, subprocess, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from lib import read_wav, align, srr
from s02_ceiling import N, HOP, oracle

BIN = '/Users/Riley/Downloads/AdditiveSynthCurrent-LongWindowFreqUsedForShortWindows/build/additive'
IN = os.path.expanduser('~/Desktop/comparingOut/inputData')
ENG = '/tmp/stage0/render'
DUR = 2.5
KMAX = 256


def dump_frames(name):
    """Run the engine with PFREQ_DEBUG and parse per-frame (start, freqs, amps)."""
    env = dict(os.environ, PFREQ_DEBUG='1')
    p = subprocess.run([BIN, f'{IN}/{name}.wav', '/tmp/x.wav', '2048', '0'],
                       capture_output=True, text=True, env=env)
    frames = []
    cur = None
    for ln in p.stderr.splitlines():
        if ln.startswith('PFRAME'):
            _, fi, st, fs, npk = ln.split()
            cur = dict(start=int(st), size=int(fs), f=[], a=[])
            frames.append(cur)
        elif ln.startswith('PF ') and cur is not None:
            _, f, a = ln.split()
            cur['f'].append(float(f)); cur['a'].append(float(a))
    for fr in frames:
        fr['f'] = np.array(fr['f']); fr['a'] = np.array(fr['a'])
    return frames


def fit_fixed_freqs(x, sr, frames, lo, hi):
    """Optimal weighted-LS amplitude/phase at the engine's own frequencies, OLA'd.

    Uses the engine's OWN frame geometry (start + frame_size, 2048 long / 256 short)
    so the comparison against the engine's render is apples-to-apples.
    """
    wcache = {}
    acc = np.zeros(len(x)); wsum = np.zeros(len(x))
    for fr in frames:
        M = fr['size']
        a0 = fr['start'] - lo
        if a0 < 0 or a0 + M > len(x):
            continue
        if M not in wcache:
            wcache[M] = (np.hanning(M + 1)[:M], np.arange(M))
        w, n = wcache[M]
        seg = x[a0:a0 + M]
        keep = (fr['f'] > 0) & (fr['f'] < sr / 2 - 1)
        f, amp = fr['f'][keep], fr['a'][keep]
        if len(f) == 0 or np.abs(seg).max() < 1e-8:
            wsum[a0:a0 + M] += w; continue
        # a frame of M samples supports at most ~M/2 independent sinusoids; cap by
        # amplitude so the LS stays well-posed and cannot simply span the frame.
        cap = min(KMAX, max(1, M // 4))
        if len(f) > cap:
            f = f[np.argsort(-amp)[:cap]]
        ang = 2 * np.pi * np.outer(n, f / sr)
        A = np.concatenate([np.cos(ang), np.sin(ang)], axis=1)
        sw = np.sqrt(w)[:, None]
        coef, *_ = np.linalg.lstsq(A * sw, seg * np.sqrt(w), rcond=1e-8)
        acc[a0:a0 + M] += (A @ coef) * w
        wsum[a0:a0 + M] += w
    safe = wsum.copy(); safe[safe < 1e-6] = 1e-6
    return acc / safe, wsum > 0.98 * np.median(wsum[wsum > 0])


if __name__ == '__main__':
    rows = []
    for name in ['FairlightAdditiveWaveTA1_C3', 'PianoSampleMono', 'SaintSaensMonoChord',
                 'Female_Sung_Line_3_48', 'choir-burst_C_major_Mono_48k',
                 'DrumLoopShort', '48kCymbal', 'take-me-out']:
        t0 = time.time()
        x, sr = read_wav(f'{IN}/{name}.wav')
        # residual-off render: the engine's TONAL model alone, which is what the
        # fixed-frequency LS fit is also modelling. Comparing against the shipped
        # render would credit the engine with its noise fill and hide the point.
        y, _ = read_wav(f'/tmp/stage0/render_nr/{name}_orig.wav')
        xa, ya = align(x, y)
        L = int(DUR * sr)
        if len(xa) > L:
            env = np.convolve(xa ** 2, np.ones(sr // 10) / (sr // 10), 'same')
            c = int(np.argmax(np.convolve(env, np.ones(L) / L, 'same')))
            lo = max(0, min(len(xa) - L, c - L // 2))
        else:
            lo = 0; L = len(xa)
        hi = lo + L
        xs, ys = xa[lo:hi], ya[lo:hi]

        frames = dump_frames(name)
        fixed, vfix = fit_fixed_freqs(xs, sr, frames, lo, hi)
        recA, vA = oracle(xs, kmax=KMAX, checkpoints=(KMAX,))
        v = vfix & vA
        e = srr(xs[v], ys[v])[0]
        ef = srr(xs[v], fixed[v])[0]
        oa = srr(xs[v], recA[KMAX][v])[0]
        rows.append(dict(name=name, engine=e, eng_freq=ef, oracle=oa))
        print(f"{name:30s} ENGINE={e:5.1f}  ENG-FREQ={ef:5.1f}  ORACLE={oa:5.1f}   "
              f"| amp/phase headroom={ef-e:5.1f} dB, freq headroom={oa-ef:5.1f} dB "
              f"[{time.time()-t0:.0f}s]", flush=True)
    json.dump(rows, open('/tmp/stage0/s04_freq_vs_amp.json', 'w'), indent=1)
