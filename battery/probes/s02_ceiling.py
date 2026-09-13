"""Stage 0.2 - oracle-fit ceiling probe.

Question: is the engine's remaining error an ESTIMATION problem (better analysis could
fix it) or a MODEL-STRUCTURE problem (the sinusoidal model cannot represent this signal
at all)?

Method: bypass the analysis stage entirely. For each frame, fit the best possible set of
K sinusoids directly to the input by matching pursuit on a 4x-oversampled frequency grid
(each atom solved in closed form against the running residual), then overlap-add exactly
the way the engine does. Frequencies are free per frame and phases are unconstrained, so
this is a strict UPPER BOUND on anything the engine's synthesis model class can achieve,
with or without a learned front end.

Sweeping K gives the shape of the ceiling:
  - SRR climbing steadily with K  -> capacity/estimation limited (more/better partials help)
  - SRR saturating low            -> STRUCTURAL: the model cannot represent this material
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from lib import read_wav, align, srr

N = 4096            # analysis/synthesis frame
HOP = 1024          # 75% overlap
ZP = 4              # frequency-grid oversampling
KS = [8, 16, 32, 64, 128, 256]
KMAX = max(KS)
DUR = 2.5           # seconds analysed per file


def frame_fit(x, w, kmax, checkpoints):
    """Matching pursuit of kmax real sinusoids on one frame.

    Returns {k: model_signal} at each checkpoint k.
    """
    n = np.arange(len(x))
    r = x.copy()
    model = np.zeros_like(x)
    out = {}
    nz = len(x) * ZP
    for k in range(1, kmax + 1):
        S = np.fft.rfft(w * r, nz)
        mag = np.abs(S)
        mag[0] = mag[-1] = 0.0
        i = int(np.argmax(mag))
        # parabolic refinement on the zero-padded magnitude grid
        if 0 < i < len(mag) - 1:
            a_, b_, c_ = mag[i - 1], mag[i], mag[i + 1]
            d = a_ - 2 * b_ + c_
            p = 0.5 * (a_ - c_) / d if d != 0 else 0.0
        else:
            p = 0.0
        f = (i + p) / nz                      # cycles/sample
        c = np.cos(2 * np.pi * f * n)
        s = np.sin(2 * np.pi * f * n)
        # exact weighted 2x2 least squares of the residual onto {cos, sin}
        wc, ws = w * c, w * s
        g11, g12, g22 = c @ wc, c @ ws, s @ ws
        r1, r2 = r @ wc, r @ ws
        det = g11 * g22 - g12 * g12
        if abs(det) < 1e-12:
            break
        A = (r1 * g22 - r2 * g12) / det
        B = (r2 * g11 - r1 * g12) / det
        comp = A * c + B * s
        r -= comp
        model += comp
        if k in checkpoints:
            out[k] = model.copy()
    for k in checkpoints:                      # kmax reached early
        if k not in out:
            out[k] = model.copy()
    return out


def oracle(x, kmax=KMAX, checkpoints=tuple(KS)):
    w = np.hanning(N + 1)[:N]   # periodic Hann -> exact COLA at hop N/4
    nfr = max(1, 1 + (len(x) - N) // HOP)
    acc = {k: np.zeros(len(x)) for k in checkpoints}
    wsum = np.zeros(len(x))
    for t in range(nfr):
        a, b = t * HOP, t * HOP + N
        if b > len(x):
            break
        seg = x[a:b]
        if np.abs(seg).max() < 1e-7:
            wsum[a:b] += w
            continue
        fits = frame_fit(seg, w, kmax, set(checkpoints))
        for k in checkpoints:
            acc[k][a:b] += fits[k] * w
        wsum[a:b] += w
    safe = wsum.copy(); safe[safe < 1e-6] = 1e-6
    # only the interior is fully covered by the Hann overlap-add; the ramp-in /
    # ramp-out regions divide by a small window sum and are not part of the model's
    # capability. Evaluate on the plateau only.
    valid = wsum > 0.98 * wsum.max()
    return {k: acc[k] / safe for k in checkpoints}, valid


def run():
  IN = os.path.expanduser('~/Desktop/comparingOut/inputData')
  ENG = '/tmp/stage0/render'
  meter = {r['name']: r for r in json.load(open('/tmp/stage0/s01_residual_meter.json'))}

  FILES = ['300hzSine_5sec', '300hzSaw_5sec', 'FairlightAdditiveWaveTA1_C3',
         'PianoSampleMono', 'SaintSaensMonoChord', 'Female_Sung_Line_3_48',
         'choir-burst_C_major_Mono_48k', 'DrumLoopShort', '48kCymbal',
         'HappyMono', 'take-me-out']

  rows = []
  for name in FILES:
    t0 = time.time()
    x, sr = read_wav(f'{IN}/{name}.wav')
    y, _ = read_wav(f'{ENG}/{name}_orig.wav')
    x, y = align(x, y)
    # pick the loudest DUR-second window so we analyse signal, not silence
    L = int(DUR * sr)
    if len(x) > L:
        env = np.convolve(x ** 2, np.ones(sr // 10) / (sr // 10), 'same')
        c = int(np.argmax(np.convolve(env, np.ones(L) / L, 'same')))
        a = max(0, min(len(x) - L, c - L // 2))
        x, y = x[a:a + L], y[a:a + L]
    recs, valid = oracle(x)
    # score engine and oracle on the identical, fully-covered sample set
    xv, yv = x[valid], y[valid]
    eng_srr, _ = srr(xv, yv)
    curve = {k: srr(xv, recs[k][valid])[0] for k in KS}
    rows.append(dict(name=name, cls=meter[name]['cls'], sr=sr,
                     engine=eng_srr, curve={str(k): curve[k] for k in KS}))
    print(f"{name:30s} engine={eng_srr:5.1f} dB | oracle " +
          ' '.join(f'K{k}={curve[k]:5.1f}' for k in KS) +
          f"  [{time.time()-t0:.0f}s]", flush=True)

  json.dump(rows, open('/tmp/stage0/s02_ceiling.json', 'w'), indent=1)


if __name__ == '__main__':
    run()
