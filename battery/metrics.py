"""Standing DSP quality metrics for the additive resynthesis regression battery.

Each function turns a perceptual complaint into a reproducible number, so that
"double voice" / "shape ripple" / "volume wobble" stop being a subjective
coin-flip. Pure functions over mono float arrays; no I/O here (see run.py).

Metrics
-------
waveform_correlation : shape fidelity vs a reference (the saw "ripple" metric)
cepstral_excess      : doubling ADDED by resynthesis vs the input (the Female
                       "double voice" metric). Reference-relative so the source's
                       own pitch rahmonics cancel -- an absolute in-band cepstral
                       peak was tested and could NOT tell a doubled voice from a
                       clean one (its rahmonics fill the same 5-40 ms band), so
                       the reference-relative form is the one that survives.
envelope_p2p         : amplitude-envelope peak-to-peak, full and band-limited
                       (the steady-tone "volume wobble" metric, and the proxy for
                       doubling BEATING under pitch shift where no clean
                       same-pitch reference exists)
residual_srr         : how much of the input the model actually captures (the
                       residual meter). Phase-sensitive, defined on every class.
residual_flatness    : is what's LEFT OVER noise-like or tonal -- i.e. did the
                       model miss noise (expected) or miss partials (a defect)?

Why residual_srr exists
-----------------------
Before it, the only phase-sensitive gate in the battery was waveform_correlation,
declared for the two saw entries. Everything else was guarded by magnitude- and
envelope-domain metrics only. That is a blind spot with a measured cost: the
largest fidelity fix in this engine's history (5ae358e, fractional-bin phase
pickup) moved 300hzSaw waveform correlation 0.504 -> 0.992 while moving magnitude
log-spectral distance by only 7.32 -> 7.10 dB. A magnitude-domain battery would
have scored that fix as noise. residual_srr is phase-sensitive by construction
(it is an error energy on the time-domain difference) and is cheap on every class.
"""
from __future__ import annotations
import numpy as np
from scipy.signal import hilbert, butter, sosfiltfilt


def _mono(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if x.ndim > 1:
        x = x.mean(axis=1)
    return x


def waveform_correlation(rendered: np.ndarray, reference: np.ndarray,
                         max_lag: int = 512) -> float:
    """Max normalized cross-correlation of two waveforms over +/- max_lag samples.

    Removes any constant synthesis delay by searching integer lags, so what is
    measured is *shape* fidelity, not alignment. 1.0 = identical shape; the saw
    target is >= 0.999 (r7 baseline ~0.996). Both signals are trimmed to equal
    length and zero-mean/unit-norm before correlating.
    """
    a = _mono(rendered)
    b = _mono(reference)
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    a = a - a.mean()
    b = b - b.mean()
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    a, b = a / na, b / nb
    best = -1.0
    for lag in range(-max_lag, max_lag + 1):
        if lag < 0:
            v = float(np.dot(a[:n + lag], b[-lag:]))
        elif lag > 0:
            v = float(np.dot(a[lag:], b[:n - lag]))
        else:
            v = float(np.dot(a, b))
        if v > best:
            best = v
    return best


def _cepstrum(x: np.ndarray) -> np.ndarray:
    y = _mono(x)
    y = y - y.mean()
    w = np.hanning(len(y))
    spec = np.fft.rfft(y * w)
    return np.fft.irfft(np.log(np.abs(spec) + 1e-12))


def cepstral_excess(rendered: np.ndarray, reference: np.ndarray, sr: int,
                    qmin_ms: float = 5.0, qmax_ms: float = 40.0) -> float:
    """Doubling ADDED by resynthesis, vs the original input (unity condition).

    Both signals share the same intrinsic pitch rahmonics, so subtracting the
    reference cepstrum from the rendered cepstrum cancels the source's own
    periodicity and leaves the echo/comb structure the OLA synthesis introduced.
    Returns the max positive in-band excess, normalized by the reference band's
    median (so it is a relative "extra doubling" score). ~0 = rendered adds no
    doubling beyond the source; large = a resynthesis-introduced double voice.
    Requires rendered and reference at the SAME pitch (unity).
    """
    cr = _cepstrum(rendered)
    cf = _cepstrum(reference)
    n = min(len(cr), len(cf))
    cr, cf = cr[:n], cf[:n]
    qmin = int(qmin_ms * 1e-3 * sr)
    qmax = min(int(qmax_ms * 1e-3 * sr), n // 2)
    if qmax <= qmin:
        return 0.0
    excess = np.abs(cr[qmin:qmax]) - np.abs(cf[qmin:qmax])
    floor = np.median(np.abs(cf[qmin:qmax])) + 1e-12
    return float(max(0.0, excess.max()) / floor)


def trajectory_jitter(x: np.ndarray, sr: int, fmin: float = 100.0, fmax: float = 4000.0,
                      n_partials: int = 10, nperseg: int = 2048, hop: int = 256,
                      vib_hz: float = 15.0) -> float:
    """Per-partial amplitude jitter of the rendered output, in dB RMS.

    Targets the "phasey/warbly, not-one-person" voice artifact, which the session
    localized to per-frame analysis MEASUREMENT JITTER in each partial (frequency
    is de-jittered in the engine, amplitude is EMA'd, phase is re-derived raw each
    frame -> the jitter surfaces as fast amplitude/phase fluctuation on each
    harmonic). Tracks the strongest partials across an STFT, high-passes each
    partial's log-amplitude trajectory above the vibrato rate (so smooth vibrato
    and envelope are removed, jitter is kept), and returns the energy-weighted RMS.

    ~0 dB on clean steady/vibrato tones; rises with reconstruction jitter. Use it
    as a within-file gate: it must DROP on voice after the Phase-B phase-smoothing
    fix while guard rails (steady, Fairlight, drums) hold.

    NOTE: the companion frequency-jitter estimate (parabolic peak tracking) was
    tried and dropped — amplitude fluctuation makes the tracked bin jump, giving
    unstable values. Amplitude jitter is the robust proxy; revisit a phase-derived
    frequency-jitter term only if real-corpus validation shows amplitude is blind.
    """
    from scipy.signal import stft, butter, sosfiltfilt
    y = _mono(x)
    f, _, Z = stft(y, fs=sr, nperseg=nperseg, noverlap=nperseg - hop, boundary=None)
    mag = np.abs(Z)
    nb, nf = mag.shape
    if nf < 32:
        return 0.0
    g = max(8, nf // 20)
    sl = slice(g, nf - g)  # drop onset/offset edge frames
    avg = mag[:, sl].mean(1)
    cand = np.where((f >= fmin) & (f <= fmax))[0]
    order = cand[np.argsort(avg[cand])[::-1]]
    peaks = []
    for b in order:
        if avg[b] < 1e-4 * avg.max():
            break
        if all(abs(b - p) >= 3 for p in peaks):
            peaks.append(b)
        if len(peaks) >= n_partials:
            break
    fr = sr / hop
    sos = butter(2, vib_hz / (fr / 2), "high", output="sos")
    jit, wts = [], []
    for b in peaks:
        atr = np.zeros(nf)
        cur = b
        for n in range(nf):
            lo, hi = max(1, cur - 2), min(nb - 2, cur + 2)
            k = lo + int(np.argmax(mag[lo:hi + 1, n]))
            cur = k
            atr[n] = mag[k, n]
        atr = atr[sl]
        am = atr.mean()
        if am <= 1e-9:
            continue
        db = 20.0 * np.log10(np.maximum(atr, 1e-4 * atr.max()) / am)
        jit.append(np.sqrt(np.mean(sosfiltfilt(sos, db) ** 2)))
        wts.append(am)
    if not wts:
        return 0.0
    wts = np.array(wts) / np.sum(wts)
    return float(np.dot(wts, jit))


def _refine_f0(x: np.ndarray, sr: int, lo: float = 80.0, hi: float = 1200.0) -> float:
    n = min(1 << 19, len(x))
    X = np.abs(np.fft.rfft(x[:n] * np.hanning(n), 1 << 19))
    fr = np.fft.rfftfreq(1 << 19, 1.0 / sr)
    m = (fr > lo) & (fr < hi)
    if not m.any():
        return 0.0
    base = int(np.where(m)[0][0] + np.argmax(X[m]))
    if base <= 0 or base >= len(X) - 1:
        return fr[base]
    a, b, c = X[base - 1], X[base], X[base + 1]
    d = a - 2 * b + c
    p = 0.5 * (a - c) / d if d != 0 else 0.0
    return float((base + p) * sr / (1 << 19))


def waveform_shape_consistency(x: np.ndarray, sr: int, nseg: int = 400,
                               periods: int = 8) -> float:
    """Does the waveform keep ONE shape over the whole file?

    Riley's complaint about the sawtooths, in his words: "there isn't a consistent
    waveform shape like in the original, which makes the audio sound like the
    amplitude gets changed pretty often." A harmonic stack whose partials drift in
    RELATIVE phase keeps every partial's own amplitude flat while the summed
    waveform morphs -- and a morphing crest factor is heard as level movement.
    Neither envelope_p2p (dominated by inter-harmonic beating on a saw) nor
    trajectory_jitter (per-partial, so blind to relative phase) can see it.

    Method: estimate f0, cut the signal into period-locked blocks using fractional
    interpolation so a non-integer period stays aligned, and correlate each block
    against the median block. 1.0 = one stable shape throughout.

    Measured Sep 2026 (input -> engine): at UNITY the engine matches the source
    (300hzSaw 0.9989 vs 0.9990), but under pitch shift it collapses --
    300hzSaw up5 0.783, down5 0.833; 440sawtooth up5 0.914, down5 0.923 -- because
    shift-mode synthesis propagates each track's phase independently instead of
    deriving it per frame. This is the gate for the harmonic-phase-coherence round.
    """
    y = _mono(x)
    f0 = _refine_f0(y, sr)
    if f0 <= 0:
        return 0.0
    P = sr / f0
    L = int(P * periods)
    if L < 16 or len(y) < L + 4:
        return 0.0
    grid = np.arange(L)
    starts = np.linspace(0, len(y) - L - 2, nseg)
    segs = []
    for s in starts:
        s2 = round(s / P) * P               # snap to an exact period multiple
        idx = s2 + grid
        i0 = np.clip(idx.astype(int), 0, len(y) - 2)
        fr = idx - i0
        segs.append(y[i0] * (1 - fr) + y[i0 + 1] * fr)
    segs = np.array(segs)
    e = np.linalg.norm(segs, axis=1)
    med = np.median(e[e > 0]) if (e > 0).any() else 0.0
    segs = segs[e > 0.3 * med] if med > 0 else segs
    if len(segs) < 5:
        return 0.0
    ref = np.median(segs, axis=0)
    rn = np.linalg.norm(ref)
    if rn <= 0:
        return 0.0
    c = (segs @ ref) / (np.linalg.norm(segs, axis=1) * rn + 1e-20)
    return float(np.median(c))


def envelope_p2p_dev(rendered: np.ndarray, reference: np.ndarray, sr: int,
                     band=(0.3, 1.0)):
    """Envelope peak-to-peak DEVIATION from the reference: |p2p(out) - p2p(in)|.

    envelope_p2p and trajectory_jitter are absolute measurements of the OUTPUT, so
    "less variation" always scores better -- even when the input genuinely has that
    variation and the engine is over-smoothing it away. Measured Sep 2026: the
    per-peak baseline sat BELOW the input on envelope and jitter for most of the
    corpus (e.g. Female jitter input 2.14 vs engine 1.58), so a change that restored
    the real modulation was scored as 49 regressions while its reference-relative
    fidelity (residual_srr) improved by 4-18 dB on every class.

    This is the same failure mode the project hit with R3 and R7, where absolute and
    ratio metrics manufactured defects that direct measurement disproved. Use the
    deviation form for any class whose input is not intrinsically flat; keep the
    absolute form only for the steady-tone gate, where "no wobble" is the truth.

    Returns (dev_full, dev_band); 0 = the output modulates exactly like the input.
    """
    x, y = _align(_mono(reference), _mono(rendered))
    rf, rb = envelope_p2p(y, sr, band)
    xf, xb = envelope_p2p(x, sr, band)
    return abs(rf - xf), abs(rb - xb)


def trajectory_jitter_dev(rendered: np.ndarray, reference: np.ndarray, sr: int,
                          **kw) -> float:
    """Per-partial amplitude-jitter DEVIATION from the reference. See
    envelope_p2p_dev for why the absolute form is misleading off steady tones."""
    x, y = _align(_mono(reference), _mono(rendered))
    return abs(trajectory_jitter(y, sr, **kw) - trajectory_jitter(x, sr, **kw))


def _align(a: np.ndarray, b: np.ndarray, max_lag: int = 8192):
    """Shift b onto a by the integer lag that maximises correlation; equal-length out."""
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    c = np.fft.irfft(np.fft.rfft(a, nfft) * np.conj(np.fft.rfft(b, nfft)), nfft)
    ml = min(max_lag, n - 1)
    cc = np.concatenate([c[-ml:], c[:ml + 1]])
    lag = int(np.argmax(cc)) - ml
    if lag > 0:
        b = np.concatenate([np.zeros(lag), b])
    elif lag < 0:
        b = b[-lag:]
    n = min(len(a), len(b))
    return a[:n], b[:n]


_BANDS = [(0.0, 200.0), (200.0, 1000.0), (1000.0, 4000.0), (4000.0, 10000.0),
          (10000.0, 24000.0)]
BAND_NAMES = ["<200", "0.2-1k", "1-4k", "4-10k", ">10k"]


def residual_srr_local(rendered: np.ndarray, reference: np.ndarray, sr: int,
                       block_ms: float = 20.0):
    """Worst-case signal-to-residual ratio, in dB, over short blocks.

    THE PERCEPTUAL GATE. In the Sep 2026 blind test the four clips a listener could
    identify had a mean whole-file `residual_srr` of 21.3 dB and the eleven they
    could not had 21.5 -- identical. Fairlight C2 at 26.8 dB was heard; both cymbals
    at 8.8 and 10.2 were not. Whole-file SRR simply does not predict audibility,
    because an artifact confined to 6% of the file is averaged into invisibility.

    What did separate the two groups was the WORST moments: mean worst-block SRR of
    1.9 dB for the detected clips against 8.9 for the rest. So report the p10 and the
    minimum over 20 ms blocks, counting only blocks where the input is actually
    present. Higher is better, as with residual_srr.

    Returns (p10_db, worst_db).
    """
    x, y = _align(_mono(reference), _mono(rendered))
    if len(x) < 16 or y @ y <= 0:
        return 0.0, 0.0
    g = (x @ y) / (y @ y)
    e = x - g * y
    B = max(1, int(block_ms * 1e-3 * sr))
    n = min(len(x), len(e)) // B
    if n < 4:
        return 0.0, 0.0
    ex = (x[: n * B].reshape(n, B) ** 2).sum(axis=1)
    ee = (e[: n * B].reshape(n, B) ** 2).sum(axis=1)
    keep = ex > ex.max() * 1e-4          # ignore silence between notes
    if keep.sum() < 4:
        return 0.0, 0.0
    loc = 10.0 * np.log10(ex[keep] / np.maximum(ee[keep], 1e-30))
    return float(np.percentile(loc, 10)), float(loc.min())


def residual_srr(rendered: np.ndarray, reference: np.ndarray, sr: int,
                 per_band: bool = False):
    """Signal-to-residual ratio in dB: 10*log10(||x||^2 / ||x - g*y||^2).

    This is the residual meter. It measures e = x - x_hat on the REAL input, so it
    reports representability directly and cannot be gamed by a model that matches
    magnitudes while scrambling phase. Higher is better; +6 dB = half the residual
    energy. g is the least-squares scalar gain, so a pure level difference (the
    output limiter) does not count as a representation failure.

    Unity condition only -- it needs `reference` at the same pitch as `rendered`.

    Reference points measured on the r7-era engine at unity: steady tone ~24 dB,
    saw ~17-18, wavetable ~15, pluck ~13, chord ~10, voice ~10, then a cliff to
    cymbal ~3-4, choir ~3.5, dense mix ~2. The cliff is the thing to watch.

    Returns srr_db, or (srr_db, [per-band srr_db]) when per_band is set.
    """
    x, y = _align(_mono(reference), _mono(rendered))
    if len(x) < 16 or y @ y <= 0:
        return (0.0, [0.0] * len(_BANDS)) if per_band else 0.0
    g = (x @ y) / (y @ y)
    e = x - g * y
    total = float(10.0 * np.log10((x @ x) / (e @ e))) if e @ e > 0 else float("inf")
    if not per_band:
        return total
    X, Y = np.fft.rfft(x), np.fft.rfft(g * y)
    fr = np.fft.rfftfreq(len(x), 1.0 / sr)
    out = []
    for lo, hi in _BANDS:
        m = (fr >= lo) & (fr < hi)
        sx = float((np.abs(X[m]) ** 2).sum()) if m.any() else 0.0
        se = float((np.abs(X[m] - Y[m]) ** 2).sum()) if m.any() else 0.0
        out.append(float(10.0 * np.log10(sx / se)) if sx > 0 and se > 0 else float("nan"))
    return total, out


def residual_flatness(rendered: np.ndarray, reference: np.ndarray,
                      nfft: int = 2048) -> float:
    """Spectral flatness of the residual, divided by that of the input.

    Diagnoses WHAT the model is failing to capture, which determines the fix:
      ratio >> 1  : the residual is much flatter than the input -> what is left is
                    noise. Expected and benign; the tonal model did its job.
      ratio ~= 1  : the residual looks like the input -> the model is not capturing
                    the material's character at all. A TONAL failure, i.e. an
                    estimation defect, not a missing noise model.

    Measured on the r7-era engine: saw ~4.3, wavetable ~5-6, pluck ~3.4 (healthy);
    choir ~1.3, cymbal ~1.05, dense mix ~1.0 (the failing classes).
    """
    x, y = _align(_mono(reference), _mono(rendered))
    if len(x) < nfft * 2:
        return 0.0
    g = (x @ y) / (y @ y) if y @ y > 0 else 1.0
    e = x - g * y

    def flat(v):
        w = np.hanning(nfft)
        hop = nfft // 2
        fr = [np.abs(np.fft.rfft(v[i:i + nfft] * w)) + 1e-12
              for i in range(0, len(v) - nfft, hop)]
        if not fr:
            return 0.0
        m = np.array(fr)
        return float(np.median(np.exp(np.log(m).mean(1)) / m.mean(1)))

    fx = flat(x)
    return float(flat(e) / fx) if fx > 0 else 0.0


def envelope_p2p(x: np.ndarray, sr: int, band=(0.3, 1.0)):
    """Amplitude-envelope peak-to-peak, full-band and band-limited.

    Returns (p2p_full, p2p_band):
      p2p_full : (max-min)/mean of the analytic envelope over the steady body
                 of the signal. Steady-tone target <= 0.05 (5%).
      p2p_band : same but after band-passing the envelope to `band` Hz, which
                 isolates the slow 0.3-1 Hz "wobble" the plan calls out from
                 fast vibrato and from onset/offset transients.
    """
    y = _mono(x)
    if len(y) < sr // 4:
        return 0.0, 0.0
    env = np.abs(hilbert(y))
    # ignore the first/last 5% so onset/offset ramps don't dominate p2p
    k = max(1, len(env) // 20)
    body = env[k:-k]
    mean = body.mean() + 1e-12
    p2p_full = float((body.max() - body.min()) / mean)

    lo, hi = band
    nyq = sr / 2.0
    p2p_band = 0.0
    if hi < nyq and lo > 0:
        sos = butter(2, [lo / nyq, hi / nyq], btype="band", output="sos")
        eb = sosfiltfilt(sos, body - body.mean())
        p2p_band = float((eb.max() - eb.min()) / mean)
    return p2p_full, p2p_band
