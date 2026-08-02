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
