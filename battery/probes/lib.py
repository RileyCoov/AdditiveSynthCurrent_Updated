"""Shared analysis helpers for the Stage-0 probes."""
import numpy as np, wave


def read_wav(path):
    w = wave.open(path, 'rb')
    n, sr, ch, sw = w.getnframes(), w.getframerate(), w.getnchannels(), w.getsampwidth()
    raw = w.readframes(n); w.close()
    assert sw == 2, f"{path}: expected 16-bit, got {sw*8}"
    d = np.frombuffer(raw, dtype='<i2').astype(np.float64) / 32768.0
    if ch > 1:
        d = d.reshape(-1, ch).mean(axis=1)
    return d, sr


def best_lag(x, y, max_lag=8192):
    """Integer lag that maximises correlation of y against x."""
    m = min(len(x), len(y))
    a, b = x[:m], y[:m]
    n = 1 << int(np.ceil(np.log2(2 * m)))
    c = np.fft.irfft(np.fft.rfft(a, n) * np.conj(np.fft.rfft(b, n)), n)
    c = np.concatenate([c[-max_lag:], c[:max_lag + 1]])
    return int(np.argmax(c)) - max_lag


def align(x, y, max_lag=8192):
    """Shift y to best match x; return equal-length (x, y)."""
    L = best_lag(x, y, max_lag)
    if L > 0:
        y = np.concatenate([np.zeros(L), y])
    elif L < 0:
        y = y[-L:]
    m = min(len(x), len(y))
    return x[:m], y[:m]


def srr(x, y, fit_gain=True):
    """Signal-to-residual ratio in dB. y is the reconstruction of x."""
    g = (x @ y) / (y @ y) if (fit_gain and y @ y > 0) else 1.0
    e = x - g * y
    if e @ e <= 0:
        return np.inf, g
    return 10 * np.log10((x @ x) / (e @ e)), g


BANDS = [(0, 200), (200, 1000), (1000, 4000), (4000, 10000), (10000, 24000)]
BAND_NAMES = ['<200', '.2-1k', '1-4k', '4-10k', '>10k']


def band_srr(x, y, sr, g=None):
    """Per-band SRR via rFFT masking of the whole signal."""
    if g is None:
        g = (x @ y) / (y @ y) if y @ y > 0 else 1.0
    X, Y = np.fft.rfft(x), np.fft.rfft(y * g)
    fr = np.fft.rfftfreq(len(x), 1 / sr)
    out = []
    for lo, hi in BANDS:
        m = (fr >= lo) & (fr < hi)
        if not m.any():
            out.append(np.nan); continue
        sx = (np.abs(X[m]) ** 2).sum()
        se = (np.abs(X[m] - Y[m]) ** 2).sum()
        out.append(10 * np.log10(sx / se) if se > 0 and sx > 0 else np.nan)
    return out, g


def stft_mag(x, N=2048, hop=None):
    hop = hop or N // 2
    w = np.hanning(N)
    idx = range(0, max(1, len(x) - N), hop)
    return np.array([np.abs(np.fft.rfft(x[i:i + N] * w)) for i in idx]) + 1e-12


def spectral_flatness(mag):
    """Per-frame flatness (geometric/arithmetic mean). 1 = white, 0 = tonal."""
    return np.exp(np.log(mag).mean(axis=1)) / mag.mean(axis=1)


def waveform_corr(x, y, max_lag=8192):
    x, y = align(x, y, max_lag)
    xa, ya = x - x.mean(), y - y.mean()
    d = np.linalg.norm(xa) * np.linalg.norm(ya)
    return float(xa @ ya / d) if d > 0 else 0.0
