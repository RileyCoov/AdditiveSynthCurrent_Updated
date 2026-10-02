#!/usr/bin/env python3
"""Phase-randomised surrogate: same short-time magnitude, no phase coherence.

The null hypothesis for "is this peak a sinusoid?". A surrogate has, by construction, no
sinusoids -- only narrowband noise with the original's spectral envelope and time envelope.
Any statistic that cannot separate the original from the surrogate is not measuring
sinusoidality.

Two deliberate differences from the ad-hoc version this project used earlier
(4 inline lines in ~/Downloads/tracksim.py):

 1. Randomise at the ANALYSIS resolution (4096 / hop 1024), not 1024 / hop 256. The engine
    analyses at 4096; randomising at 1024 destroys coherence only at a 46.9 Hz scale and
    leaves ~47 Hz-wide narrowband structures that a 4096 analysis still resolves as peaks.
    That makes the null too easy to distinguish for the wrong reason. `--nfft 1024` keeps
    the old behaviour for comparability.
 2. Report per-band magnitude preservation. An STFT with i.i.d. random phase is INCONSISTENT:
    no signal has that STFT, so ISTFT returns a least-squares compromise whose short-time
    magnitude is a smoothed version of the target. Nobody had checked how much. If the
    surrogate is not level-matched per band, every false-alarm rate calibrated on it is wrong.

Usage:
  python3 surrogate.py <in.wav> [--seeds 20] [--nfft 4096] [--hop 1024] [--out DIR]
"""
import sys, os, argparse, struct
import numpy as np

SR_EXPECT = 48000

def rdwav(p):
    b = open(p, 'rb').read()
    assert b[:4] == b'RIFF' and b[8:12] == b'WAVE', p
    i, fmt, data = 12, None, None
    while i + 8 <= len(b):
        cid = b[i:i+4]; sz = struct.unpack('<I', b[i+4:i+8])[0]
        if cid == b'fmt ':  fmt = b[i+8:i+8+sz]
        elif cid == b'data': data = b[i+8:i+8+sz]
        i += 8 + sz + (sz & 1)
    af, ch, sr, _, ba, bits = struct.unpack('<HHIIHH', fmt[:16])
    if af == 3:     x = np.frombuffer(data, '<f4' if bits == 32 else '<f8')
    elif bits == 16: x = np.frombuffer(data, '<i2').astype(np.float64) / 32768.0
    elif bits == 32: x = np.frombuffer(data, '<i4').astype(np.float64) / 2147483648.0
    elif bits == 24:
        a = np.frombuffer(data, np.uint8).reshape(-1, 3).astype(np.int32)
        v = a[:, 0] | (a[:, 1] << 8) | (a[:, 2] << 16)
        x = np.where(v & 0x800000, v - 0x1000000, v).astype(np.float64) / 8388608.0
    else: raise ValueError(f'{p}: af={af} bits={bits}')
    x = x.astype(np.float64)
    if ch > 1: x = x.reshape(-1, ch).mean(1)
    return sr, x

def wrwav(p, x, sr):
    x = np.clip(x, -1.0, 1.0)
    d = (x * 32767).astype('<i2').tobytes()
    hdr = (b'RIFF' + struct.pack('<I', 36 + len(d)) + b'WAVEfmt ' +
           struct.pack('<IHHIIHH', 16, 1, 1, sr, sr * 2, 2, 16) + b'data' +
           struct.pack('<I', len(d)))
    open(p, 'wb').write(hdr + d)

def stft(x, n, hop):
    w = np.hanning(n + 1)[:-1]
    nfr = 1 + max(0, (len(x) - n)) // hop
    Z = np.empty((n // 2 + 1, nfr), complex)
    for i in range(nfr):
        Z[:, i] = np.fft.rfft(x[i*hop:i*hop+n] * w)
    return Z, w

def istft(Z, n, hop, w, L):
    """Weighted overlap-add. Divides by the summed squared window, which is the
    least-squares inverse for an inconsistent STFT.

    The window sum decays to ~0 in the first and last window, so dividing by it raw
    produces enormous edge samples -- which then crush the entire file through any
    subsequent peak normalisation (this cost one debugging round: every band read exactly
    -46 dB, a uniform offset that looked like a banding bug and was an edge-division bug).
    Floored at a fraction of its own maximum, as the engine's residual OLA does."""
    y = np.zeros(L + n); ws = np.zeros(L + n)
    for i in range(Z.shape[1]):
        seg = np.fft.irfft(Z[:, i], n)
        y[i*hop:i*hop+n]  += seg * w
        ws[i*hop:i*hop+n] += w * w
    return (y / np.maximum(ws, 0.1 * ws.max()))[:L]

def make(x, n, hop, seed):
    Z, w = stft(x, n, hop)
    rng = np.random.default_rng(seed)
    Zs = np.abs(Z) * np.exp(2j * np.pi * rng.random(Z.shape))
    return istft(Zs, n, hop, w, len(x))

BANDS = [(0,500),(500,2000),(2000,5000),(5000,8000),(8000,12000),(12000,16000),(16000,24000)]

def band_mag_dev(x, y, n=4096, hop=1024):
    """Per-band RMS deviation of short-time MAGNITUDE, measured at the ANALYSIS resolution --
    the quantity that has to survive for the surrogate to be a fair null."""
    Zx, _ = stft(x, n, hop); Zy, _ = stft(y, n, hop)
    m = min(Zx.shape[1], Zy.shape[1])
    f = np.fft.rfftfreq(n, 1.0 / SR_EXPECT)
    out = []
    for lo, hi in BANDS:
        sel = (f >= lo) & (f < hi)
        if not sel.any(): out.append((lo, hi, float('nan'), float('nan'))); continue
        ax = np.abs(Zx[sel, :m]); ay = np.abs(Zy[sel, :m])
        # overall band level difference, and per-frame scatter around it
        ex = (ax ** 2).sum(); ey = (ay ** 2).sum()
        lvl = 10 * np.log10(max(ey, 1e-30) / max(ex, 1e-30))
        px = (ax ** 2).sum(0); py = (ay ** 2).sum(0)
        ok = px > px.max() * 1e-6
        scat = np.std(10 * np.log10(np.maximum(py[ok], 1e-30) / np.maximum(px[ok], 1e-30))) if ok.any() else float('nan')
        out.append((lo, hi, lvl, scat))
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('wav'); ap.add_argument('--seeds', type=int, default=20)
    ap.add_argument('--nfft', type=int, default=4096); ap.add_argument('--hop', type=int, default=1024)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    sr, x = rdwav(a.wav)
    assert sr == SR_EXPECT, f'expected {SR_EXPECT}, got {sr}'
    base = os.path.basename(a.wav)[:-4]
    outdir = a.out or os.path.join(os.path.dirname(os.path.abspath(a.wav)), 'surrogate')
    os.makedirs(outdir, exist_ok=True)
    print(f'{base}: {len(x)} samples, randomising at nfft={a.nfft} hop={a.hop}, {a.seeds} seeds')
    devs = []
    for s in range(a.seeds):
        y = make(x, a.nfft, a.hop, s)
        # match overall level so the engine's max-based detection threshold sees the same scale
        y *= np.sqrt(np.mean(x ** 2) / max(np.mean(y ** 2), 1e-30))
        pk = np.abs(y).max()
        if pk > 0.98: y *= 0.98 / pk
        wrwav(os.path.join(outdir, f'{base}__surr_n{a.nfft}_s{s:02d}.wav'), y, sr)
        devs.append(band_mag_dev(x, y, a.nfft, a.hop))
    print(f'  wrote {a.seeds} files to {outdir}/')
    print(f'\n  PER-BAND MAGNITUDE PRESERVATION (mean over {a.seeds} seeds)')
    print(f'  {"band Hz":>14s} {"level diff dB":>14s} {"per-frame scatter dB":>21s}   verdict')
    for bi, (lo, hi, _, _) in enumerate(devs[0]):
        lv = np.nanmean([d[bi][2] for d in devs]); sc = np.nanmean([d[bi][3] for d in devs])
        v = 'ok' if abs(lv) < 1.0 and sc < 3.0 else ('LEVEL OFF' if abs(lv) >= 1.0 else 'SCATTERED')
        print(f'  {lo:6d}-{hi:6d} {lv:+14.2f} {sc:21.2f}   {v}')
