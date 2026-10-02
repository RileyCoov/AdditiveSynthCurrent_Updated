#!/usr/bin/env python3
"""Two musical-noise metrics for comparing renders of the SAME source.

T -- temporal spikes ("blips"). A bin more than 10 dB above its OWN median over +/-5 frames
     (live bins within 70 dB of the max). T_blips uses STFT 1024 / hop 256 and is the metric that
     matched Riley's choice 6/6 on 2026-09-28 (toy renders) and 6/6 on 2026-09-29 (blind pairs).
     It is a CLICK detector: a line that persists for more than ~5 frames becomes its own
     median and never registers, so T is blind to steady whistling tones by construction.
     T_added uses STFT 2048 / hop 256 and additionally drops spikes that have a spike at the
     source frequency (output frequency mapped back through the shift ratio); it is the version
     that placed Riley's four marks at the 96-99.6th percentile. The two resolutions are kept
     exactly as they were when each was validated -- unifying them changed T_blips on 1985 down
     from 451.0 to 237.4, which would silently break comparison with every earlier result.

F -- spectral lines. A local frequency maximum more than 10 dB above the median of the 31
     surrounding bins in the SAME frame (STFT 1024 / hop 256). This sees a narrow line however
     long it lasts. "Invented" lines are those with no source line within +/-1 bin of the
     mapped frequency. F_tonal = invented lines that form a connected run (+/-1 bin drift) of
     at least KP frames (12 = 64 ms, more than two independent window lengths); F_fluct = the
     rest. Adapted from the 2026-10-01 external plan, section 4, with one deliberate change:
     that code maps EVERY source local maximum through the ratio. In a dense mix roughly every
     third bin is a local maximum, so after the +/-1-bin dilation nearly every output bin counts
     as "mapped" and F would read near zero on exactly the files that matter. Here the source
     side uses the same thresholded line definition as the output side.

CHANCE FLOORS matter more than anything else here. Pure white noise scored against an
independent white noise source gives T_added = 243/s -- higher than 1985's down-shift
baseline -- because a noisy output's random fluctuations never coincide with the source's.
So absolute T counts are NOT counts of invented components; only same-file A/B differences are
meaningful. `chance()` measures each floor; F_tonal's should be ~0 because chance rarely keeps a
>10 dB line in one bin for 12 consecutive frames.

Usage:
  python3 spikes.py <source.wav> <render.wav> <ratio>        # ratio 1.4983 up, 0.6674 down
  python3 spikes.py --chance
"""
import sys, struct
import numpy as np
from scipy.ndimage import median_filter, label, binary_dilation, find_objects

SR = 48000
KP = 12


def rdwav(p):
    """Robust WAV reader: 16/24/32-bit int and 32/64-bit float, any chunk order. The engine
    writes a header scipy.io.wavfile rejects, and lib.read_wav handles only 16-bit."""
    b = open(p, 'rb').read()
    assert b[:4] == b'RIFF' and b[8:12] == b'WAVE', p
    i, fmt, data = 12, None, None
    while i + 8 <= len(b):
        cid = b[i:i+4]; sz = struct.unpack('<I', b[i+4:i+8])[0]
        if cid == b'fmt ':   fmt = b[i+8:i+8+sz]
        elif cid == b'data': data = b[i+8:i+8+sz]
        i += 8 + sz + (sz & 1)
    af, ch, sr, _, _, bits = struct.unpack('<HHIIHH', fmt[:16])
    if af == 3:       x = np.frombuffer(data, '<f4' if bits == 32 else '<f8').astype(np.float64)
    elif bits == 16:  x = np.frombuffer(data, '<i2').astype(np.float64) / 32768.0
    elif bits == 32:  x = np.frombuffer(data, '<i4').astype(np.float64) / 2147483648.0
    elif bits == 24:
        a = np.frombuffer(data, np.uint8).reshape(-1, 3).astype(np.int32)
        v = a[:, 0] | (a[:, 1] << 8) | (a[:, 2] << 16)
        x = np.where(v & 0x800000, v - 0x1000000, v).astype(np.float64) / 8388608.0
    else:
        raise ValueError(f'{p}: format {af}, {bits} bits')
    if ch > 1: x = x.reshape(-1, ch).mean(1)
    return x, sr


def _stft_db(x, n, hop):
    w = np.hanning(n + 1)[:-1]
    nfr = 1 + (len(x) - n) // hop
    S = np.empty((n // 2 + 1, nfr))
    for i in range(nfr):
        S[:, i] = np.abs(np.fft.rfft(x[i*hop:i*hop+n] * w))
    return 20 * np.log10(S + 1e-12)


def _norm(x):
    return x / (np.abs(x).max() + 1e-12)


# ---------------------------------------------------------------- T: temporal spikes
T_N, T_HOP = 2048, 256        # T_added
TB_N = 1024                   # T_blips (the 6/6-validated resolution)

def _t_mask(L):
    n = L.shape[1]; med = np.empty_like(L)
    for t in range(n):
        a, b = max(0, t - 5), min(n, t + 6)
        med[:, t] = np.median(L[:, a:b], axis=1)
    live = L > (L.max() - 70.0)
    return ((L - med) > 10.0) & live


# ---------------------------------------------------------------- F: spectral lines
F_N, F_HOP, F_MED = 1024, 256, 31

def _f_mask(L):
    med = median_filter(L, size=(F_MED, 1), mode='nearest')
    S = np.zeros_like(L, bool)
    S[1:-1] = (L[1:-1] > L[:-2]) & (L[1:-1] > L[2:]) & ((L[1:-1] - med[1:-1]) > 10.0)
    return S


def _map_source(S_src, ratio, K):
    """Mark output bins that a source line maps to (source bin k -> output bin k*ratio)."""
    m = np.zeros((K, S_src.shape[1]), bool)
    ks, ts = np.nonzero(S_src)
    km = np.round(ks * ratio).astype(int); ok = km < K
    m[km[ok], ts[ok]] = True
    return binary_dilation(m, structure=np.ones((3, 1)))


def _tonal(inv, kp=KP):
    bridged = binary_dilation(inv, structure=np.ones((3, 1)))   # allow +/-1 bin drift
    lab, nl = label(bridged, structure=np.ones((3, 3)))
    keep = np.zeros(nl + 1, bool)
    for i, s in enumerate(find_objects(lab), 1):
        if s is not None and (s[1].stop - s[1].start) >= kp:
            keep[i] = True
    return keep[lab] & inv


class Source:
    """Cache a source file's masks once; score any number of renders against it."""
    def __init__(self, path):
        self.x, sr = rdwav(path); assert sr == SR, (path, sr)
        self.LT = _stft_db(_norm(self.x), T_N, T_HOP); self.BT = _t_mask(self.LT)
        self.LF = _stft_db(_norm(self.x), F_N, F_HOP); self.SF = _f_mask(self.LF)

    def score(self, path, ratio, return_masks=False):
        y, sr = rdwav(path); assert sr == SR, (path, sr)
        # T
        LT = _stft_db(_norm(y), T_N, T_HOP); BT = _t_mask(LT)
        nT = min(BT.shape[1], self.BT.shape[1]); durT = nT * T_HOP / SR
        fT = np.fft.rfftfreq(T_N, 1 / SR)
        sb = np.clip(np.round(fT / ratio / (SR / T_N)).astype(int), 0, len(fT) - 1)
        validT = ((fT / ratio) < SR / 2 * 0.95)[:, None]
        t_added = BT[:, :nT] & ~self.BT[sb, :][:, :nT] & validT
        # F
        LF = _stft_db(_norm(y), F_N, F_HOP); SF = _f_mask(LF)
        nF = min(SF.shape[1], self.SF.shape[1]); durF = nF * F_HOP / SR
        mapped = _map_source(self.SF[:, :nF], ratio, SF.shape[0])
        inv = SF[:, :nF] & ~mapped
        ton = _tonal(inv)
        LB = _stft_db(_norm(y), TB_N, T_HOP); BB = _t_mask(LB)
        out = dict(T_blips=BB.sum() / (BB.shape[1] * T_HOP / SR),
                   T_added=t_added.sum() / durT,
                   F_lines=SF[:, :nF].sum() / durF,
                   F_tonal=ton.sum() / durF,
                   F_fluct=(inv & ~ton).sum() / durF)
        if return_masks:
            return out, dict(t_added=t_added, f_tonal=ton, f_inv=inv,
                             fT=fT, fF=np.fft.rfftfreq(F_N, 1 / SR))
        return out


def chance(seconds=5, seed=1):
    """Floors: white noise scored against an INDEPENDENT white noise source at ratio 1."""
    import tempfile, os
    rng = np.random.default_rng(seed)
    d = tempfile.mkdtemp()
    def w(p, x):
        x = np.clip(x, -1, 1); b = (x * 32767).astype('<i2').tobytes()
        open(p, 'wb').write(b'RIFF' + struct.pack('<I', 36 + len(b)) + b'WAVEfmt ' +
            struct.pack('<IHHIIHH', 16, 1, 1, SR, SR * 2, 2, 16) + b'data' +
            struct.pack('<I', len(b)) + b)
    a, b = os.path.join(d, 'a.wav'), os.path.join(d, 'b.wav')
    w(a, 0.1 * rng.standard_normal(SR * seconds)); w(b, 0.1 * rng.standard_normal(SR * seconds))
    return Source(a).score(b, 1.0)


def fmt(d):
    return ('T_blips {T_blips:7.1f}  T_added {T_added:7.1f}  |  F_lines {F_lines:7.1f}  '
            'F_tonal {F_tonal:6.2f}  F_fluct {F_fluct:7.1f}').format(**d)


if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] == '--chance':
        print('chance floor (white noise vs independent white noise):'); print(' ', fmt(chance()))
    elif len(sys.argv) == 4:
        print(fmt(Source(sys.argv[1]).score(sys.argv[2], float(sys.argv[3]))))
    else:
        print(__doc__)
