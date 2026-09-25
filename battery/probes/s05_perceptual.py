"""Stage B1 - a perceptually weighted residual measure.

Every metric in this project is UNWEIGHTED, and that is why four consecutive rounds of
measurable improvement were inaudible (docs 5b.1) and why up-shift now measures BETTER
than down-shift while listeners report the reverse (docs 5c).

Two measures here, because shifted audio has no waveform reference:

  unity_weighted_srr   : residual e = x - g*y, but scored per Bark band against the
                         input's own masked threshold and weighted by an equal-loudness
                         curve, so an error at 3 kHz counts for more than the same error
                         at 200 Hz.
  shifted_weighted_dev : the output's short-time Bark spectrum against the input's
                         Bark spectrum TRANSPOSED by the shift ratio (docs 5c.1: band
                         comparisons must transpose). Magnitude-domain, so weaker than a
                         residual, but it is the only reference shifted audio has.

GATE (this is a falsification test, not a feature): the measure must reproduce judgements
we already have. It must rank the four clips a listener identified in round 1 worse than
the eleven they could not, and it must rank up-shift worse than down-shift on the drum
loop and piano. If it cannot, it is not the instrument and should not be trusted.
"""
import numpy as np

# Bark band edges (Zwicker), Hz
BARK = [20,100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,2320,2700,
        3150,3700,4400,5300,6400,7700,9500,12000,15500,20000]


def _elc_weight_db(f):
    """Equal-loudness-ish weighting: how audible is a given frequency, dB, peak near 3-4 kHz.
    A smooth approximation of the 60-phon contour inverted; exact values do not matter
    much, the shape does."""
    f = np.maximum(np.asarray(f, dtype=float), 20.0)
    # A-weighting magnitude (dB), a standard and defensible stand-in
    f2 = f * f
    ra = (12194.0**2 * f2**2) / ((f2 + 20.6**2) * np.sqrt((f2 + 107.7**2) * (f2 + 737.9**2)) * (f2 + 12194.0**2))
    return 20.0 * np.log10(ra) + 2.0


def _bark_frames(x, sr, n=2048, hop=512):
    w = np.hanning(n)
    nf = max(1, 1 + (len(x) - n) // hop)
    f = np.fft.rfftfreq(n, 1.0 / sr)
    idx = [np.where((f >= lo) & (f < hi))[0] for lo, hi in zip(BARK[:-1], BARK[1:])]
    out = np.zeros((nf, len(idx)))
    for i in range(nf):
        S = np.abs(np.fft.rfft(x[i * hop:i * hop + n] * w)) ** 2
        for b, ix in enumerate(idx):
            out[i, b] = S[ix].sum() if len(ix) else 0.0
    centres = np.array([0.5 * (lo + hi) for lo, hi in zip(BARK[:-1], BARK[1:])])
    return out, centres


def _spread(B):
    """Simple Bark spreading function -> masking threshold from the signal itself.
    Slopes: +25 dB/Bark on the low side of a masker, -10 dB/Bark on the high side
    (maskers hide content ABOVE them far more than below). This is what a weighting
    curve alone cannot capture, and it is why broadband material hides its own error:
    a cymbal's residual sits under the cymbal's own mask."""
    nb = B.shape[1]
    k = np.arange(nb)
    # matrix of spreading gains, rows = masker band, cols = masked band
    d = k[None, :] - k[:, None]
    slope = np.where(d >= 0, -10.0 * d, 25.0 * d)      # dB, d<0 means below the masker
    G = 10.0 ** (slope / 10.0)
    return B @ G


def _absolute_threshold_power(c, full_scale_spl=90.0):
    """Absolute threshold of hearing per Bark band, as a digital POWER level.

    Without this a quiet error in an EMPTY band counts as audible when it is in fact
    below audibility -- which is exactly what happened on the 300 Hz sine: masking
    alone ranked it the worst file in the corpus while no listener could hear it at
    42.7 dB SRR. Terhardt's approximation of the threshold in dB SPL, referred to
    digital full scale at `full_scale_spl`.
    """
    f = np.asarray(c, dtype=float) / 1000.0
    spl = 3.64 * f ** -0.8 - 6.5 * np.exp(-0.6 * (f - 3.3) ** 2) + 1e-3 * f ** 4
    return 10.0 ** ((spl - full_scale_spl) / 10.0)


def noise_to_mask(x, y, sr, offset_db=4.0):
    """Noise-to-mask ratio, dB. LOWER is better; <0 means the error is below the
    masked threshold and should be inaudible. The standard perceptual-codec measure."""
    m = min(len(x), len(y)); x, y = x[:m], y[:m]
    g = (x @ y) / (y @ y) if y @ y > 0 else 1.0
    e = x - g * y
    Bx, c = _bark_frames(x, sr)
    Be, _ = _bark_frames(e, sr)
    thr = _spread(Bx) * 10.0 ** (-offset_db / 10.0)
    # the threshold is the LOUDER of the masked threshold and absolute audibility
    nrm = Bx.sum(axis=1).max()
    atp = _absolute_threshold_power(c) * max(nrm, 1e-30) / max(Bx.shape[1], 1)
    thr = np.maximum(thr, atp[None, :])
    active = Bx.sum(axis=1) > Bx.sum(axis=1).max() * 1e-4
    if active.sum() < 4:
        return 0.0, 0.0
    nmr = Be[active] / np.maximum(thr[active], 1e-30)
    per_frame = 10.0 * np.log10(np.maximum(nmr.mean(axis=1), 1e-30))
    return float(10.0 * np.log10(max(nmr.mean(), 1e-30))), float(np.percentile(per_frame, 95))


def unity_weighted_srr(x, y, sr):
    """Bark-band, equal-loudness-weighted signal-to-residual ratio, dB. Higher is better."""
    m = min(len(x), len(y)); x, y = x[:m], y[:m]
    g = (x @ y) / (y @ y) if y @ y > 0 else 1.0
    e = x - g * y
    Bx, c = _bark_frames(x, sr)
    Be, _ = _bark_frames(e, sr)
    w = 10.0 ** (_elc_weight_db(c) / 10.0)          # power weighting
    active = Bx.sum(axis=1) > Bx.sum(axis=1).max() * 1e-4
    if active.sum() < 4:
        return 0.0
    num = (Bx[active] * w).sum()
    den = (Be[active] * w).sum()
    return float(10.0 * np.log10(num / max(den, 1e-30)))


def shifted_weighted_dev(x, y, sr, semitones):
    """Weighted Bark deviation of shifted output from the input's TRANSPOSED spectrum, dB.
    Lower is better. Reported as rms over active frames and the p95 worst frame."""
    ratio = 2.0 ** (semitones / 12.0)
    Bx, c = _bark_frames(x, sr)
    By, _ = _bark_frames(y, sr)
    nf = min(len(Bx), len(By)); Bx, By = Bx[:nf], By[:nf]
    # transpose the input's Bark profile: band at centre c should move to c*ratio
    tgt = np.zeros_like(Bx)
    for b, cc in enumerate(c):
        dest = cc * ratio
        j = np.searchsorted(c, dest)
        if 0 < j < len(c):
            f0, f1 = c[j - 1], c[j]
            a = (dest - f0) / max(f1 - f0, 1e-9)
            tgt[:, j - 1] += Bx[:, b] * (1 - a)
            tgt[:, j] += Bx[:, b] * a
        elif j == 0:
            tgt[:, 0] += Bx[:, b]
        # content transposed past the top band has nowhere to go: dropped, correctly
    w = 10.0 ** (_elc_weight_db(c) / 10.0)
    active = (tgt.sum(axis=1) > tgt.sum(axis=1).max() * 1e-4)
    if active.sum() < 4:
        return 0.0, 0.0
    lt = 10 * np.log10(tgt[active] + 1e-20)
    ly = 10 * np.log10(By[active] + 1e-20)
    # remove a single global gain so level is not counted as an error
    off = np.average(ly - lt, weights=np.broadcast_to(w, ly.shape))
    d = np.abs(ly - lt - off)
    per_frame = np.sqrt(np.average(d ** 2, axis=1, weights=w))
    return float(np.sqrt(np.average(d ** 2, weights=np.broadcast_to(w, d.shape)))), float(np.percentile(per_frame, 95))


def shifted_nmr(x, y, sr, semitones, offset_db=4.0):
    """Masked deviation of shifted output from the input's TRANSPOSED Bark spectrum, dB.

    Shifted audio has no waveform reference, so this is magnitude-domain: build the
    target by transposing the input's Bark profile (docs 5c.1 - band comparisons must
    transpose), then score the per-band shortfall/excess against the target's own
    masked threshold plus absolute audibility. Lower is better.
    """
    ratio = 2.0 ** (semitones / 12.0)
    Bx, c = _bark_frames(x, sr)
    By, _ = _bark_frames(y, sr)
    nf = min(len(Bx), len(By)); Bx, By = Bx[:nf], By[:nf]
    tgt = np.zeros_like(Bx)
    for b, cc in enumerate(c):
        dest = cc * ratio
        j = np.searchsorted(c, dest)
        if 0 < j < len(c):
            f0, f1 = c[j - 1], c[j]
            a = (dest - f0) / max(f1 - f0, 1e-9)
            tgt[:, j - 1] += Bx[:, b] * (1 - a)
            tgt[:, j] += Bx[:, b] * a
        elif j == 0:
            tgt[:, 0] += Bx[:, b]
    active = tgt.sum(axis=1) > tgt.sum(axis=1).max() * 1e-4
    if active.sum() < 4:
        return 0.0, 0.0
    # single global gain so level is not scored as an error
    gt = tgt[active].sum() / max(By[active].sum(), 1e-30)
    Bys = By[active] * gt
    T = tgt[active]
    thr = _spread(T) * 10.0 ** (-offset_db / 10.0)
    nrm = T.sum(axis=1).max()
    atp = _absolute_threshold_power(c) * max(nrm, 1e-30) / max(len(c), 1)
    thr = np.maximum(thr, atp[None, :])
    dev = (np.sqrt(np.maximum(Bys, 0)) - np.sqrt(np.maximum(T, 0))) ** 2
    nmr = dev / np.maximum(thr, 1e-30)
    per_frame = 10.0 * np.log10(np.maximum(nmr.mean(axis=1), 1e-30))
    return float(10.0 * np.log10(max(nmr.mean(), 1e-30))), float(np.percentile(per_frame, 95))
