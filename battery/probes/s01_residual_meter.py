"""Stage 0.1 - residual meter.

Measures e = x - x_hat on the REAL input for every corpus file:
  - overall signal-to-residual ratio (SRR, dB)  = how much of the signal the model captures
  - per-band SRR
  - residual spectral flatness                  = is what's left noise-like or tonal
  - passthrough fraction                        = how much of the "reconstruction" is
    literally the original waveform pasted back in by the unity transient blend
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from lib import (read_wav, align, srr, band_srr, stft_mag, spectral_flatness,
                 waveform_corr, BAND_NAMES)

IN = os.path.expanduser('~/Desktop/comparingOut/inputData')
FULL = '/tmp/stage0/render'      # engine as shipped (blend + stochastic residual)
NORES = '/tmp/stage0/render_nr'  # residual fill off -> blend regions are exact copies

CLASS = {
    '300hzSine_5sec': 'steady', '300hzSaw_5sec': 'saw', '440sawtooth': 'saw',
    'FairlightAdditiveWaveTA1_C2': 'wavetable', 'FairlightAdditiveWaveTA1_C3': 'wavetable',
    'PianoSampleMono': 'pluck', 'SaintSaensMonoChord': 'chord',
    'Female_Sung_Line_3_48': 'voice', 'choir-burst_C_major_Mono_48k': 'choir',
    'DrumLoopShort': 'perc', '48kCymbal': 'cymbal', 'Out48kCymbal': 'cymbal',
    '1985': 'mix', 'HappyMono': 'mix', 'take-me-out': 'mix',
}

rows = []
for name in sorted(CLASS, key=lambda n: (CLASS[n], n)):
    x, sr = read_wav(f'{IN}/{name}.wav')
    y, _ = read_wav(f'{FULL}/{name}_orig.wav')
    yn, _ = read_wav(f'{NORES}/{name}_orig.wav')
    x_a, y_a = align(x, y)
    xn_a, yn_a = align(x, yn)

    # passthrough mask: with the stochastic fill off, the unity transient blend
    # writes the ORIGINAL samples verbatim -> |x - y| is ~0 there.
    m = min(len(xn_a), len(yn_a))
    d = np.abs(xn_a[:m] - yn_a[:m])
    scale = np.abs(xn_a[:m]).max()
    pass_mask = d < (1e-4 * scale)
    # only count it where the signal is actually present
    pass_mask &= np.abs(xn_a[:m]) > (1e-3 * scale)
    pass_frac = float(pass_mask.mean())

    s_all, g = srr(x_a, y_a)
    bands, _ = band_srr(x_a, y_a, sr, g)

    # SRR restricted to NON-passthrough samples = the honest additive number
    keep = ~pass_mask
    if keep.sum() > 1000:
        s_syn, _ = srr(xn_a[:m][keep], yn_a[:m][keep])
    else:
        s_syn = np.nan

    e = x_a - g * y_a
    fl_e = float(np.median(spectral_flatness(stft_mag(e))))
    fl_x = float(np.median(spectral_flatness(stft_mag(x_a))))

    rows.append(dict(name=name, cls=CLASS[name], sr=sr, srr=s_all, srr_synth_only=s_syn,
                     gain=g, pass_frac=pass_frac, bands=bands,
                     flat_resid=fl_e, flat_input=fl_x,
                     wcorr=waveform_corr(x, y)))

hdr = (f"{'file':30s} {'class':9s} {'SRR':>6s} {'SRRsyn':>7s} {'pass%':>6s} {'wcorr':>6s} "
       + ' '.join(f'{b:>6s}' for b in BAND_NAMES) + f" {'flatE':>6s} {'flatX':>6s}")
print(hdr); print('-' * len(hdr))
for r in rows:
    b = ' '.join(f'{v:6.1f}' if np.isfinite(v) else '    --' for v in r['bands'])
    ss = f"{r['srr_synth_only']:7.1f}" if np.isfinite(r['srr_synth_only']) else '     --'
    print(f"{r['name']:30s} {r['cls']:9s} {r['srr']:6.1f} {ss} "
          f"{100*r['pass_frac']:6.1f} {r['wcorr']:6.3f} {b} {r['flat_resid']:6.3f} {r['flat_input']:6.3f}")

json.dump(rows, open('/tmp/stage0/s01_residual_meter.json', 'w'), indent=1)
print('\nSRR = signal-to-residual ratio, dB (higher = more of the signal is captured).')
print('SRRsyn = SRR excluding samples the unity transient blend pastes from the original.')
print('pass%  = share of active samples that ARE the pasted original.')
