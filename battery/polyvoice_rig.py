#!/usr/bin/env python3
"""Polyphonic-vibrato rig — the gate for the choir gap (docs 4d.4 item 1).

`vibrato_rig.py` covers ONE voice with vibrato, which pitch-synchronous analysis
already fixes (bacfed2) by warping time so the single f0 is constant. A choir has no
single f0: several voices, each with its own vibrato rate and phase, so the warp
cannot engage and correctly refuses to. The measured consequence is choir sitting
~13 dB below its oracle ceiling with a TONAL residual (flatness 0.076) -- missing or
mis-estimated partials, not missing noise.

This rig separates the two possible causes:

  POLY_STEADY  three steady voices        -> is plain polyphony the problem?
  ONE_VIB      one voice with vibrato     -> the case pitch-sync already handles
  POLY_VIB     three voices, INDEPENDENT
               vibrato rates and phases   -> the choir case

If POLY_STEADY scores well and POLY_VIB badly, the defect is per-partial motion under
polyphony and the fix is per-track demodulation (estimate each track on its own
smoothed frequency trajectory instead of warping the whole signal).

Metrics: amp_energy (reconstructed energy per known partial / input), shape_corr
(windowed best-lag waveform correlation) and srr_db (signal-to-residual on the
delay-aligned output, the same meter the battery uses).

Usage:
    python3 battery/polyvoice_rig.py
    python3 battery/polyvoice_rig.py --binary build/additive --args "0 0 0.0 0 0 -1 0"
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile, wave
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from metrics import waveform_correlation

SR = 48000
DUR = 3.0
NHARM = 12
# A major triad, the shape of choir-burst_C_major: (f0, vib rate Hz, vib depth, vib phase)
VOICES = [(220.00, 5.5, 0.030, 0.0),
          (277.18, 6.3, 0.025, 1.7),
          (329.63, 4.8, 0.035, 3.0)]


def _write(path, x):
    x = x / max(1e-9, np.max(np.abs(x))) * 0.7
    with wave.open(path, "wb") as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(SR)
        w.writeframes((np.clip(x, -1, 1) * 32767).astype("<i2").tobytes())


def _read(path):
    with wave.open(path, "rb") as w:
        return np.frombuffer(w.readframes(w.getnframes()), "<i2").astype(float) / 32768.0


def _voice(t, f0, rate, depth, phase, vibrato):
    """One harmonic stack; every partial follows the SAME f0 contour, as a real voice does."""
    contour = 1 + depth * np.sin(2 * np.pi * rate * t + phase) if vibrato else np.ones_like(t)
    y = np.zeros_like(t)
    for k in range(1, NHARM + 1):
        y += (1.0 / k) * np.sin(2 * np.pi * np.cumsum(f0 * k * contour) / SR)
    return y


def make_signals(tmp):
    t = np.arange(int(DUR * SR)) / SR
    sigs = {
        "POLY_STEADY": sum(_voice(t, *v, vibrato=False) for v in VOICES),
        "ONE_VIB":     _voice(t, *VOICES[0], vibrato=True),
        "POLY_VIB":    sum(_voice(t, *v, vibrato=True) for v in VOICES),
    }
    out = {}
    for name, y in sigs.items():
        p = os.path.join(tmp, f"{name}.wav"); _write(p, y); out[name] = p
    return out


def partials(name):
    vs = [VOICES[0]] if name == "ONE_VIB" else VOICES
    return [f0 * k for f0, _, _, _ in vs for k in range(1, NHARM + 1)]


def amp_energy(inp, out, freqs):
    xi, xo = _read(inp), _read(out)
    def e(y, f, wid=30):
        Y = np.abs(np.fft.rfft(y * np.hanning(len(y)))); b = int(round(f * len(y) / SR))
        return np.sqrt((Y[max(0, b - wid):b + wid] ** 2).sum())
    return float(np.mean([e(xo, f) / e(xi, f) for f in freqs if e(xi, f) > 0]))


def shape_corr(inp, out, win=2048, hop=1024):
    xi, xo = _read(inp), _read(out); n = min(len(xi), len(xo)); xi, xo = xi[:n], xo[:n]
    cs = [waveform_correlation(xo[s:s + win], xi[s:s + win], max_lag=64)
          for s in range(win, n - 2 * win, hop) if np.sqrt(np.mean(xi[s:s + win] ** 2)) > 1e-3]
    return float(np.mean(cs)) if cs else 0.0


def srr_db(inp, out, max_lag=4096):
    """Signal-to-residual after removing the engine's bulk delay and fitting one gain."""
    xi, xo = _read(inp), _read(out); n = min(len(xi), len(xo)); xi, xo = xi[:n], xo[:n]
    c = np.correlate(xo[: 4 * SR // 2], xi[: 4 * SR // 2], "full")
    lag = int(np.argmax(np.abs(c))) - (len(xi[: 4 * SR // 2]) - 1)
    lag = int(np.clip(lag, -max_lag, max_lag))
    if lag > 0:   xo2, xi2 = xo[lag:], xi[: n - lag]
    elif lag < 0: xo2, xi2 = xo[: n + lag], xi[-lag:]
    else:         xo2, xi2 = xo, xi
    m = min(len(xi2), len(xo2)); xi2, xo2 = xi2[:m], xo2[:m]
    if xo2 @ xo2 <= 0: return float("nan")
    e = xi2 - ((xi2 @ xo2) / (xo2 @ xo2)) * xo2
    return float(10 * np.log10((xi2 @ xi2) / max(e @ e, 1e-30)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default="build/additive")
    ap.add_argument("--block", default="2048")
    ap.add_argument("--args", default="0 0 0.0 0 0 -1 0")
    ap.add_argument("--keep", metavar="DIR", help="write the rig wavs and renders here")
    a = ap.parse_args()
    extra = a.args.split()
    tmp = a.keep or tempfile.mkdtemp()
    os.makedirs(tmp, exist_ok=True)
    sigs = make_signals(tmp)
    print(f"{'signal':<14}{'amp_energy':>12}{'shape_corr':>12}{'srr_db':>10}")
    for name, inp in sigs.items():
        out = os.path.join(tmp, name + "_out.wav")
        subprocess.run([a.binary, inp, out, a.block, *extra], check=True, stdout=subprocess.DEVNULL)
        print(f"{name:<14}{amp_energy(inp, out, partials(name)):>12.3f}"
              f"{shape_corr(inp, out):>12.3f}{srr_db(inp, out):>10.2f}")
    if a.keep:
        print(f"\n  wavs kept in {tmp}")
    print("\n  POLY_STEADY good + POLY_VIB bad  -> per-partial motion under polyphony,")
    print("  i.e. the per-track demodulation case. ONE_VIB is the pitch-sync control.")


if __name__ == "__main__":
    main()
