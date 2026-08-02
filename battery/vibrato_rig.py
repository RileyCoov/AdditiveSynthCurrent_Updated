#!/usr/bin/env python3
"""Controlled vibrato rig — the cloud-testable gate for the voice artifact.

The Female "phasey/double voice" was localized (see docs/voice-artifact-and-plan.md)
to the analysis mis-measuring partials that MOVE under vibrato. This rig makes that
measurable without ears or real recordings: it synthesizes known equal-amplitude
harmonics, steady and with vibrato, renders them through the engine, and reports

  amp_energy  : per-partial reconstructed energy / input   (1.0 = amplitude kept)
  shape_corr  : windowed best-lag waveform correlation vs input (1.0 = phase/shape kept)

Baselines (current engine, amp_mode 0, residual off):
  STEADY   amp~0.97  shape~1.00   (fine)
  VIBRATO  amp~0.70  shape~0.47   (broken: amplitude smeared, relative phase scrambled)

A real fix (chirp-aware / reassigned / pitch-synchronous phase, and a
vibrato-robust amplitude) must push the VIBRATO row toward 1.0/1.0 here BEFORE
spending an ear pass on the real voice. Energy-integrated amplitude (amp_mode 1)
raised amp to 0.93 on this rig but regressed on the real voice (noise pickup) —
so this rig is necessary but not sufficient; always confirm on the real Female too.

Usage:
    python3 battery/vibrato_rig.py                 # default engine (mode 0)
    python3 battery/vibrato_rig.py --args "0 1"    # extra CLI args after block (e.g. semi synth_mode)
    python3 battery/vibrato_rig.py --binary build/additive
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile, wave
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from metrics import waveform_correlation

SR = 48000
FREqS = [150 * k for k in range(1, 25)]   # 150..3600 Hz, well separated at 4096


def _write(path, x):
    x = x / max(1e-9, np.max(np.abs(x))) * 0.7
    with wave.open(path, "wb") as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(SR)
        w.writeframes((np.clip(x, -1, 1) * 32767).astype("<i2").tobytes())


def _read(path):
    w = wave.open(path, "rb")
    return np.frombuffer(w.readframes(w.getnframes()), "<i2").astype(float) / 32768.0


def make_signals(tmp):
    t = np.arange(3 * SR) / SR
    steady = sum((1.0 / len(FREqS)) * np.sin(2 * np.pi * f * t) for f in FREqS)
    vibf = 1 + 0.03 * np.sin(2 * np.pi * 5.5 * t)
    vib = sum((1.0 / len(FREqS)) * np.sin(2 * np.pi * np.cumsum(f * vibf) / SR) for f in FREqS)
    ps, pv = os.path.join(tmp, "steady.wav"), os.path.join(tmp, "vib.wav")
    _write(ps, steady); _write(pv, vib)
    return ps, pv


def amp_energy(inp, out):
    xi, xo = _read(inp), _read(out)
    def e(y, f, wid=25):
        Y = np.abs(np.fft.rfft(y * np.hanning(len(y)))); b = int(round(f * len(y) / SR))
        return np.sqrt((Y[max(0, b - wid):b + wid] ** 2).sum())
    return float(np.mean([e(xo, f) / e(xi, f) for f in FREqS if e(xi, f) > 0]))


def shape_corr(inp, out, win=2048, hop=1024):
    xi, xo = _read(inp), _read(out); n = min(len(xi), len(xo)); xi, xo = xi[:n], xo[:n]
    cs = []
    for s in range(win, n - 2 * win, hop):
        a = xi[s:s + win]
        if np.sqrt(np.mean(a ** 2)) < 1e-3:
            continue
        cs.append(waveform_correlation(xo[s:s + win], a, max_lag=64))
    return float(np.mean(cs)) if cs else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default="build/additive")
    ap.add_argument("--block", default="2048")
    ap.add_argument("--args", default="0 0 0.0 0 0 -1 0",
                    help="CLI args after <block>: semi synth_mode residScale unityDedup residMode residHp ampMode")
    a = ap.parse_args()
    extra = a.args.split()
    with tempfile.TemporaryDirectory() as tmp:
        ps, pv = make_signals(tmp)
        print(f"{'signal':<10}{'amp_energy':>12}{'shape_corr':>12}")
        for name, inp in [("STEADY", ps), ("VIBRATO", pv)]:
            out = os.path.join(tmp, name + "_out.wav")
            subprocess.run([a.binary, inp, out, a.block, *extra], check=True, stdout=subprocess.DEVNULL)
            print(f"{name:<10}{amp_energy(inp, out):>12.3f}{shape_corr(inp, out):>12.3f}")
        print("\n  target: VIBRATO -> ~1.0 / ~1.0 (steady is already there). "
              "Confirm any fix on the real voice too.")


if __name__ == "__main__":
    main()
