#!/usr/bin/env python3
"""Transient rig — one hit at t=0.75 s in digital silence, so pre-hit energy is pre-echo.

Busy percussion makes every envelope metric ambiguous: "before this hit" is also "after
the previous hit", and a slow amplitude release reads as pre-echo. That ambiguity sent
one round of this project's transient analysis down a wrong path (docs 4d.7). Here the
input is exactly zero before the hit, so any output energy there is the engine's own.

Usage:
    python3 battery/transient_rig.py --gen  OUTDIR             # write the 4 rig wavs
    python3 battery/transient_rig.py --score OUTDIR BIN        # render + score
"""
import argparse, os, subprocess, sys, wave
import numpy as np

SR = 48000
T_HIT = 0.75
FILES = ['kick', 'snare', 'click', 'note']


def _write(path, x):
    x = np.clip(x, -1, 1)
    with wave.open(path, 'wb') as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(SR)
        w.writeframes((x * 32767).astype('<i2').tobytes())


def _read(path):
    with wave.open(path, 'rb') as w:
        n, ch = w.getnframes(), w.getnchannels()
        x = np.frombuffer(w.readframes(n), dtype='<i2').astype(float) / 32768.0
    return x.reshape(-1, ch).mean(1) if ch > 1 else x


def generate(outdir):
    os.makedirs(outdir, exist_ok=True)
    n = int(1.5 * SR); t = np.arange(n) / SR
    rng = np.random.default_rng(0)
    def decay(tau):
        e = np.zeros(n); m = t >= T_HIT
        e[m] = np.exp(-(t[m] - T_HIT) / tau); return e
    dt = np.maximum(t - T_HIT, 0)
    kick = np.sin(2 * np.pi * (60 + 40 * np.exp(-dt / 0.02)) * dt) * decay(0.12)
    snare = 0.7 * rng.standard_normal(n) * decay(0.08) + 0.5 * np.sin(2 * np.pi * 200 * dt) * decay(0.15)
    click = rng.standard_normal(n) * decay(0.002)
    note = np.sin(2 * np.pi * 440 * dt) * decay(0.35)
    _write(f'{outdir}/rig_kick.wav', 0.8 * kick)
    _write(f'{outdir}/rig_snare.wav', 0.6 * snare / np.abs(snare).max())
    _write(f'{outdir}/rig_click.wav', 0.8 * click / np.abs(click).max())
    _write(f'{outdir}/rig_note.wav', 0.7 * note)
    return [f'{outdir}/rig_{f}.wav' for f in FILES]


def score(outdir, binary, shifts=(0, 7, -7)):
    generate(outdir)
    print(f"{'file':8s} {'shift':>6s} {'pre-echo':>10s} {'onset delay':>12s}")
    print("(pre-echo = peak output level before the hit, dB below the hit's own peak;")
    print(" the input is digital silence there, so lower is better. unity should be ~-240.)")
    rows = []
    for f in FILES:
        src = f'{outdir}/rig_{f}.wav'
        for sh in shifts:
            dst = f'{outdir}/rig_{f}_s{sh}.wav'
            subprocess.run([binary, src, dst, '2048', str(sh)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            y = _read(dst)
            i0 = int(T_HIT * SR); pk = np.abs(y).max()
            pre = np.abs(y[:i0 - int(0.005 * SR)]).max()
            pre_db = 20 * np.log10(pre / pk + 1e-12)
            delay = (int(np.argmax(np.abs(y) > 0.1 * pk)) - i0) / SR * 1000
            print(f"{f:8s} {sh:>6d} {pre_db:9.1f} {delay:+11.1f} ms")
            rows.append(dict(file=f, shift=sh, pre_echo_db=pre_db, onset_delay_ms=delay))
    return rows


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--gen', metavar='OUTDIR')
    ap.add_argument('--score', nargs=2, metavar=('OUTDIR', 'BIN'))
    a = ap.parse_args()
    if a.gen:
        print('\n'.join(generate(a.gen)))
    elif a.score:
        score(a.score[0], a.score[1])
    else:
        ap.print_help(); sys.exit(1)
