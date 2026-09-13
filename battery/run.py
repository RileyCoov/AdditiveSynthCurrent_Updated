#!/usr/bin/env python3
"""Standing regression battery for the additive resynthesis engine (Phase 0).

Builds the portable CLI, renders every corpus entry across every pitch
condition, scores the results with the metrics in metrics.py, and prints a
pass/fail scorecard. Optionally diffs against a saved baseline scorecard so a
fix for one file can't silently regress another (the guard-rail check).

This is the deliverable that ends the whack-a-mole: a fixed objective every
change must pass, reproducible on any machine (the engine is portable C++).

Usage:
    python3 battery/run.py --manifest battery/corpus.json
    python3 battery/run.py --manifest battery/corpus.json --out r7.json
    python3 battery/run.py --manifest battery/corpus.json --baseline r7.json

Manifest (JSON) shape: see corpus.example.json. Paths are resolved relative to
the manifest file's directory.

Metric applicability:
    saw_corr        unity (0 semi) only; needs a reference wav (default: input)
    cepstral_excess unity (0 semi) only; needs a reference wav (default: input)
    env_p2p         all conditions (steady-tone wobble; also the shift beating
                    proxy for doubling, where no clean same-pitch reference
                    exists)
    residual        unity (0 semi) only; needs a reference wav (default: input).
                    Emits residual_srr_db (+ per-band) and residual_flat_ratio.
                    Phase-sensitive and defined on EVERY class -- declare it on
                    every entry. See metrics.py for why the battery needed a
                    phase-sensitive metric outside the two saw entries.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, tempfile, wave
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from metrics import (waveform_correlation, cepstral_excess, envelope_p2p,
                     trajectory_jitter, residual_srr, residual_flatness, BAND_NAMES,
                     envelope_p2p_dev, trajectory_jitter_dev)


def read_wav(path: str) -> tuple[np.ndarray, int]:
    with wave.open(path, "rb") as w:
        sr = w.getframerate()
        n = w.getnframes()
        ch = w.getnchannels()
        raw = np.frombuffer(w.readframes(n), dtype="<i2").astype(np.float64) / 32768.0
    if ch > 1:
        raw = raw.reshape(-1, ch).mean(axis=1)
    return raw, sr


def build(repo_root: str) -> None:
    bdir = os.path.join(repo_root, "build")
    subprocess.run(["cmake", "-S", repo_root, "-B", bdir, "-DCMAKE_BUILD_TYPE=Release"],
                   check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["cmake", "--build", bdir], check=True, stdout=subprocess.DEVNULL)


def render(binary: str, inp: str, out: str, block: int, semi: int, mode: int,
           joint: int = 0) -> None:
    # positional CLI: in out block semi synth resid dedup residMode residHp amp
    #                 phase pitchSync joint
    cmd = [binary, inp, out, str(block), str(semi), str(mode)]
    if joint:
        cmd += ["1", "0", "0", "-1", "0", "0", "1", str(joint)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def check_unique_names(entries: list[dict]) -> None:
    """Entry names key the scorecard and the baseline diff, so duplicates silently
    make one file's measurement get compared against another file's baseline.
    (Both Fairlight entries were named 'FairlightAdditiveWaveTA1' until Sep 2026,
    which is exactly what happened.) Fail loudly instead."""
    seen = set()
    dupes = sorted({e["name"] for e in entries if e["name"] in seen or seen.add(e["name"])})
    if dupes:
        raise SystemExit(f"manifest error: duplicate entry names {dupes} -- "
                         "names must be unique, they key the baseline diff")


def score_entry(entry: dict, base: str, binary: str, block: int,
                conditions: list[int], tmp: str, mode: int, joint: int = 0) -> list[dict]:
    inp = os.path.join(base, entry["path"])
    ref_path = os.path.join(base, entry["reference"]) if entry.get("reference") else inp
    wants = set(entry["metrics"])
    thr = entry.get("thresholds", {})
    rows = []
    for semi in conditions:
        out = os.path.join(tmp, f"{entry['name']}_{semi:+d}.wav")
        render(binary, inp, out, block, semi, mode, joint)
        y, sr = read_wav(out)
        vals = {}
        if semi == 0 and "saw_corr" in wants:
            ref, _ = read_wav(ref_path)
            vals["saw_corr"] = waveform_correlation(y, ref)
        if semi == 0 and "cepstral_excess" in wants:
            ref, _ = read_wav(ref_path)
            vals["cepstral_excess"] = cepstral_excess(y, ref, sr)
        if "env_p2p" in wants:
            full, band = envelope_p2p(y, sr)
            vals["env_p2p_full"] = full
            vals["env_p2p_band"] = band
        if "jitter" in wants:
            vals["jitter_db"] = trajectory_jitter(y, sr)
        # Reference-relative forms. Unity only (they need the same-pitch input).
        # These are the ones to gate on for anything that is not a steady tone --
        # see envelope_p2p_dev in metrics.py.
        if semi == 0 and "env_p2p_dev" in wants:
            ref, _ = read_wav(ref_path)
            df, db_ = envelope_p2p_dev(y, ref, sr)
            vals["env_p2p_dev_full"] = df
            vals["env_p2p_dev_band"] = db_
        if semi == 0 and "jitter_dev" in wants:
            ref, _ = read_wav(ref_path)
            vals["jitter_dev_db"] = trajectory_jitter_dev(y, ref, sr)
        if semi == 0 and "residual" in wants:
            ref, _ = read_wav(ref_path)
            total, bands = residual_srr(y, ref, sr, per_band=True)
            vals["residual_srr_db"] = total
            for nm, bv in zip(BAND_NAMES, bands):
                if bv == bv:                      # skip NaN (empty band)
                    vals[f"residual_srr_{nm}"] = bv
            vals["residual_flat_ratio"] = residual_flatness(y, ref)
        for metric, value in vals.items():
            # full-envelope p2p is only meaningful for near-sinusoidal steady
            # tones; on harmonic/vibrato material the Hilbert envelope beats
            # between harmonics, so gate voice-class wobble on env_p2p_band.
            ok = True
            if metric == "saw_corr" and "saw_corr_min" in thr:
                ok = value >= thr["saw_corr_min"]
            elif metric == "cepstral_excess" and "cepstral_excess_max" in thr:
                ok = value <= thr["cepstral_excess_max"]
            elif metric == "env_p2p_full" and "env_p2p_full_max" in thr:
                ok = value <= thr["env_p2p_full_max"]
            elif metric == "env_p2p_band" and "env_p2p_band_max" in thr:
                ok = value <= thr["env_p2p_band_max"]
            elif metric == "jitter_db" and "jitter_max" in thr:
                ok = value <= thr["jitter_max"]
            elif metric == "residual_srr_db" and "residual_srr_min" in thr:
                ok = value >= thr["residual_srr_min"]
            rows.append({"entry": entry["name"], "class": entry.get("class", ""),
                         "semi": semi, "metric": metric, "value": value, "pass": ok})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--binary", default=None, help="override CLI path (default build/additive)")
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--out", default=None, help="write scorecard JSON here")
    ap.add_argument("--baseline", default=None, help="scorecard JSON to diff against")
    ap.add_argument("--regress-tol", type=float, default=0.02,
                    help="relative worsening vs baseline that counts as a regression")
    ap.add_argument("--regress-abs", type=float, default=0.02,
                    help="absolute slack added to the tolerance so near-zero metrics "
                         "don't false-alarm (e.g. steady 0.0007->0.0014 is +96%% but noise)")
    ap.add_argument("--synth-mode", type=int, default=0,
                    help="engine: 0 = OLA (default), 1 = oscillator bank (MQ)")
    ap.add_argument("--joint-mode", type=int, default=0,
                    help="amp/phase estimation: 0 = per-peak (default), "
                         "1 = joint LS in analysis, 2 = joint LS at final track freqs")
    args = ap.parse_args()

    base = os.path.dirname(os.path.abspath(args.manifest))
    repo_root = os.path.dirname(base)
    with open(args.manifest) as f:
        man = json.load(f)
    block = int(man.get("block_size", 2048))
    conditions = [int(s) for s in man.get("conditions_semitones", [-5, 0, 5])]
    binary = args.binary or os.path.join(repo_root, man.get("binary", "build/additive"))

    if not args.no_build:
        build(repo_root)
    if not os.path.exists(binary):
        print(f"ERROR: binary not found: {binary}", file=sys.stderr)
        return 2

    check_unique_names(man["entries"])

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for entry in man["entries"]:
            rows.extend(score_entry(entry, base, binary, block, conditions, tmp,
                                    args.synth_mode, args.joint_mode))

    # scorecard
    print(f"\n{'entry':<16}{'class':<10}{'semi':>5}  {'metric':<18}{'value':>12}  result")
    print("-" * 74)
    n_fail = 0
    for r in rows:
        flag = "PASS" if r["pass"] else "FAIL"
        if not r["pass"]:
            n_fail += 1
        print(f"{r['entry']:<16}{r['class']:<10}{r['semi']:>+5d}  "
              f"{r['metric']:<18}{r['value']:>12.5f}  {flag}")

    # guard-rail regression check vs baseline
    n_regress = 0
    if args.baseline:
        with open(args.baseline) as f:
            bl = {(x["entry"], x["semi"], x["metric"]): x["value"] for x in json.load(f)}
        # higher-is-better for saw_corr; lower-is-better for the rest
        higher_better = {"saw_corr"} | {m for m in
                         (r["metric"] for r in rows)
                         if m.startswith("residual_srr") or m == "residual_flat_ratio"}
        print("\nGuard-rail diff vs baseline:")
        for r in rows:
            key = (r["entry"], r["semi"], r["metric"])
            if key not in bl:
                continue
            old, new = bl[key], r["value"]
            # scale the relative tolerance by |old| so the test behaves correctly for
            # zero and negative baselines (a dB value can legitimately be negative;
            # old*(1-tol) moves the threshold the WRONG way when old < 0).
            slack = abs(old) * args.regress_tol + args.regress_abs
            worse = (new < old - slack) if r["metric"] in higher_better \
                else (new > old + slack)
            if worse:
                n_regress += 1
                print(f"  REGRESS {r['entry']:<14} {r['semi']:+d} {r['metric']:<16} "
                      f"{old:.5f} -> {new:.5f}")
        if n_regress == 0:
            print("  none")

    if args.out:
        with open(args.out, "w") as f:
            json.dump(rows, f, indent=2)
        print(f"\nwrote scorecard: {args.out}")

    print(f"\n{len(rows)} metrics | {n_fail} threshold failures | {n_regress} regressions")
    return 1 if (n_fail or n_regress) else 0


if __name__ == "__main__":
    sys.exit(main())
