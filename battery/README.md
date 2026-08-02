# Regression battery (Phase 0)

A standing, reproducible objective for the additive resynthesis engine. Every
change to the DSP must pass this before it's judged by ear — it turns perceptual
complaints ("double voice", "shape ripple", "volume wobble") into fixed numbers
so a fix for one file can't silently regress another.

The engine is portable C++, so this runs on any machine (Linux/CI or the Mac),
not just where the reference renders live.

## Run

```sh
# builds build/additive via CMake, renders the corpus x pitch conditions, scores it
python3 battery/run.py --manifest battery/corpus.json

# save a baseline scorecard (e.g. the current r7 behavior)
python3 battery/run.py --manifest battery/corpus.json --out r7.json

# later: fail if any metric regresses vs the baseline (the guard-rail check)
python3 battery/run.py --manifest battery/corpus.json --baseline r7.json

# Phase-1 A/B: render everything with the oscillator bank (MQ) and diff vs the
# OLA baseline r7.json -> the doubling/wobble metrics should improve, guards hold
python3 battery/run.py --manifest battery/corpus.json --synth-mode 1 --baseline r7.json
```

Exit code is non-zero if any threshold fails or any metric regresses, so it
drops straight into CI. Requires `numpy` and `scipy`. `--synth-mode` selects the
engine (0 = OLA, 1 = oscillator bank); the regression check has both a relative
tolerance (`--regress-tol`) and an absolute slack (`--regress-abs`) so near-zero
metrics don't false-alarm.

## Corpus

Copy `corpus.example.json` to `corpus.json` and point `path` at your WAVs
(relative to the manifest). Commit a small representative subset (a few seconds
each) so the battery is reproducible in CI; keep the full ~15-source corpus on
the Mac. Each entry declares its material `class`, the `metrics` to compute, and
per-entry `thresholds`.

## Metrics (`metrics.py`)

| Metric | What it catches | Notes |
|---|---|---|
| `saw_corr` (`waveform_correlation`) | Saw shape ripple / OLA seam | Max normalized cross-correlation vs a reference over +/- lag (removes constant delay). Unity only. Target >= 0.999. **Note:** at G0 both engines sat ~0.990 on real saws, i.e. this is analysis-limited, not OLA-limited — see below. |
| `env_p2p_full` / `env_p2p_band` | Steady-tone volume wobble; shift beating | Analytic-envelope peak-to-peak, full and band-limited to 0.3-1 Hz. All conditions. Gate steady tones on `env_p2p_full_max`; gate voice-class slow wobble on `env_p2p_band_max`. |

`cepstral_excess(rendered, reference)` also exists in `metrics.py` but is **not a
gate** — see finding 1.

### Findings from validating the metrics against the real corpus (not assumed)

1. **The Female "double voice" is not offline-measurable before a fix exists.**
   Three cepstral approaches were tried and each failed on real material: an
   absolute in-band peak (a voice's own pitch rahmonics fill the 5-40 ms band), a
   reference-relative excess (`cepstral_excess` — its excess sat at the pitch
   quefrencies 5-16 ms, not the 21.3 ms OLA hop, and scored good-sounding files
   *higher* than the doubled voice), and a hop-quefrency prominence (Female ~ a
   pure-tone control). Root cause: the doubling is a *shift-mode* artifact and
   measuring it needs a doubling-free reference at the shifted pitch — which does
   not exist until the oscillator bank produces one. So the Female gate is the
   **Phase-1 A/B**: render the same shifted file OLA vs MQ (`--synth-mode 1`) and
   confirm the doubling drops. (`env_p2p_band` on voice tracks slow wobble but is
   confounded with musical dynamics, so it is tracked, not gated.)
2. **`env_p2p_full` is only meaningful for near-sinusoidal steady tones.** On
   harmonic/vibrato material the Hilbert envelope beats between harmonics.

## Gate G0 (as run)

- **Saw** (`saw_corr`) and **steady** (`env_p2p_full`) are validated absolute
  gates. On the real corpus both saws failed 0.999 (~0.990) and steady was clean.
- **Female doubling** has no pre-fix gate (finding 1); validated by Phase-1 A/B.
- All other files are **tracked** (no hard gate) for the regression diff.
- `saw_corr ~0.990 in both engines` means the saw shape error is **analysis-
  limited**, not an OLA overlap artifact — so MQ is expected to help the doubling/
  wobble (frequency-moving cases) far more than the steady saw. If the saw still
  rings after MQ, the fix is in analysis, not synthesis.

## Engine A/B (Phase 1)

`--synth-mode` selects the engine (0 = OLA, 1 = oscillator bank / MQ). The MQ path
renders one continuous oscillator per track instead of overlap-adding per-frame
constant-frequency copies, so the OLA comb (saw ripple, Female doubling, wobble)
has no overlap to form. Lock the OLA baseline (`--out r7.json`), then diff MQ
against it (`--synth-mode 1 --baseline r7.json`).

## Pitch conditions

`conditions_semitones` (default `[-5, 0, 5]`) and the engine mode are passed via
CLI (`additive in.wav out.wav <block> <semitones> <synth_mode>`), so all
conditions and both engines render from one binary with no recompile.
