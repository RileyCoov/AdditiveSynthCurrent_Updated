# Stage-0 ceiling probes

One-off diagnostic probes, not part of the standing battery. They answer "how much
headroom is left, and is it estimation or model structure?" — the question that decides
what to work on next. Results and interpretation: `docs/engine-direction-2026-09.md` §4b.

Run from the repo root with the engine built (`cmake -S . -B build && cmake --build build`).
They read the 15-file corpus from `~/Desktop/comparingOut/inputData` and renders from
`/tmp/stage0/render{,_nr}`; regenerate those with:

```sh
mkdir -p /tmp/stage0/render /tmp/stage0/render_nr
for f in ~/Desktop/comparingOut/inputData/*.wav; do b=$(basename "$f" .wav)
  ./build/additive "$f" /tmp/stage0/render/${b}_orig.wav    2048 0      # as shipped
  ./build/additive "$f" /tmp/stage0/render_nr/${b}_orig.wav 2048 0 0 0  # residual off
done
```

| script | what it measures |
|---|---|
| `lib.py` | shared wav IO, delay alignment, SRR, per-band SRR, flatness |
| `s01_residual_meter.py` | per-class residual meter + how much of the unity output is the pasted original |
| `s02_ceiling.py` | Ceiling A: best K-sinusoid per-frame fit vs K (the model-class upper bound) |
| `s03_trackceiling.py` | Ceiling B attempt (phase-continuous MQ via Adam). **Does not converge** — kept as the recorded negative result, do not read its numbers as a ceiling |
| `s04_freq_vs_amp.py` | splits the engine-to-ceiling gap into frequency error vs amplitude/phase error; needs `PFREQ_DEBUG` |

`s02` and `s04` are self-validating: `s02`'s probe scores a pure sine at 98 dB and white
noise at 1.1 dB (K=32), and `s04` has a null control (a 512-partial basis at random
frequencies fits white noise to only 4.8 dB). Re-run those sanity checks if you change them.

Engine instrumentation they depend on, both env-gated and inert when unset:
`PCOUNT_DEBUG` (partials rendered per frame) and `PFREQ_DEBUG` (per-frame track
frequencies + amplitudes), in the synthesis loop in `main.cpp`.
