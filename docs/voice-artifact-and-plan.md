# Voice "phasey double-voice" — findings and next plan

Consolidation of the diagnostic session. The goal was the Female "double voice";
this records what was built, what was proven, and the evidence-driven plan for the
actual root. Written so a future session (fresh clone) can resume without re-deriving.

## Project orientation (read first)

**What it is.** A true additive analysis/resynthesis engine: it decomposes audio
into sinusoidal partials (per-partial frequency, amplitude, phase) that can be
modified (pitch shift, etc.) and resynthesized. Pure portable C++ in
`AdditiveSynthClean/`; ships as an Xcode project but also builds headless.

**Hard constraint.** It must stay *decomposable* — every output sample comes from
summed sinusoids (+ a stochastic/true residual). No passthrough of the original
(that is why hp300 true-phase residual, though it sounded great, is unusable: it
overlays the original signal above 300 Hz).

**Build & run (headless, any OS):**
```
cmake -S . -B build && cmake --build build
./build/additive <in.wav> <out.wav> <block> [semi] [synth_mode] [residScale] \
                 [unityDedup] [residMode] [residHp] [ampMode]
```
CLI args (positional; defaults preserve legacy behavior):
1-3 in / out / block · 4 `semi` pitch shift · 5 `synth_mode` 0=OLA 1=MQ · 6
`residScale` ×residual gain (0=off) · 7 `unityDedup` · 8 `residMode` 0=stochastic
1=true-phase(unity) · 9 `residHp` residual HP override · 10 `ampMode` 0=parabolic
1=energy. (Battery/rig use only up to arg 5; args 6-10 are diagnostics/A-B.)

**Signal path & key code (all `AdditiveSynthClean/main.cpp` unless noted):**
- Analysis, two tiers: main FFT `ANALYSIS_SIZE=4096` (~L57); LF tier
  `LF_ANALYSIS_SIZE=16384` below `lf_cutoff_hz=200`. Peak pick → `detect_peaks` →
  `parabolic_interpolation` (~L338, amplitude=peak height, freq=parabolic) →
  phase read w/ fractional-bin correction (5ae358e) → PV instantaneous-frequency.
- Tracking: `find_best_match_peak_HZ`, per-track median-of-3 freq de-jitter (~L1052),
  amplitude EMA. **Phase is NOT smoothed** (re-derived per frame at unity).
- Synthesis: OLA per-frame constant-freq render + Hann overlap-add (default), or
  continuous-phase oscillator bank `mq_synthesize` (~L100, `synth_mode=1`).
- Residual (~L1852): stochastic spectral-deficit noise fill above `residual_hp_hz`
  =2500 (`residMode=0`), or true-phase `(original−model)` (`residMode=1`, unity).
- Output: COLA-hole bridge, windowed limiter, PCM16 WAV out.

**Tooling (committed, `battery/`):** `run.py` (corpus battery, scorecard +
baseline diff, `--synth-mode`), `metrics.py` (waveform_correlation,
trajectory_jitter, envelope_p2p, cepstral_excess), **`vibrato_rig.py` (the gate for
the voice fix — see below)**. Requires numpy+scipy.

**Working method that finally worked:** controlled synthetic rigs with a metric
that tracks the artifact, validated in-cloud, THEN confirmed by ear on the real
Female. Ear is the final authority — three metrics (cepstral, jitter) gave numbers
that disagreed with the ear; only the rig's shape-corr tracks this artifact.

## Assets built this session (durable)

- **Portable headless build** — `CMakeLists.txt` + `<cstdint>`/`<cmath>` fixes. The
  engine compiles and runs on Linux/CI, not just Xcode. `cmake -S . -B build && cmake --build build`.
- **Standing regression battery** — `battery/` (`metrics.py`, `run.py`, manifest).
  Renders corpus × pitch conditions × engine, scorecard + baseline-diff, CI-ready.
  Locked OLA baseline = `r7.json`.
- **Oscillator-bank (MQ) synthesis** behind `synth_mode` — `mq_synthesize()` in
  `main.cpp`. `synth_mode=0` (default) is byte-identical to the legacy OLA engine.
- **Diagnostic CLI args** — `additive in out <block> [semitones] [synth_mode]
  [residual_scale] [unity_dedup]`. Defaults preserve behavior.

## What was proven about the artifact (by bisection)

The Female "double voice" is a **phasey / hollow / warbly, not-quite-one-person**
quality — same pitch, a very short (<10 ms) smear, "like an effect on the voice."
Present at unity, worse under pitch shift.

| Hypothesis | Test | Result |
|---|---|---|
| OLA overlap comb (the plan's #1) | MQ oscillator bank (`synth_mode=1`) | **Not it.** No audible change; comb metric barely moved. Hypothesis falsified. |
| Pitch-shift path | render at unity (0 semi) | **Not it.** Artifact present at unity. |
| Duplicate/close partials | unity dedup (`unity_dedup=1`) | **Not it.** No audible change. |
| Stochastic residual | residual off (`residual_scale=0`) | **Minor.** Small, real improvement — keep residual gain lower, not the cause. |
| Vibrato smear (window too long) | `ANALYSIS_SIZE 4096 → 2048` | **Opposite:** shorter window = *faster, worse* warble. |

**Conclusion.** Shorter window → worse ⇒ the artifact is **per-frame analysis
measurement JITTER** in each partial's estimate, reconstructed audibly — not vibrato
smear, not synthesis, not shift, not dedup. It lives in the analysis core.

Key corroborating code fact: at unity the tonal reconstruction already de-jitters
**frequency** (median-of-3, `main.cpp` ~L1052) and smooths **amplitude** (EMA), but
**phase is re-derived raw every frame** (`phase0 = wrap(phase + inc*offset)`), then
overlap-added. Phase is the one un-smoothed per-partial parameter, so **phase jitter
is the prime suspect**. This also explains why MQ didn't help: MQ-unity matches the
same raw per-frame phases, inheriting the same jitter.

Constraint (from the session): window length is a **trade** — longer averages jitter
but smears transients (drum loop). So a global longer window is rejected; the fix must
cut jitter **without** costing temporal resolution.

Side conclusions:
- **MQ is a lateral move**, not an improvement — keep it behind `synth_mode`, default 0.
- **Saw shape is analysis-limited** too (~0.990 in both engines), not an OLA artifact.

## New plan — analysis-core jitter reduction

### Phase A — make jitter measurable (do first)
Add a per-partial **trajectory-jitter metric** to the battery. For a sustained voiced
region: fit/remove each track's smooth trend (freq vibrato ramp, amp envelope), then
measure frame-to-frame residual variance of freq, amplitude, and **phase** (phase
residual after removing the expected `∫2πf`). Unlike the doubling metric (which three
attempts could not isolate), jitter is directly measurable per track — this finally
gives a voice-relevant gate. Confirm it is high on Female, low on Fairlight/steady.

### Phase B — reduce jitter, cheapest first, each A/B'd by ear + Phase-A metric
1. **Smoothed phase reconstruction (highest-signal, phase is the un-smoothed one).**
   At unity, stop re-deriving raw phase per frame. Reconstruct phase by integrating
   the (already de-jittered) frequency trajectory, slowly re-anchoring to measured
   phase — a leaky phase-locked integrator: smooth like shift-mode propagation, but
   drift-corrected so waveform shape isn't lost. This is the MQ `match_phase=false`
   path plus a slow phase-correction term; the code is mostly present.
2. **Extend de-jitter to amplitude/phase, not just frequency.** The median-of-3 is
   freq-only; apply an analogous spike-robust filter to per-track amplitude and phase
   residual below `traj_median_max_hz`.
3. **Better per-frame estimator** — reduce jitter at the source (reassigned-spectrum
   frequency/phase; improve the parabolic phase read).
4. **Adaptive analysis window** — long window for steady/voiced regions, short for
   transients, keyed off the existing transient detector. This is the "real" fix if
   smoothing plateaus, and it directly honors the transient constraint above. Higher
   effort.

### Phase C — validate & gate
Phase-A jitter metric drops on voice; guard rails (drums/cymbals/Fairlight/steady)
hold in the battery; ear confirms less phasiness. Gate every change on the battery.

## DEFINITIVE ROOT CAUSE (controlled experiment)

Later work localized it precisely. Feeding the engine **known equal-amplitude
partials** (pure tonal model, residual off) and measuring reconstructed band
energy:

- **Steady partials:** energy out/in ≈ **0.97** (on-bin partials 1.00, off-bin
  ~0.95). The tonal model is accurate on steady tones. OLA and MQ synth are
  identical here — the loss is not in synthesis.
- **Same partials WITH 3% / 5.5 Hz vibrato:** energy out/in collapses to **0.70
  mean**, wildly uneven per partial (0.375 … 1.12).

**Conclusion:** the artifact is the **analysis** (4096-pt / 85 ms window) failing
to measure partials that are *moving* under vibrato — they smear across bins, the
peak-picker mis-reads each partial's amplitude, and resynthesis rebuilds the voice
with wrong, fluctuating amplitudes → the phasey/warbly/"not-one-person" timbre.
This is the Gabor limit (bottleneck #3), triggered by vibrato. It explains why no
window length helped (longer smears more; shorter is noisier) and why hp300
true-phase residual sounded right (it overlays the real signal, bypassing the
mis-measured tonal model — but it is near-passthrough, not true additive synthesis,
so it is unusable for a decomposable engine).

The Female "double voice" under pitch shift is the same root, amplified: the
mis-measured moving-partial amplitudes/phases are then propagated through the
shift path.

### Revised fix plan — vibrato-aware analysis (the real lever)
Steady partials already reconstruct at ~1.0, so the target is measuring MOVING
partials correctly. Options, cheapest first, each A/B'd on the Mac + the
`jitter`/energy metrics:
1. **Energy-integrated amplitude estimate.** Replace the parabolic peak-height
   amplitude with the integrated main-lobe energy around each tracked partial, so
   a partial smeared by vibrato keeps its true amplitude. Bounded change to the
   analysis; prototype-test on the vibrato equal-harmonics rig (does 0.70 → ~1.0?).
   Risk: energy bleed between close partials — guard by the partial's resolved
   bandwidth.
2. **Chirp-aware / reassigned estimation.** Estimate each partial's frequency
   *slope* within the window and correct the amplitude/phase for the smearing
   (spectral reassignment). More accurate, more work.
3. **Vibrato-demodulated (pitch-synchronous) analysis.** Track f0, time-warp to
   remove the vibrato, analyze the now-steady partials, un-warp. The classic
   high-quality monophonic-voice approach; largest effort.

Gate on the controlled vibrato rig (energy out/in → 1.0) + the `jitter` metric
(rendered → original) + the ear. Keep it fully additive (no passthrough).

## Energy-integrated amplitude (amp_mode=1): TRIED, FAILED on real voice

Implemented energy-integrated amplitude (sum main-lobe power over ±6 bins →
equivalent-sinusoid amplitude) to counter the vibrato peak under-read. On the
clean synthetic rig it worked (round-trip 0.70 → 0.93 vibrato, 0.999 steady),
but on the **real** Female it sounded clearly worse — jagged waveform, amplitude
now fluctuating frame-to-frame. Cause: a real voice's spectrum has noise/breath/
partial skirts, and the ±bin energy sum scoops that up, so per-partial amplitude
jitters with the noise. The peak estimator (amp_mode=0) ignores that surround and
stays smooth. **amp_mode=1 overfit the rig; it is a dead end on real material.**
It remains opt-in and amp_mode=0 is byte-identical, so nothing regressed. Also
notable: the `jitter` metric DROPPED for the worse-sounding render — it does not
capture this artifact either.

### Sharpened hypothesis: it's PHASE under vibrato, not amplitude
Ruled out on the real voice now: OLA comb (MQ), shift path (present at unity),
duplicate partials, stochastic residual (minor), window length (both worse),
amplitude estimation (worse). The only thing that sounded correct was hp300 true-
phase residual — which substitutes the ORIGINAL's phase above 300 Hz. Every
amplitude-side fix failed. This points to **per-partial phase accuracy under
vibrato**: a moving partial's phase, measured once per 85 ms frame, is wrong, so
the harmonics' relative phases (waveform shape) are off → "phasey/doubled, same
pitch," worse under shift (which integrates phase error). MQ didn't help because
it matched those same wrong measured phases.

Next real lever (a dedicated effort, not a quick knob): recover accurate
per-partial phase/frequency for moving partials — spectral reassignment, or
pitch-synchronous / vibrato-demodulated analysis — validated by the EAR (metrics
have repeatedly failed to capture this artifact) on the real voice, kept fully
additive. Consider a phase-accuracy rig: resynthesize a KNOWN vibrato harmonic
signal and compare per-partial reconstructed phase to ground truth.

## CONFIRMED: phase (not amplitude) — chirp-biased analysis under vibrato

Built the vibrato rig (battery/vibrato_rig.py). Reconstruction fidelity vs input:
- STEADY:  amp 0.97, waveform-shape corr 1.00.
- VIBRATO: amp 0.70, waveform-shape corr **0.47**.

Shape corr 0.47 (< the 0.70 amplitude alone would give) means the harmonics'
RELATIVE PHASES are scrambled under vibrato = the "phasey/double voice." MQ was
slightly WORSE (0.44), so it is not the synthesis: both engines faithfully render
wrong MEASURED phases. Under vibrato each partial is a chirp within the 85 ms
window; the phase estimate assumes a stationary sinusoid (the 5ae358e fix only
corrects the fractional-bin offset, not the chirp) -> biased phase.

**The fix is chirp-aware / reassigned phase estimation** (estimate each partial's
frequency slope within the window and correct its phase), developed against the
rig's VIBRATO shape corr (0.47 -> ~1.0) and then confirmed on the real Female by
ear. Amplitude is secondary and must not use raw energy integration (noise pickup).

## NEXT: chirp-aware phase estimation — implementation spec

Goal: make the analysis measure the phase (and frequency) of a partial that is a
linear **chirp** within the window, so the harmonics' relative phases are right
under vibrato. Gate: `python3 battery/vibrato_rig.py` VIBRATO `shape_corr`
0.47 → ~1.0 (steady 1.00 must not regress), THEN confirm on the real Female by ear.

Method — **spectral reassignment** (Auger–Flandrin), applied in the main analysis
tier where `analysis_frame`/`analysis_mag`/`analysis_phase_store` are built:
1. Precompute two extra length-`ANALYSIS_SIZE` windows once:
   - `th[n] = ((n − (N−1)/2)) * hann[n]` (time-weighted, in samples)
   - `dh[n]` = discrete derivative of `hann` (e.g. central difference).
2. Per frame, FFT the same samples windowed by `th` and by `dh` → complex `Xth`,
   `Xdh` (two more RealFFTs; reuse the RealFFT packing re=idx k, im=idx N−k).
3. For each peak bin k (complex `X[k]` from the existing Hann FFT), with
   `p = X[k]`, `d = |X[k]|^2`:
   - reassigned freq (rad/sample): `w_hat = w_k − Im(Xdh[k]·conj(p)) / d`
   - reassigned time (samples): `t_hat = Re(Xth[k]·conj(p)) / d`
   - corrected phase at frame center: `phi_hat = arg(p) + w_hat * t_hat`
   Use `w_hat` as the frequency and `phi_hat` as `peak.phase` (replacing the
   current stationary estimate). Amplitude: keep the parabolic peak (amp_mode=0);
   do NOT use energy integration (noise pickup — proven bad on real voice).
4. Sign/convention caveat: exact signs depend on the FFT/derivative conventions.
   The rig is self-correcting — implement, run it; if VIBRATO shape_corr rises
   toward 1.0, conventions are right; if it drops, flip the sign on `t_hat`/`w_hat`.
5. Gate behind a CLI arg (e.g. arg 11 `phase_mode`, default 0 = current) so it is
   A/B-able and byte-identical off, like every other change here.

If reassignment plateaus below ~1.0 on the rig, escalate to pitch-synchronous /
vibrato-demodulated analysis (track f0, time-warp out the vibrato, analyze steady,
un-warp) — larger effort, highest quality for monophonic voice.

Also verify the fix helps **pitch shift** (the "double voice" is worse there): the
corrected phase/freq feed the shift path, so it should carry.

## Housekeeping (do before next session)
- **Revert `ANALYSIS_SIZE` to 4096** (2048 was strictly worse). — DONE (already 4096).
- **Lower `residual_noise_gain`** (was 1.3) — e.g. ~0.6 — to bank the residual-off win;
  consider making it per-material later. — DEFERRED: it perturbs phase_mode=0 output, so
  keep it out of the reassignment A/B; do it as its own eared change.
- Keep `synth_mode` default 0. MQ stays available for A/B.

## RESULT: first-order spectral reassignment does NOT fix it (Aug 2)

Implemented reassignment behind a new opt-in CLI arg 11 `phaseMode` (bitmask: 1=phase,
2=frequency, 3=both; default 0 = byte-identical to before). Built headless, gated on
`vibrato_rig.py`. Verdict — **it does not work; frequency was never the problem:**

| phase_mode | VIBRATO shape_corr | vs baseline 0.469 |
|---|---|---|
| 0 (off) | 0.469 | — |
| 1 (reassigned phase) | 0.206 | worse |
| 2 (reassigned freq)  | 0.419 | ~no-op (slightly worse) |
| 3 (both) | 0.208 | worse |

`RA_DEBUG` dump (VIBRATO, frame 30, true partials 150·k Hz) shows **`w_hat` matches the
existing parabolic/PV estimate to ~1 Hz** — the analysis already measures frequency well,
so the reassigned-frequency path is a no-op. The reassigned **time** `t_hat` is large and
its sign varies per harmonic (−40, +19, −84, +70, −88, +137, −63, +193 samples across the
first 8 partials), so the phase correction `phi + w_hat·t_hat` is a **multi-radian
per-partial term that scrambles relative phase in either sign** (both signs → ~0.21). This
is why it degrades instead of helping.

**Root reason the linear model fails here:** an 85 ms (4096-pt) analysis window spans
~half a 5.5 Hz vibrato cycle (0.085·5.5 ≈ 0.47 cycle), so within one window each partial
*curves* rather than chirps linearly — a first-order (and even a second-order linear-chirp)
reassignment model is a poor fit. The per-partial `t_hat` spread confirms every harmonic
sits at a different point in its sweep inside the window; a single per-frame (freq, amp,
phase) snapshot cannot represent that. This is the Gabor limit, and reassignment does not
lift it for a window this long relative to the modulation.

The `phaseMode` code + `RA_DEBUG` harness are kept opt-in (mirroring `amp_mode=1`'s
"validated dead end, kept for A/B") and are byte-identical when off. Not committed.

### Redirected recommendation → pitch-synchronous / vibrato-demodulated analysis
The remaining lever is to make the partials **stationary before analysis**: track f0,
time-warp the vibrato out, analyze the now-steady spectrum (where the existing stationary
phase read already gives shape_corr ~1.0, per the STEADY row), then un-warp for synthesis.
This is the classical high-quality monophonic-voice approach and the plan's stated
escalation — a substantial, separate build. Gate remains `vibrato_rig.py` VIBRATO row → 1.0,
then ear on the real Female.
