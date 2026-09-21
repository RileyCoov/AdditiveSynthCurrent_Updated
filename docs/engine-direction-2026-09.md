# Engine direction — evidence review and plan (Sep 2026)

Written in response to the attached research report (`compass_artifact_wf-eb2aa66d…md`)
and the research brief behind it (Q14/Q15: "can we build a structured, synthesis-oriented
latent representation that beats hand-designed FFT analysis?").

**Nothing in the engine has been changed.** This document is the analysis and the plan.

---

## 0. Verdict in four sentences

1. The attached report answers a **narrower question** than the one asked, and two of its
   five concrete DSP recommendations are things this codebase has already tried and
   rejected by ear. Use it as a literature / licensing / evaluation reference, not as the
   work plan.
2. The central research question ("neural latent between DSP and embeddings") is
   **technically sound in general and the wrong move for this engine right now** — and I
   can show that with measurements from this repo's own outputs, not from the literature.
3. The representation the question describes — *interpretable harmonic + transient + noise
   components, multi-resolution, synthesis-oriented, editable* — is achievable here
   **classically**, and is in fact the single highest-value piece of work remaining
   (Stage 2 below). The engine already has ~70% of it.
4. Before committing to any of it, run **one cheap experiment (Stage 0.2)** that decides
   whether the remaining gap is an *estimation* problem (fixable by better analysis,
   learned or not) or a *model-structure* problem (only fixable by adding components).
   That experiment needs no training data, no license audit, and no C++ changes.

---

## 1. What the attached report actually answers

The report is a good piece of work on the question *"should we replace the classical
spectral-peak validity test with a learned classifier?"* Its answer is "no, build the
classical ensemble first," and that answer is well-evidenced.

**Take these findings as load-bearing — they are correct and directly relevant:**

| Finding | Why it matters here |
|---|---|
| **Circularity**: never train on your own engine's resynthesis | This repo has already lived this failure twice — see §2.3. This is the most important single sentence in the report. |
| **Residual meter is the highest-value / lowest-risk item** | Agreed, and it is Stage 0.1 below. `e = x − x̂` is measured on the real input, so it cannot lie about representability. |
| **Evaluation = per-class regression CI + blind A/B, not absolute quality** | This is exactly what `battery/` already is. The report validates the existing approach. |
| **Licence audit** (GuitarSet / NSynth / Slakh2100 / MusicNet-orig / Vocadito are CC-BY; MDB-stem-synth, MedleyDB, MAESTRO, MERT are NC) | Only matters at Stage 3, but it's the correct list and worth keeping. |
| **Deployment**: a sub-1M-param model needs no ONNX Runtime; hand-rolled inference is practical | Correct, and it removes the entire binary-size / ODR / notarisation risk if Stage 3 ever happens. |

**Discount these — they are literature-derived and locally falsified:**

| Report recommendation | What this repo already found |
|---|---|
| §8(a) "vibrato double-voice → use reassignment / derivative method" | **Tried and it made things worse.** `phaseMode` arg 11, Aug 2: vibrato rig `shape_corr` 0.469 → **0.206**. Root reason documented in `docs/voice-artifact-and-plan.md`: an 85 ms window spans ~0.47 of a 5.5 Hz vibrato cycle, so each partial *curves* rather than chirps linearly; a first-order model is a bad fit. `RA_DEBUG` showed `w_hat` already agreed with the existing estimate to ~1 Hz — frequency was never the problem. |
| §8(e) "brightness under shift → Röbel–Rodet true-envelope formant-preserving shift" | **Prototyped and rejected by ear** (Jul, `jul3Work_formant/`): metrics looked great (Female up5 formant centroid 4583 vs input 4592 vs 6086 un-preserved), Riley: *"100% digresses and makes it worse."* Röbel–Rodet's cepstral estimator is genuinely better than the sparse-partial envelope that was tried, so this is a *warning*, not a proof — but it is not a free win. |
| §8(c) "close bass partials → ESPRIT/MUSIC on a narrow LF band" | Largely already solved classically: `LF_ANALYSIS_SIZE = 16384` tier below 200 Hz (`main.cpp:88`), with a stationarity gate, LF coast path and mutual-nearest assignment (commit `510108e`). SaintSaëns pedal partials now ±2 dB in all conditions; spurious 50–110 Hz junk +18..+32 dB → ≤ +8 dB. Remaining LF nit is *gate churn*, not resolution. |

**And note the omission:** the report never mentions pitch-synchronous / vibrato-demodulated
analysis, which is the fix that actually worked on the problem its §8(a) addresses (rig
0.469 → 0.999; real Female 0.706 → 0.995; Riley: *"sounds like a single singer" for the
first time*). A literature survey that recommends the approach that failed here and omits
the one that succeeded is a survey to read critically.

**Most importantly:** the report does **not** answer Q14/Q15. It never proposes the
architecture, the latent contents, the temporal-scale question, the losses, or the minimal
experiment. Section 4 below fills that gap.

---

## 2. State of the engine

### 2.1 What is built

Pure portable C++, offline file→file, ~2,640 lines in `AdditiveSynthClean/main.cpp`.

```
WAV in
  └─ pitch-sync warp (auto-gated: clarity ≥0.80 ∧ span ≥30¢ → singers only)   main.cpp:611–815
      └─ multi-resolution STFT analysis
           · main tier    4096 Hann  (11.72 Hz/bin, ~47 Hz main lobe)          main.cpp:57
           · LF tier     16384 Hann  (<200 Hz, stationarity-gated flux<0.10)   main.cpp:88
           · transient tier 2048/256 + Edler transition windows                main.cpp:34–51
           · reassignment windows th/dh — built, currently only feeding the
             (retired) phaseMode path                                          main.cpp:68–86
      └─ peak detection → quality filter (floor 60 dB, HF ramp +12 dB,
           amplitude-aware sidelobe exclusion 12 dB/±4 bins)                   main.cpp:305
      └─ parabolic interpolation + fractional-bin phase correction + PV inst-freq
                                                                               main.cpp:367
      └─ tracking: greedy Hz match / mutual-nearest (LF), birth-confirm 2,
           rebirth credit, coast, median-of-3 freq de-jitter <2 kHz,
           asymmetric amp EMA + slope/streak-gated fast release                main.cpp:400–456
      └─ synthesis: per-frame OLA constant-freq (default)                      main.cpp:~2100
           or continuous-phase oscillator bank (synth_mode=1)                  main.cpp:151
      └─ pitch shift: propagated phase, dedupe, freq EMA, rebirth phase
           inheritance, anti-alias roll-off, steady-tone freq lock
      └─ residual: magnitude-deficit random-phase fill above 2500 Hz;
           under shift, deficit vs a parallel UNSHIFTED tonal render, added
           unshifted                                                           main.cpp:2250+
      └─ COLA hole bridge → windowed look-ahead limiter → PCM16 out
  └─ pitch-sync un-warp
```

Hard product constraint (from `docs/voice-artifact-and-plan.md`): **every output sample must
come from summed sinusoids + a modelled residual.** No passthrough of the original. This
constraint is what makes several otherwise-attractive options unusable, and it is the single
most important thing any new architecture must respect.

### 2.2 What helped (committed, ear-validated)

| Commit | Change | Measured effect |
|---|---|---|
| `4d68638` | Asymmetric amp EMA + HF detection sensitivity | HF recovered up to **+12 dB** (SaintSaëns −12.8 → −1.0 dB) |
| `4d68638` | Magnitude-domain stochastic residual | Cymbal attack **37 ms → 2.7 ms** |
| `5480a0c` | COLA overlap-add hole bridge | Real impulsive click on every drum onset, gone |
| `a104867` | Shift-mode dedupe + freq EMA + rebirth phase | Fairlight warble line **41.65% → 0.53%** |
| `510108e` | 16384-pt LF analysis tier | SaintSaëns rumble +18..+32 dB → ≤ +8 dB |
| `de6b65b` | Windowed look-ahead limiter | Levels recovered **+0.6..+4.3 dB** across limited files |
| `ea07722` | Shift-mode amp release | Female shifted ghost tail +2.1 → **+0.1 dB** |
| `e5c2b41` | Median-of-3 trajectory de-jitter | Female f0 jerk 29.6 → 17.7 |
| **`5ae358e`** | **Fractional-bin phase pickup fix** | **saw waveform corr 0.802 → 0.996 — biggest single win** |
| `a1c3745` | PV inst-freq trust + steady-tone lock | 440saw down5 band-envelope p2p 36% → 19% |
| `bacfed2` / `365cc98` | **Pitch-synchronous analysis, auto-gated** | **rig 0.469 → 0.999; real Female 0.706 → 0.995** |

### 2.3 What hurt, and the two meta-lessons

| Attempt | Outcome |
|---|---|
| MQ oscillator bank for unity | Lateral. No audible change; *slightly worse* on doubling (0.44 vs 0.47) — it faithfully renders the same wrong measured phases. |
| Formant-preserving shift (partial envelope) | Worse by ear despite excellent centroid metrics. Fully reverted. |
| R6 frequency-shifted noise residual | Worse by ear. Shifted noise sounds wrong. |
| **`amp_mode=1` energy-integrated amplitude** | **Rig 0.70 → 0.93. Real voice clearly worse** — the ±bin energy sum scoops up breath/noise skirts. |
| `phaseMode` spectral reassignment | Rig 0.469 → 0.206. Worse. |
| Fix-3 ghost classification | Metric didn't move; inaudible; dropped. |
| Harmonic phase lock | Killed the metric (15%/24% → 0%/0%) but cancelled by ear. |
| R3 steady-tone AM, R7 bass octave jitter | **Both were false positives from ratio-based metrics.** Not real defects. |
| Relative match tolerance | Destabilised dense-material tracking (−1.5 dB). |
| `ANALYSIS_SIZE` 2048 | Strictly worse. |

> **Meta-lesson 1 — metrics mislead until validated against the ear.** R3, R7 and R6 were
> all chased on the strength of numbers that turned out to be measuring nothing. The
> battery exists because of this.
>
> **Meta-lesson 2 — a fix that improves a synthetic rig can be worse on real material.**
> `amp_mode=1` and the harmonic lock both did this. **A trained model is a rig-overfitting
> machine**; this history is the strongest local argument for extreme caution about
> learning, and for the report's circularity warning.

### 2.4 Open problems, with their measured numbers

| # | Problem | Status / number | Nature |
|---|---|---|---|
| **A** | Gabor limit on *moving* partials | Solved for monophonic voice via pitch-sync. **Unsolved for choir / polyphony / vibrato instruments.** | Estimation |
| **B** | Relative-phase drift under shift | 300saw down5 band-envelope p2p ~**15%** (target <5%). Independent per-track phase propagation. | Estimation / synthesis |
| **C** | Re-attack under-capture | 2nd crash onset **−3.7 dB**; Fairlight C2 second attack **−4 dB** notch at ~120 ms | Model structure |
| **D** | Cymbal chirp under shift | Ridge glide **3700 Hz vs input 2208** — noise is being tracked as tonal, so it transposes | **Model structure** |
| **E** | Dense mixes | Parked at "additive ceiling", no number | Unknown — Stage 0.2 answers this |
| **F** | Formant/envelope preservation under shift | Absent; one attempt rejected | Model structure |
| **G** | Residual is not a decomposition | Magnitude-deficit fill above 2500 Hz only | **Model structure** |
| **H** | Threshold brittleness | e.g. `lf_max_flux=0.10` vs Fairlight C3 flux p50 = 0.091 → **20 gate flips / 151 frames** | Control logic |

Note the split: **A and B are estimation problems. C, D, F, G are model-structure problems.**
That distinction drives everything below.

---

## 3. Verdict on the central research question

> *Can arbitrary audio be mapped into a structured, multi-resolution, synthesis-oriented
> latent consisting of interpretable harmonic/transient/noise components plus a learned
> residual, reconstructing faithfully while supporting pitch shift, timbral manipulation,
> and component-level editing?*

**In general: yes — the DDSP/S+T+N literature establishes the pieces.**
**For this engine, now: the learned-residual part is the wrong move, and I can show it
from this repo's own outputs.** Four pieces of local evidence.

### 3.1 A magnitude-domain loss is blind to this engine's biggest win — measured

Comparing renders of the same file before and after the phase-fidelity fix (`5ae358e`),
against the input:

| File | render | magnitude LSD vs input | waveform corr vs input |
|---|---|---|---|
| 300hzSaw | `jul3Work_r3` (pre-fix) | 7.315 dB | **0.504** |
| 300hzSaw | `jul3Work_r5` (post-fix) | **7.103 dB** | **0.992** |
| 440sawtooth | `jul3Work_r3` (pre-fix) | 5.183 dB | **0.798** |
| 440sawtooth | `jul3Work_r5` (post-fix) | **4.824 dB** | **0.990** |

The fix that took the sawtooth from visibly-not-a-sawtooth to correct — the largest single
quality gain in this engine's history — is worth **0.2–0.4 dB** on a magnitude objective
and **~0.5 in absolute waveform correlation.**

DDSP-family models are trained on **multi-scale STFT magnitude loss**, which is phase-blind
by construction. Such a training objective would have been *indifferent* to `5ae358e`, would
never have found it, and would happily trade it away for a 0.4 dB LSD improvement elsewhere.
**Any learned system for this engine that uses a magnitude-only reconstruction loss is
disqualified before it starts.** This is not a literature argument; it is this repo's data.

### 3.2 The mel-bottleneck ceiling is *below* current engine fidelity — measured

The `jul3Work_ml/` BigVGAN probe on disk (Jul 5) lets us measure the ceiling of a neural
decoder that takes a mel spectrogram. Running the *original input* through mel → BigVGAN
(full-band 44.1 kHz; band energies preserved, so this is not bandwidth loss):

| File | additive engine vs input | **input → mel → BigVGAN vs input** |
|---|---|---|
| 48kCymbal | 6.54 dB LSD | **8.35 dB** |
| DrumLoopShort | 6.54 dB | **7.35 dB** |
| Female_Sung_Line_3 | 6.75 dB | **7.82 dB** |
| SaintSaensMonoChord | 6.87 dB | **7.57 dB** |

A strong modern neural vocoder, given the *clean original*, reconstructs it **worse than the
current additive engine reconstructs it**. Any architecture that routes the signal through a
mel-resolution bottleneck starts below where we already are. (Caveat: LSD is crude and
BigVGAN synthesises its own phase; but the comparison is consistent across all four files
and the direction is not marginal.)

### 3.3 Pitch shift is a *parametric* operation by requirement, and a learned residual has no defined behaviour under it

The engine's shift rule — sines transpose, noise/reverb does **not**, transients translate in
time — was not chosen for elegance. It was arrived at by three failed attempts:

* R6 frequency-shifted noise residual → worse by ear;
* Female "ghost echo" → diagnosed as the reverb tail being transposed, fixed by keeping it unshifted;
* `residual_hp_hz_shift` 300 Hz → injected unshifted voice-shaped noise at her partials, audibly a second voice.

A learned residual is a vector with no transposition semantics. To make it sound right under
shift you would have to hand-impose exactly the rule above — at which point the learning has
bought nothing in the place where it would matter most.

### 3.4 The decomposability constraint kills the interesting part

"Every output sample comes from summed sinusoids + a modelled residual" is the product
premise. The hp300 true-phase residual sounded great and is unusable precisely because it is
near-passthrough. A *learned* residual that reconstructs faithfully is, by information
content, on the same spectrum — the better it reconstructs, the more it is carrying the
original signal rather than describing it. **Reconstruction fidelity and decomposability are
in direct tension, and the learned residual sits on the wrong side of it.**

### 3.5 What this leaves

The honest positive answer to Q15 is: **the structured, multi-resolution, interpretable
harmonic/transient/noise representation is exactly right — and you build it with DSP.**
That is Stage 2. The "plus a learned residual" clause is the part that fails, for the four
reasons above. And the one place learning is genuinely defensible here is not the signal
path at all — it is the **control path** (problem H), where the engine currently makes
discrete strategy decisions from hand-tuned scalar thresholds that demonstrably fight each
other. That is Stage 3(a).

---

## 4. Direct answers to the 14 questions

| # | Question | Answer |
|---|---|---|
| 1 | **What stays from the C++ engine** | The whole analysis core: multi-resolution STFT (16384 LF / 4096 main / 2048-256 transient), parabolic + PV instantaneous frequency (tuning is already 0.4 ¢ median), the fractional-bin phase correction, tracking with birth/coast/death, OLA synthesis, the shift math + anti-alias roll-off, COLA bridge, windowed limiter, the pitch-sync wrapper, and `battery/`. This is a mature, ear-validated stack — roughly 30 committed fixes, most of which no learned system would rediscover. |
| 2 | **What gets replaced** | (a) The binary threshold stack in `filter_peaks_by_quality` → a **graded** peak-quality score. (b) Greedy `find_best_match_peak_HZ` → one global assignment with a non-stationary cost. (c) The magnitude-deficit residual → a **true S+T+N decomposition**. (d) Per-file-tuned global scalars → material-conditioned parameter sets. |
| 3 | **What gets learned** | Ideally nothing, initially. At most, later: a **strategy/material classifier** in the control path, and possibly a peak-validity refinement — both only after a labelled gap is demonstrated. Never a residual or vocoder in the signal path. |
| 4 | **What stays interpretable** | Everything. There is no reason for anything here to be opaque. |
| 5 | **What the latent contains** | The parameter stream itself: per-partial tracks `(id, f, a, φ, confidence)`; a transient list `(time, per-band gain, decay)`; a noise component as time-varying band envelopes; plus a scalar **representability** measure per frame/band. No opaque vector. |
| 6 | **Frame-by-frame / continuous / multi-scale** | **Multi-scale and event-driven — which it already is, and which is better than the alternative.** A fixed 50 Hz frame latent (the DDSP/HuBERT default) would be a *downgrade*: it throws away the 16384-pt LF resolution and the 256-pt transient resolution that fixed the rumble and the cymbal attacks. |
| 7 | **Training data** | Only relevant at Stage 3. Commercially clean: GuitarSet, NSynth, Slakh2100, MusicNet (original), Vocadito, filtered FSD50K, plus this project's own 15-file corpus and **deliberately synthesised hard negatives** (known sidelobes, cancelling partial pairs, single-frame phantoms). **Never engine output** — `amp_mode=1` is the local proof of what that produces. |
| 8 | **Losses** | Never magnitude-only (§3.1). Use waveform-domain / phase-aware reconstruction + per-partial parameter loss against classical estimates + the battery metrics as hard constraints. Add a relative-phase (waveform-shape) term explicitly. |
| 9 | **How to train** | Small supervised heads only. Parameterise as a **delta from the classical estimate**, initialised so that zero output ≡ today's engine, so the worst case is "no change" rather than "new failure modes". |
| 10 | **How pitch shifting works** | Unchanged and parametric: sinusoids scale by `2^(n/12)` with anti-alias roll-off; noise stays put; transients translate in time; formants optionally re-applied. Empirically validated by three failed attempts to do otherwise (§3.3). |
| 11 | **How a user manipulates it** | Track-level (select/mute/retune/re-envelope partials) and component-level (harmonic / transient / noise balance). Plus a Melodyne-style **material indicator** per file/region driven by the residual meter. Per-frame confidence stays an expert diagnostic, off by default — the report's §7 is right about this. |
| 12 | **How to evaluate "is the learned representation better"** | The existing per-class battery + blind A/B vs the current engine, **with at least one phase-sensitive metric per class**. Any evaluation using magnitude-only metrics is invalid here — §3.1 is the proof. Add the oracle ceiling (below) as the denominator. |
| 13 | **How to detect the residual compensating for an inadequate synthesis model** | **The oracle-fit probe (Stage 0.2).** Optimise the synthesis parameters directly against the input with no analysis stage at all. Whatever residual remains at the oracle optimum is *structural model inadequacy*, by construction. |
| 14 | **Smallest viable prototype** | The oracle-fit probe — no training, no dataset, no C++ change, ~1–2 days. It is also decision-critical: it tells you whether the remaining gap is estimation (learnable in principle) or structure (not). |

---

## 4b. STAGE 0 RESULTS (run 2026-09-12) — what the probes actually found

Stage 0 has been run. Scripts in `/tmp/stage0/` (`lib.py`, `s01_residual_meter.py`,
`s02_ceiling.py`, `s03_trackceiling.py`, `s04_freq_vs_amp.py`); durable pieces landed in
`battery/`. Two env-gated instrumentation hooks added to `main.cpp` (`PCOUNT_DEBUG`,
`PFREQ_DEBUG`), verified byte-identical on five files when unset.

### 4b.1 Residual meter — the representability table

`e = x − x̂` at unity on the current engine. SRR = signal-to-residual ratio (higher =
more of the signal captured; +6 dB = half the residual energy). `flat` = residual
spectral flatness ÷ input flatness.

| class | file | SRR | flat | reading |
|---|---|---|---|---|
| steady | 300hzSine | **23.9** | 234 | residual is pure noise — healthy |
| saw | 300hzSaw / 440saw | **18.0 / 17.1** | 4.2 / 1.6 | healthy |
| wavetable | Fairlight C2 / C3 | **15.2 / 15.3** | 6.0 / 4.9 | healthy |
| pluck | Piano | **12.9** | 3.5 | healthy |
| chord | SaintSaëns | **10.1** | 2.1 | |
| voice | Female | **9.8** | 1.5 | |
| perc | DrumLoop | **4.7** | 1.2 | ← cliff |
| cymbal | 48kCymbal / Out48k | **3.1 / 4.0** | 1.04 / 1.03 | |
| choir | choir-burst | **3.5** | 1.31 | |
| mix | 1985 / Happy / take-me-out | **4.1 / 2.2 / 1.8** | 1.06 / 0.92 / 0.98 | |

Two things fall out immediately:

**(a) There is a cliff, and `flat ≈ 1` marks it.** On every failing class the residual has
the *same spectral character as the input* — the model is not capturing the material at
all. And on choir the residual is **tonal** (absolute flatness 0.044, close to the input's
0.034): what is missing there is *partials*, not noise. That is an estimation defect, not a
missing noise model.

**(b) The unity output is partly passthrough, and it is load-bearing.** At transient regions
`main.cpp:2234` crossfades the **original waveform** back into the output, ungated by any
knob. Measured share of active samples that are literally the original:

| file | passthrough | SRR as shipped | SRR excluding pasted samples |
|---|---|---|---|
| 1985 | **32.4%** | 4.1 | **2.4** |
| DrumLoop | **27.8%** | 4.7 | **1.4** |
| HappyMono | **22.8%** | 2.2 | **1.4** |
| Fairlight C3 | 6.9% | 15.3 | 14.9 |
| most tonal files | <2% | — | ~unchanged |

On percussive material the "additive reconstruction" is substantially the original signal.
This is the sharpest answer to Q13 in the whole exercise: it is not merely that the residual
compensates for an inadequate synthesis model — the *input itself* is doing the
compensating, and it directly contradicts the decomposability premise.

### 4b.2 Oracle ceiling — no class is structurally exhausted

Per-frame matching-pursuit fit of K sinusoids directly to the input, frequencies free,
overlap-added like the engine. A strict upper bound on the synthesis model class.
Probe validated first: pure sine **98 dB**, three static harmonics **95 dB**, white noise
**1.1 dB** at K=32 (capacity-limited, as it must be).

| file | engine | K8 | K16 | K32 | K64 | K128 | K256 |
|---|---|---|---|---|---|---|---|
| 300hzSine | 50.3 | 75.0 | 75.5 | 76.6 | 78.5 | 81.7 | 86.7 |
| 300hzSaw | 19.7 | 11.5 | 14.5 | 17.9 | 24.2 | **73.7** | 77.1 |
| Fairlight C3 | 15.3 | 26.5 | 30.8 | 34.5 | 38.5 | 43.3 | 48.3 |
| Piano | 12.8 | 13.3 | 19.0 | 26.2 | 31.7 | 37.3 | 43.5 |
| SaintSaëns | 10.5 | 13.0 | 19.4 | 25.4 | 30.9 | 35.8 | 40.8 |
| Female | 9.8 | 11.0 | 13.0 | 15.4 | 18.3 | 21.1 | 24.8 |
| **choir** | **3.5** | 11.2 | 14.7 | 18.4 | 22.8 | 28.6 | 34.8 |
| DrumLoop | 6.4 | 8.5 | 10.7 | 12.6 | 14.9 | 18.2 | 23.6 |
| **48kCymbal** | **5.2** | 1.8 | 2.7 | 4.1 | 6.2 | **9.1** | 13.1 |
| HappyMono | 2.3 | 7.3 | 9.5 | 11.1 | 12.9 | 15.1 | 18.5 |
| take-me-out | 1.7 | 5.8 | 7.2 | 8.7 | 10.7 | 13.1 | 16.5 |

**Nothing saturates.** The sinusoidal model class is not exhausted on any class in the
corpus — including dense mixes, which have been parked for months as "the additive
ceiling." They are not at the ceiling; they are 11–16 dB below it.

The one genuine structural case is **cymbals**: the engine (5.2) is close to the oracle
(9.1 at K=128), and the ceiling itself is low (13.1 even at K=256). Spend model-structure
effort there, not estimation effort.

**Budget control.** Measured via `PCOUNT_DEBUG`, mean partials rendered per long frame, and
how many are within 40 dB of the loudest (the audible budget):

| file | rendered | within 40 dB |
|---|---|---|
| Fairlight C3 | 213 | **14** |
| Piano | 166 | 36 |
| SaintSaëns | 842 | 43 |
| choir | 651 | 103 |
| DrumLoop | 427 | 205 |
| Female | 933 | 239 |
| HappyMono | 1270 | 545 |
| 48kCymbal | 894 | 634 |
| take-me-out | 1184 | 754 |

The engine renders 166–1270 partials per frame and is beaten by an oracle using **8 to 32**
on nearly every class. At matched budget it sits **13–25 dB below the ceiling everywhere**.
So the gap is not the partial budget, not the model class, and not the representation.

### 4b.3 The mechanism: inter-partial interference bias

Decisive test. Hold the engine's **own** frequencies for each frame (dumped via
`PFREQ_DEBUG`, engine's own frame geometry) and solve only amplitude and phase, two ways —
independently per partial (what the engine does) versus jointly by least squares.
Engine baseline here is the residual-off render, so all three are tonal-model-only.

| file | ENGINE | INDEP-LS | JOINT-LS | indep gain | **joint extra** |
|---|---|---|---|---|---|
| Fairlight C3 | 15.2 | 8.7 | 31.0 | −6.5 | **+22.3** |
| Piano | 12.8 | 11.0 | 28.2 | −1.8 | **+17.2** |
| SaintSaëns | 10.5 | 11.8 | 30.3 | +1.3 | **+18.5** |
| choir | 3.5 | 12.9 | 29.1 | +9.5 | **+16.1** |
| DrumLoop | 4.5 | 2.5 | 13.9 | −2.0 | **+11.4** |
| take-me-out | 1.7 | 8.3 | 23.0 | +6.6 | **+14.7** |

Null control, so this is not a degrees-of-freedom artefact: fitting **white noise at random
frequencies** with the largest basis used (K=512, 50% DOF in a 2048 frame) reaches only
**4.76 dB** (K=256 → 2.25, K=128 → 1.10). The joint gains of 11–22 dB are real signal.

**Diagnosis.** At a 4096-point Hann the main lobe is ~47 Hz wide. The engine estimates each
partial's amplitude and phase *independently* from its own spectral peak, so in any dense
spectrum neighbouring partials contaminate each other's peak height and phase. That
information is entangled; **no better per-peak reader can recover it.** A joint least-squares
solve over the already-detected frequency set removes it by construction.

This retro-explains most of the project's negative results in one stroke:

* **`amp_mode=1` (energy-integrated amplitude) failed** because it is still an independent
  per-peak estimator, just with a wider aperture — it scooped up *more* neighbour and noise
  energy. Not a bad idea badly tuned; the wrong class of fix.
* **`phaseMode` (spectral reassignment) failed** for the same reason — also per-peak.
* **`5ae358e` was the biggest win** because it removed the one *systematic* component of the
  independent phase read (the fractional-bin residue). Joint estimation removes the rest.
* **The vibrato rig's per-partial amplitude spread (0.375 … 1.12 on equal-amplitude
  harmonics)** is textbook interference between overlapping main lobes.
* **Pitch-sync works** partly because warping concentrates each partial's energy, which
  *reduces* main-lobe overlap.
* **Choir is the worst file in the corpus** because it is both dense (651 partials) and
  moving — maximum interference.

### 4b.4 A relevant negative result: the end-to-end optimisation does not converge

Ceiling B (fully phase-continuous MQ model, frequency and amplitude breakpoints, phase by
integration) was attempted by Adam on a waveform loss, initialised from the matching-pursuit
fit. It **collapsed to ~0 dB on all five files** — far below even the engine. That is an
optimisation failure, not a ceiling, and it is not reported as one.

It is worth recording because it is precisely the obstacle any "learn the synthesis
parameters end-to-end" approach has to clear: gradient descent through an integrated-phase
oscillator bank against a time-domain loss is badly non-convex in frequency (Hayes, Saitis &
Fazekas, ICASSP 2023), and here it failed on real material even from a good initialisation.
It does not prove such training is impossible — DDSP systems mitigate with multi-scale
spectral losses and careful parameterisation — but it does mean that path is not free, and
§3.1 already showed the usual mitigation (magnitude loss) is blind to this engine's
most important defect.

### 4b.5 Harness defects found and fixed

Three real defects, all in the regression harness itself:

1. **`r7.json` is not reproducible from any commit, and the guard rail was vacuous.** The
   committed `battery/corpus.json` declared only `jitter`; `r7.json` contains only
   `saw_corr` / `cepstral_excess` / `env_p2p`. The key sets are **disjoint**, so
   `--baseline r7.json` compared an empty intersection and could only ever print
   "0 regressions." Worse, the r7-era commit (`acf0b18`) rejects the pitch-shift CLI
   argument that `run.py` passes — so `r7.json` was locked from an **uncommitted** tree.
   *Implication:* recent "corpus regression vs r7 = 0 regressions" claims (including
   `365cc98`'s) were not testing anything. I found **no evidence of a hidden regression** —
   pitch-sync was independently re-verified output-neutral on five non-singers — but the
   check was not capable of finding one.
2. **Duplicate entry names.** Both Fairlight entries were named `FairlightAdditiveWaveTA1`,
   so C2 and C3 collided on the `(entry, semi, metric)` scorecard key and each file's
   measurement was compared against the *other* file's baseline. Present in `r7.json` too.
   Fixed, plus a loud `check_unique_names` guard so it cannot recur silently.
3. **The regression test broke on negative values.** `old*(1−tol)` moves the threshold the
   wrong way when `old < 0`, so a −1.07 dB band SRR was flagged as a regression against its
   own identical baseline. Replaced with `abs(old)*tol + abs_slack`.

Landed in `battery/`: `residual_srr` (+ per-band) and `residual_flatness` in `metrics.py`,
declared on **every** entry; manifest restored to `[-5, 0, +5]` with the full metric set and
thresholds; new reproducible baseline `battery/baseline_head.json` (**242 metrics**,
self-consistent, 0 regressions on re-run, 6 threshold failures = the known pre-existing two
`saw_corr` + four Female rows). The engine was separately confirmed deterministic.

### 4b.6 Verdict

* The representation is **not** the problem. The model class has 13–25 dB of headroom on
  every class in the corpus.
* The problem is **parameter estimation**, and specifically **independent per-peak amplitude
  and phase estimation** in the presence of overlapping main lobes. Worth **+11 to +22 dB**,
  recoverable by a convex linear solve.
* The only genuinely structural case is **cymbals** (ceiling ~13 dB), plus the
  **unity passthrough** that is currently standing in for a transient model.
* This *strengthens* §3's conclusion. The largest available win is linear algebra, not
  learning; and the end-to-end optimisation a learned front end would require failed to
  converge here.

---

## 4c. STAGE 1 RESULTS (run 2026-09-13) — joint amp/phase, built and measured

Implemented in `main.cpp` behind CLI arg 13 `jointMode` (0 = off, byte-identical default;
verified on four files). `joint_amp_phase()` + `hann_dtft()` sit above `main()`; the analytic
window DTFT was checked against brute-force sums to ~1e-12 before being wired in.

### 4c.1 The decisive implementation detail: solve at the frequencies you will render

The first version (`jointMode 1`) solved inside the analysis block, at the frequencies
measured there. It gained almost nothing — **+0.95 dB mean**, against a predicted +11 to +22.

Knob sweeps ruled out the solver: band 8 vs 40 bins, regularisation 1e-3 vs 1e-6, and 60 vs
400 CG iterations all produced **identical results to two decimal places**. Bypassing the
amplitude EMA moved it by ~0.1–0.6 dB. The solve was converging; the loss was downstream.

Two causes, found by bisection:

1. **The analysis window is twice the synthesis frame.** Fitting the engine's own frequencies
   over 4096 samples instead of 2048 costs **3–23 dB** (choir −14.0, DrumLoop −22.6) — an
   85 ms stationary-sinusoid fit is simply a worse model than a 43 ms one.
2. **Tracking moves the frequencies after the solve.** The median-of-3 de-jitter, the
   sub-200 Hz smoothing and LF replacement all run *after* analysis. An amplitude and phase
   fitted at `f` is wrong when rendered at `f'`; measured p90 frame-to-frame movement is
   5–7 Hz, which over a 2048-sample frame is **85–111° of phase error**. That alone caps
   achievable SRR in the low single digits, which is exactly where mode 1 sat.

`jointMode 2` re-solves once more at the very end, using each frame's **final** track
frequencies — the ones about to be synthesised. That is the whole difference:

| file | mode 0 | mode 1 | **mode 2** |
|---|---|---|---|
| choir | 3.48 | 4.05 | **17.01** |
| SaintSaëns | 10.11 | 10.94 | **23.65** |
| Piano | 12.87 | 15.17 | **29.85** |
| Fairlight C3 | 15.29 | 16.29 | **29.45** |
| take-me-out | 1.77 | 2.12 | **19.87** |
| Female | 9.76 | 10.29 | **21.88** |

**Lesson worth keeping: estimate on the basis you synthesise.** Any future estimator change
has to land after the tracker, or the tracker has to stop moving frequencies.

### 4c.2 Result vs the Stage-0 ceiling

`residual_srr_db` at unity, full corpus, versus the oracle ceiling from §4b.2:

| file | base | joint | gain | ceil K128 | headroom closed |
|---|---|---|---|---|---|
| take-me-out | 1.77 | **16.52** | +14.75 | 13.1 | 130% |
| HappyMono | 2.19 | **15.73** | +13.54 | 15.1 | 105% |
| Female | 9.76 | **19.46** | +9.71 | 21.1 | 85% |
| PianoSampleMono | 12.87 | **25.38** | +12.51 | 37.3 | 51% |
| choir | 3.48 | **15.47** | +11.99 | 28.6 | 48% |
| SaintSaëns | 10.11 | **21.86** | +11.75 | 35.8 | 46% |
| Fairlight C3 | 15.29 | **23.95** | +8.66 | 43.3 | 31% |
| DrumLoopShort | 4.73 | **8.92** | +4.19 | 18.2 | 31% |
| 48kCymbal | 3.12 | **4.61** | +1.49 | 9.1 | 25% |
| 300hzSaw | 18.03 | **22.59** | +4.56 | 73.7 | 8% |
| 1985 | 4.08 | **16.40** | +12.32 | — | — |
| Fairlight C2 | 15.22 | **26.78** | +11.56 | — | — |

**Better on 15/15 files. Median 46% of the ceiling headroom closed.** The two mixes exceed
the K=128 oracle because the engine legitimately renders far more partials than 128 (754 and
545 significant). The Stage-0 diagnosis was correct and the predicted magnitude was real.

`saw_corr` also improved on both saws (0.9918 → 0.9970, 0.9903 → 0.9938), and
`cepstral_excess` — the doubling metric — improved on 10/12 including Female (4.35 → 3.04).

### 4c.3 A second harness bug: absolute metrics punished the fix

The runner initially reported **49 regressions** while fidelity improved 4–18 dB everywhere.
Cause: `envelope_p2p` and `trajectory_jitter` measure the **output alone**, so "less
variation" always scores better — even when the input genuinely has that variation.

Measured against the input:

| file | env_p2p: input / base / joint | jitter: input / base / joint |
|---|---|---|
| Female | 8.22 / 7.75 / **7.94** | 2.14 / 1.58 / **1.96** |
| choir | 3.73 / 4.12 / **3.72** | 1.32 / 1.07 / **1.18** |
| Out48kCymbal | 7.70 / 7.10 / **7.52** | 1.39 / 1.06 / **1.11** |
| Fairlight C3 | 2.79 / 2.87 / **2.82** | 1.23 / 0.89 / **0.93** |

On 6 of 7 files both metrics moved **closer to the input**. The per-peak baseline was
*over-smoothing*; the battery was rewarding it for that. This is the same failure mode as R3
and R7 — an absolute metric manufacturing a defect that direct measurement disproves.

Fixed by adding `envelope_p2p_dev` and `trajectory_jitter_dev` (deviation from the reference)
in `metrics.py`, declared on every entry. On the corrected metrics joint mode is better on
**13/15** (`env_p2p_dev_full`), **13/15** (`env_p2p_dev_band`) and **10/15** (`jitter_dev_db`).
Keep the absolute form only for the steady-tone gate, where "no wobble" is the truth.

### 4c.4 Open risks — the honest caveats

* **Not yet heard.** Every number here is objective. This project's record is that metrics
  have overstated R3, R6 and R7, and that `amp_mode=1` improved a rig while sounding worse.
  A 45-file listening set is rendered to `~/Desktop/comparingOut/jul3Work_joint/`
  (`_orig`/`_up5`/`_down5`). **This is the gate.**
* **Degrees of freedom / noise fitting.** On the densest material the basis is large relative
  to the window: take-me-out is 1184 partials = 2368 unknowns against 4096 samples (58% DOF).
  The Stage-0 null control says a basis that size fits *white noise at random frequencies* to
  ~4.8 dB, so on the mixes perhaps 5 dB of the gain is the solve absorbing noise into
  partials. That is the classic musical-noise/birdie mechanism, and it is the first thing to
  listen for on 1985 / HappyMono / take-me-out. Mitigation if it is audible: cap the joint
  basis by amplitude (partials within N dB of the loudest) rather than solving all of them.
* **Minority objective regressions**, all small but real: Fairlight C2 `cepstral_excess`
  13.18 → 20.23 (the largest single regression anywhere), choir 5.60 → 6.45,
  Out48kCymbal `env_p2p_dev_band` 0.40 → 0.61, and `jitter_dev` worse on cymbal, drums,
  Fairlight C2, Happy and piano.
* **Short frames are untouched.** The pass only covers `LONG_SIZE` frames on the 4096 phase
  reference, so transient regions get none of the benefit — visible as DrumLoop's modest
  +4.19 and cymbals' +1.5. That is Stage 2 territory, not a defect of this change.
* **Cost:** ~0.6 s per file. Negligible.

### 4c.5 Verdict

Stage 1 delivers what Stage 0 predicted: the engine's dominant error was inter-partial
interference in per-peak amplitude/phase estimation, and a convex joint solve recovers
9–15 dB on most of the corpus with no learned component. It stays off by default pending
the listening test.

**Stage 2 is now clearly the right follow-up**, and Stage 1 sharpens it: the classes that
gained least are exactly the ones Stage 0 flagged as structure-limited (cymbals +1.5,
drums +4.2), and those are the transient/noise cases. Do Stage 2 next.

---

## 4d. RE-BASELINE (run 2026-09-20) — where the headroom is after joint LS + harmonic lock

Same probes (`s01`, `s02`), same corpus, engine at `74d68d4` (joint LS default, pitch-sync
auto, harmonic lock). Purpose: re-read the ceiling table now that the estimation fix is in,
before committing Stage 2 effort.

### 4d.1 Residual meter, then vs now (unity SRR, dB)

| file | Sep 12 | Sep 20 | pass% now | SRR excl. passthrough |
|---|---|---|---|---|
| 300hzSine | 23.9 | **24.1** | 18.4 | 25.9 |
| 300hzSaw / 440saw | 18.0 / 17.1 | **22.6 / 19.1** | 6.3 / 1.7 | 24.2 / 19.2 |
| Fairlight C2 / C3 | 15.2 / 15.3 | **26.8 / 23.9** | 1.8 / 7.1 | 26.7 / 23.6 |
| Piano | 12.9 | **25.4** | 4.4 | 25.1 |
| SaintSaëns | 10.1 | **21.9** | 0.4 | 21.8 |
| Female | 9.8 | **19.5** | 1.9 | 19.5 |
| choir | 3.5 | **15.5** | 0.5 | 15.5 |
| DrumLoop | 4.7 | **8.9** | 28.5 | 5.1 |
| 48kCymbal / Out48k | 3.1 / 4.0 | **4.6 / 5.4** | 2.4 / 2.2 | 5.4 / 6.2 |
| 1985 / Happy / take-me-out | 4.1 / 2.2 / 1.8 | **16.4 / 15.7 / 16.5** | 32.7 / 23.0 / 0.3 | 14.4 / 14.7 / 16.6 |

The cliff is gone on every tonal class and on the mixes. What remains low is percussive
(DrumLoop, cymbals) — and §4d.3 shows that number is not what it looks like.

### 4d.2 Oracle ceiling, then vs now (engine SRR on the probe's interior region)

| file | engine Sep 12 | engine Sep 20 | audible budget | oracle at ~that K | K256 |
|---|---|---|---|---|---|
| 300hzSine | 50.3 | **65.2** | few | 75.0 (K8) | 86.7 |
| 300hzSaw | 19.7 | **36.0** | ~100 | 73.7 (K128) | 77.1 |
| Fairlight C3 | 15.3 | **23.9** | 14 | 30.8 (K16) | 48.3 |
| Piano | 12.8 | **25.3** | 36 | 31.7 (K64) | 43.5 |
| SaintSaëns | 10.5 | **23.2** | 43 | 30.9 (K64) | 40.8 |
| Female | 9.8 | **20.3** | 239 | 24.8 (K256) | 24.8 |
| choir | 3.5 | **15.5** | 103 | 28.6 (K128) | 34.8 |
| DrumLoop | 6.4 | **18.5** | 205 | 23.6 (K256) | 23.6 |
| 48kCymbal | 5.2 | **10.9** | 634 | 13.1 (K256) | 13.1 |
| HappyMono | 2.3 | **16.4** | 545 | 18.5 (K256) | 18.5 |
| take-me-out | 1.7 | **17.3** | 754 | 16.5 (K256) | 16.5 |

Readings:

* **The dense mixes and cymbals have reached the K=256 oracle.** The probe no longer bounds
  them; a K≥1024 oracle would be needed to know their remaining headroom (slope ~+3 dB per
  doubling of K suggests ~5 dB). They are no longer the priority.
* **Remaining matched-budget headroom, ranked:** choir **~13 dB** (residual still tonal:
  missing/mis-estimated partials, many voices with independent vibrato — the per-track
  demodulation case, §5.1-OLD 1.3), Fairlight C3 **~7 dB to K16** (14 audible partials and
  the oracle beats it with 16 — a time-varying wavetable; estimation under spectral motion),
  Piano / SaintSaëns **~7–8 dB**, Female **~4.5 dB**, saw **large but inaudible** (already
  0.999 shape).
* Cymbals: engine 10.9 vs 13.1 at K256 — still the structural case, but the gap is 2 dB, not
  the 8 dB it was. A noise model remains the right tool; it is no longer urgent for SRR.

### 4d.3 The percussive "transient problem" is mostly a FILE-START defect

The probe's interior region excludes the first/last ~2% of samples. On DrumLoop that
region scores 18.5 dB while the whole file scores 8.9; on 48kCymbal 10.9 vs 4.6. Locating
the residual in 50 ms blocks:

| file | whole-file SRR | first 50 ms: share of residual | share of signal | block SRR |
|---|---|---|---|---|
| DrumLoop | 8.9 | **90%** | 15% | 1.0 |
| 48kCymbal | 4.6 | **82%** | 30% | 0.3 |
| HappyMono | 15.7 | **22%** | 0.8% | 1.4 |

And the engine's output level in 10 ms blocks from t=0, as a ratio to the input:
DrumLoop `0.08 0.03 0.09 0.96 1.06`, Happy `0.13 0.09 0.12 0.52 0.76 0.95`, 300hzSine
`0.00 0.15 0.72 0.99 1.00`. **The engine renders almost nothing for the first ~30–40 ms of
every file.** Frame 0 already has a rect-fade-to-Hann OLA window (`main.cpp:2654`), so it
is not the synthesis ramp; it is consistent with `peak_birth_confirm_frames = 2` — every
track is newborn at frame 0 and nothing is confirmed until the third hop (~43 ms at hop
1024) — with no transient-region passthrough covering it because there is no prior frame
for the onset detector to compare against.

On files that begin on a hit (DrumLoop, both cymbals) this single defect is 80–90% of the
total residual and the reason the percussive classes sat at the bottom of every table
since July. It is audible (the first hit of a loop is the one the listener notices), and it
is a small fix, not a transient model.

**Consequence for the plan:** fix the file-start birth latency first (§5, Step 1a). Then
re-measure DrumLoop and cymbals. Only what remains after that is transient-model work.

---

### 4d.4 Step 1a landed (2026-09-21) — file-start credit, and the ceiling out to K=1024

`bbc8d6a`. Riley's ear verdict: approved ("everything seems to be in effect"). Whole-file
unity SRR after: sine 42.7, 300saw 35.1, 440saw 28.5, DrumLoop 13.8, cymbals 8.8 / 10.2,
Happy 16.5; everything else unchanged. What remains at the file start is now the onset
itself: DrumLoop's first 50 ms still holds 70% of its residual (block SRR 7.0), the cymbal's
54%, because frame 0 cannot qualify as a transient (`transientNegotiationTactics` requires
`f > 0`) and so is rendered with 4096-frame sinusoids. That is the next file-start step.

Oracle extended to K=1024 (interior region, engine with credit on):

| file | engine | K32 | K128 | K256 | K512 | K1024 | audible budget | gap at ~budget |
|---|---|---|---|---|---|---|---|---|
| 300hzSine | 65.2 | 76.6 | 81.7 | 86.7 | 92.4 | 98.6 | few | ~10 (inaudible) |
| 300hzSaw | 36.0 | 17.9 | 73.7 | 77.1 | 82.5 | 89.4 | ~100 | ~38 (shape already 0.999) |
| 440saw (vibrato) | 29.5 | 17.3 | 23.3 | 26.6 | 30.9 | 42.2 | ~100 | **0 — at ceiling** |
| Fairlight C2 | 27.0 | 38.2 | 48.4 | 53.7 | 59.4 | 65.2 | 14 | **~11** |
| Fairlight C3 | 23.9 | 34.5 | 43.3 | 48.3 | 53.4 | 58.6 | 14 | **~11** |
| Piano | 25.3 | 26.2 | 37.3 | 43.5 | 49.6 | 55.6 | 36 | ~7 |
| SaintSaëns | 23.2 | 25.4 | 35.8 | 40.8 | 44.9 | 48.7 | 43 | ~8 |
| Female | 20.3 | 15.4 | 21.1 | 24.8 | 30.0 | 36.8 | 239 | ~4.5 |
| choir | 15.5 | 18.4 | 28.6 | 34.8 | 40.2 | 45.5 | 103 | **~13** |
| DrumLoop | 18.5 | 12.6 | 18.2 | 23.6 | 32.0 | 40.0 | 205 | ~5 (+ file-start onset) |
| 48kCymbal / Out48k | 10.9 / 11.0 | 4.1 / 4.6 | 9.1 / 10.0 | 13.1 / 14.5 | 19.0 / 20.6 | 27.2 / 28.7 | 634 | ~8 nominal, see note |
| 1985 | 16.5 | 11.3 | 17.7 | 21.7 | 27.1 | 34.0 | ~500 | ~8–10 |
| HappyMono | 16.4 | 11.1 | 15.1 | 18.5 | 24.1 | 32.1 | 545 | ~8 |
| take-me-out | 17.3 | 8.7 | 13.1 | 16.5 | 22.4 | 31.2 | 754 | ~5–8 |

**Note on K≥512:** 1024 sinusoids in a 4096 frame is 50% of the DOF; the Stage-0 null control
showed a 512-partial basis fits white noise to 4.8 dB. So the K512/K1024 columns on noisy
material (cymbals, drums, mixes) are partly noise-fitting, not a sinusoidal ceiling, and the
"gap" there overstates what a better sinusoidal estimator could recover. The cymbal's real
sinusoidal ceiling is still the structural one; a noise component is the right tool.

**Reading, ranked by recoverable headroom:**
1. **choir ~13 dB** — residual still tonal. Many voices with independent vibrato: the
   per-track demodulation case (§5.1-OLD 1.3).
2. **Fairlight C2/C3 ~11 dB** — the oracle beats the engine with 16 partials where the
   engine renders 213. A time-varying wavetable: estimation under spectral motion. Cheap
   to investigate — it is a synthetic file, the truth is knowable.
3. **Piano / SaintSaëns ~7–8 dB** — dense stationary partials; likely closely-spaced pairs
   (SaintSaëns) and inharmonic partials with beating (piano). Candidate for per-band
   high-resolution estimation.
4. **Mixes ~5–10 dB** — but discounted by the noise-fit caveat, and they hold 20–33%
   passthrough at unity (Happy, 1985). The honest number needs the transient model.
5. **Female ~4.5 dB** — close.
6. **Percussive file-start onset** — DrumLoop/cymbal first 50 ms is still 54–70% of their
   residual; frame 0 as transient is the next small step.
7. **At ceiling:** sine, 440saw, 300saw (perceptually).

### 4d.5 External listening feedback (Sep 21) — "the drum loop took a turn for the worse"

A second listener on the Sep 21 set: SaintSaëns and Female much better; DrumLoop worse;
suggested comparing unity renders of DrumLoop/Happy "to see what the transients are doing".

Measured (per-onset, 5 ms envelope, RMS-matched, envelope-aligned; 10 DrumLoop onsets):

| DrumLoop | onset peak vs input (mean / worst) | pre-onset energy vs input (mean / worst) | attack delay |
|---|---|---|---|
| unity, any round | +0.3 / −0.0 dB | +0.0 / +0.3 dB | 1.5 ms |
| r7 up5 (July) | +0.9 / −6.7 | +5.2 / +11.2 | 5 ms |
| joint (Sep 13) up5 | −0.8 / −7.5 | **+8.1 / +22.1** | 7 ms |
| Sep 21 up5 | −0.1 / −8.1 | +7.7 / +19.0 | 4 ms |
| Sep 21 down5 | −0.9 / −9.9 | +6.1 / +21.6 | 6 ms |

**The feedback is valid, and it is shift-only.** At unity the drum envelope is essentially
perfect — because 28% of DrumLoop's unity samples are the pasted original (§4b.1b). Under
shift there is no passthrough, so transients are rebuilt from 4096-frame sinusoids: the
hit's energy is smeared up to 85 ms before the onset (+6–8 dB mean, +19–22 dB before the
1.62 s hit), the weakest hits are 8–10 dB short, and attacks land 4–7 ms late. The listener's
instinct to compare unity was exactly right: unity hides it.

**Attribution (A/B on Sep 21 binary, DrumLoop up5):** joint LS off restores pre-onset to
+5.2 mean / +13.1 worst (≈ r7). Harmonic lock, pitch-sync and the file-start credit change
nothing. So joint estimation, which is worth +10 dB on every tonal class, made shifted drum
smear ~2.5 dB worse on average and ~6 dB at the worst hit: the joint solve absorbs the
transient's energy into stationary sinusoids more completely than the per-peak read did,
and a stationary sinusoid spreads that energy across the whole frame.

**Tried and rejected:** generalising the retroactive birth credit to every onset (render
each track from its first observation). Pre-onset got worse (+7.7 → +9.1), attack
unchanged (+4.0 → +4.5 ms). Confirmation latency is not the mechanism; the frame length is.

**Consequence:** the transient component (plan Stage 2 item 2) is promoted to next. Its
minimal form is directly targeted at this: under shift, in transient regions, add the
unity residual `x − x̂_unity` **unshifted** (transients translate, they do not transpose)
in place of the sinusoidal transient, exactly as the stochastic residual fill already does
for noise (4d85b5f) but time-domain and phase-true, gated to onset regions. Gate: shifted
DrumLoop pre-onset mean < +2 dB, worst < +6, peaks within ±3 dB; same on Happy/1985/
take-me-out; tonal files byte-identical (no transient regions).

## 5. The plan

Ordered so that each stage is independently shippable and each gate can stop the next.

### Stage 0 — Truth instrumentation (days, zero risk, pure upside)

**0.1 Residual meter.** Compute `e = x − x̂` per frame per band → residual energy fraction +
residual spectral flatness. The engine already builds a parallel unshifted model
(`synth_unity`, `main.cpp:1774`), so this works at unity *and* under shift. Ship first as an
engineering diagnostic; later as the per-file/per-region "additive suitability" indicator.
Non-circular by construction. This is the report's best idea and I agree with it fully.

**0.2 Oracle-fit ceiling probe — the decision-critical experiment.**
In Python/torch, reimplement only the *synthesis* math (sum of sinusoids + the existing noise
fill) and optimise the parameters **directly against the input waveform**, bypassing analysis
entirely. Run on the vibrato rig, a saw, the Female, a cymbal, and a dense mix.

* If oracle-fit reaches near-perfect on a class → the synthesis model is adequate there and
  the whole gap is **estimation** → Stage 1 is the right work, and a learned estimator is at
  least *conceivable* later.
* If oracle-fit still fails on a class (my prediction: cymbals and dense mixes) → the
  **model is structurally short of components** → no encoder, learned or otherwise, can fix
  it → Stage 2 is mandatory and Stage 3 is moot for those classes.

This answers Q13 and Q14 at once, gives problem **E** its first real number, and costs days.

**0.3 Battery hardening.** Add `waveform_correlation` and residual-fraction as **tracked
metrics for every class**, not just saws. Rationale is §3.1: the battery currently has one
phase-sensitive metric and it covers two files, so most of the corpus is guarded only by
metrics that would have missed the biggest fix ever made here. Wire to CI; re-lock baseline.

**Gate:** 0.2 produces a per-class ceiling table. That table decides Stage 1 vs Stage 2 priority.

---

> **REVISED after Stage 0 (2026-09-12).** Stage 0 is done (§4b). Its results reorder what
> follows: the old Stage 1 was about better *detection and tracking*, but the data says the
> engine's frequencies are already fine and detection/tracking is **not** the binding
> constraint. Stage 1 is replaced by **§5.1-NEW joint amplitude/phase estimation**, which is
> now the single highest-value piece of work in the project. The old Stage 1 items are
> demoted to §5.1-OLD (do later, or not at all). Stage 2 survives but narrows.

### Stage 1-NEW (do this next) — joint amplitude/phase estimation

**The change.** Keep peak detection, tracking, the LF tier, pitch-sync and the frequency
estimates exactly as they are. Replace the *independent per-peak* amplitude and phase read
with a **per-frame joint least-squares solve** over the frequencies already detected for
that frame: build the cos/sin basis at those frequencies, solve for the coefficients against
the windowed frame, and take amplitude and phase from the solution.

**Why this and not the peak-quality ensemble.** Measured: +11 to +22 dB SRR on every class
(§4b.3), versus a frequency-side headroom of 0–17 dB that a detection/tracking overhaul would
chase at much higher risk. It also explains and supersedes two failed experiments
(`amp_mode=1`, `phaseMode`) rather than repeating their mistake.

**Feasibility.** Convex, closed-form, no iteration, no training. Partials only interact
within a few bins of a 4096-point Hann, so the normal-equations matrix is **banded** —
sort by frequency, and cost is near-linear in partial count rather than cubic. This is an
offline renderer, so even a modest constant factor is affordable.

**Known risks, with mitigations already present in the engine:**
* *Ill-conditioning* when two tracks sit at nearly the same frequency. The engine already has
  dedupe machinery (`shift_dedup_bins`) and a near-equal margin; add Tikhonov regularisation
  and reuse the dedupe for the joint solve's basis.
* *It changes unity output everywhere* — this cannot be a byte-identical change, so it is
  ear-gated, not metric-gated alone. Given the project's history (metrics have overstated
  R3/R6/R7, and the harmonic-lock round was cancelled by ear), expect to A/B it per class.
* *The 40 dB-down partials.* Fairlight renders 213 partials of which 14 matter. Including
  199 junk basis vectors in a joint solve invites overfitting noise into "partials." Cap the
  basis by amplitude (the null control quantifies the cost: ≤4.8 dB at 50% DOF).

**Gate.** `residual_srr_db` in the battery, per class, against `battery/baseline_head.json`:
choir 3.5 → >20, take-me-out 1.8 → >12, SaintSaëns 10.1 → >25, Fairlight C3 15.3 → >28,
no class regressed, `saw_corr` ≥ current 0.992/0.990. Then ear, per class.

### Stage 1-OLD (demoted) — detection and tracking

**1.1 Graded peak-quality ensemble.** Replace the binary threshold cascade in
`filter_peaks_by_quality` (`main.cpp:305`) with a fused score: Rodet SLM (main-lobe
cross-correlation), PV-vs-parabolic **disagreement** (free — both are already computed),
cross-frame fractional-offset variance, and reassignment agreement. **The `analysis_th` /
`analysis_dh` windows for this are already built** (`main.cpp:68–86`) from the retired
`phaseMode` experiment — reassignment failed as a *phase correction* but is well-suited as a
*validity feature*. Score feeds birth/coast/kill instead of a hard accept/reject.
*Targets:* dense-mix birdies (problem E), threshold brittleness (H).
*Risk:* threshold changes have destabilised dense tracking before — guard hard on SaintSaëns
(over-rejection canary) and Happy/1985 (under-rejection canary).

**1.2 Global track assignment.** Replace greedy `find_best_match_peak_HZ` (`main.cpp:400`)
with one Hungarian/linear-assignment step over all peaks and tracks, with a non-stationary
cost (frequency **slope**, not just frequency distance), generalising the mutual-nearest
matcher the LF tier already uses (`main.cpp:427`). *Targets:* track stealing, vibrato track
splitting, octave hops. *Note:* keep the relative-tolerance lesson in mind — relative
tolerance alone failed; a proper global cost is a different thing.

**1.3 Per-track demodulated analysis — the generalisation of pitch-sync (problem A).**
Pitch-sync works because a monophonic voice has one f0; choir and polyphony have none, which
is why the auto-gate correctly refuses to engage on them. The generalisation is to demodulate
**each track against its own smoothed frequency trajectory** before estimating its amplitude
and phase, rather than warping the whole signal. This is the real frontier for A and the
natural successor to the biggest recent win. Extend `vibrato_rig.py` to a two-voice /
polyphonic rig first so there is a gate before any code.

**1.4 Harmonic phase coherence (problem B) — low priority.** The trial killed the metric
(15%/24% → 0%/0%) but Riley cancelled the round by ear. Revisit only if saw-shift wobble
resurfaces in listening. The patch is archived at
`jul3Work_analysis/harmonic_lock_trial.patch`.

**Gate:** battery green with no class regressed; blind A/B on affected classes.

---

### Stage 2 — Sines + Transients + Noise (narrowed by Stage 0, still required)

> **Narrowed.** Stage 0 showed most classes are estimation-limited, not structure-limited, so
> S+T+N is no longer "the main event" for fidelity. Two specific jobs remain, and both are
> real:
>
> 1. **Cymbals are the one structural case** — engine 5.2 dB vs a ceiling of 9.1 (K=128) /
>    13.1 (K=256). No estimator fixes that; only a noise component does. This is also the
>    open cymbal-chirp-under-shift problem (glide 3700 vs input 2208), same root.
> 2. **Replace the unity passthrough with a real transient model.** The engine currently
>    pastes the original waveform over 23–32% of active samples on percussive material
>    (§4b.1b). That is the single largest violation of the decomposability premise in the
>    codebase, and it is invisible in every metric because pasting the original scores
>    perfectly. Removing it will *lower* measured SRR on drums before a transient model
>    raises it again — expect and accept that.
>
> Run Stage 2 **after** Stage 1-NEW: joint estimation changes what the residual contains, so
> fitting a noise model first would fit it to the wrong residual.

#### Stage 2 detail

**This is the "structured interpretable representation" the research question asks for, built
classically.** Problems C, D, F and G are all the same root: the engine has a sinusoidal
model and a *noise patch*, not a decomposition.

Concretely, today: noise is tracked as tonal (that is why the cymbal glides 3700 vs 2208 Hz
under shift), transients are rebuilt from confirmed sinusoids (they can't be — hence the
−3.7 dB re-attack), and the residual is a magnitude-deficit fill above 2500 Hz with no
identity of its own.

Target decomposition:
* **Sines** — the existing track model. Transpose under shift.
* **Transients** — explicit detection (Daudet / Verma-Meng style), stored as `(time, per-band
  gain, decay)`. **Translate in time under shift; do not transpose.** Directly addresses C.
* **Noise** — time-varying band envelopes fitted to `x − sines − transients` across the *full*
  band, not just above 2500 Hz. **Stays at pitch under shift.** Directly addresses D and G.
* **Envelope** — optionally a true-envelope estimate for formant preservation (F), behind a
  switch, given the prior ear rejection.

Payoffs beyond the artifacts: it gives the user component-level editing (Q11), it makes the
shift rule structural rather than a set of gated heuristics, and it makes the residual meter
meaningful per-component.

**Gate:** cymbal glide under shift approaches input (3700 → <2600); 2nd-crash onset deficit
−3.7 → <−1 dB; no tonal class regresses; ear.

---

### Stage 3 — The learning decision (only after Stage 1-NEW and Stage 2)

> **Stage 0 verdict on this stage: the case for learning got weaker, not stronger.** The
> largest measured win in the project (+11–22 dB) is a convex linear solve, not a learned
> function. The gap was never pattern recognition. And the optimisation that an end-to-end
> learned front end would have to solve failed to converge here from a good initialisation
> (§4b.4). So Stage 3 stays where §3 put it — the control path — and it now has to wait
> behind two classical stages that are known to be worth more.
>
> **Revisit the question after Stage 1-NEW lands.** That is the honest trigger: once joint
> estimation is in, re-run the ceiling probe. If a class is *still* far from its ceiling and
> the remaining error is not explained by noise/transient structure, that residual gap is the
> first genuinely defensible target for a learned estimator — and by then there will be a
> reproducible per-class baseline to measure it against, which there was not before today.

Ranked by risk. Do not skip ahead.

**(a) Strategy / material classifier — lowest risk, most defensible.**
Not in the signal path. It picks *which validated code path* runs: LF tier on/off (replacing
`lf_max_flux=0.10`, which flips 20×/151 frames on Fairlight C3), transient sensitivity,
pitch-sync engage (replacing the `clarity ≥0.80 ∧ span ≥30¢` scalar gate), residual strategy.
Its worst case is bounded by the worst existing path — it cannot invent a new artifact.
Supervision is cheap and non-circular: the *battery outcome* of each strategy per region is
the label. This is the Melodyne pattern the report's §7 describes, and it fits this codebase
better than anything else learned.

**(b) Peak-validity refinement.** Only if Stage 1.1's ensemble leaves a measurable, audible
gap against a few hundred hand-labelled real peaks. The report's threshold is right: if the
classical ensemble is already >~95% precision/recall, **do not ship a model.**

**(c) Parameter-correction encoder** — the minimal experiment from the brief, corrected:
supervise from classical estimates + hard negatives + real labels (never engine output);
phase-aware loss (never magnitude-only); and parameterised as a **delta** that is the identity
at zero. Only if 0.2 said the gap is estimation-limited *and* Stage 1 plateaued short of the
oracle ceiling.

**Never:** a learned residual in the signal path (§3.3, §3.4), or a neural vocoder at the
output (§3.2).

**Deployment, if (a)/(b)/(c) ever ships:** sub-1M params, hand-rolled inference. No ONNX
Runtime, no LibTorch — the report is right that a general runtime buys nothing at this size
and costs binary-size/ODR/notarisation risk.

---

### Stage 4 — Release gate

Per-class battery CI green, no class regressed; blind A/B vs the current engine per class;
MUSHRA with ~10 experienced listeners only if this goes commercial. Frame the gate as *"does
it beat the current engine on every class"*, not *"is it good"*.

---

## 6. Revert / retire list

| Item | Recommendation |
|---|---|
| `amp_mode=1` (energy-integrated amplitude) | **Retire.** Validated dead end on real material, documented in `docs/voice-artifact-and-plan.md`. Dead code in a hot path. |
| `phaseMode` reassignment **phase correction** | **Retire the application** (0.469 → 0.206). **Keep the `analysis_th`/`analysis_dh` machinery** — Stage 1.1 reuses it as a validity feature. |
| `residual_mode=1` (true-phase residual) | **Keep as diagnostic only, clearly marked unusable** — it violates decomposability (near-passthrough). |
| `synth_mode=1` (MQ oscillator bank) | **Keep, but do not invest.** Never won on unity, but it is the only continuous-phase synthesiser available and problem B lives there. |
| Per-file-tuned scalar thresholds | **Do not revert — supersede** in Stage 3(a). They are ear-validated; the problem is that they are global constants doing a per-material job. |

Nothing in the committed history should be reverted. The engine's trajectory has been
sound; the losses were all caught and reverted at the time.

---

## 6b. What is in the working tree after Stage 0 (uncommitted)

| file | change | output-affecting? |
|---|---|---|
| `AdditiveSynthClean/main.cpp` | `PCOUNT_DEBUG` + `PFREQ_DEBUG` env-gated dumps in the synthesis loop | **No** — verified byte-identical on 5 files when unset |
| `battery/metrics.py` | `residual_srr` (+ per-band), `residual_flatness`, `_align`, docstring on why a phase-sensitive metric was needed | tooling |
| `battery/run.py` | wire up `residual`; `check_unique_names` guard; fix regression test for negative values | tooling |
| `battery/corpus.json` | restored `[-5,0,+5]` + full metric set + thresholds; `residual` on every entry; **de-duplicated the two Fairlight names** | tooling |
| `battery/corpus.example.json` | `residual` added | tooling |
| `battery/baseline_head.json` | **new** reproducible baseline at HEAD (242 metrics) | tooling |
| `battery/probes/` | **new** — the Stage-0 probes + README | tooling |
| `docs/engine-direction-2026-09.md` | this document | docs |

Also still present from before: the uncommitted **Catmull-Rom resampler** in the pitch-sync
warp. Unrelated to Stage 0; still wants a deliberate commit-or-discard.

`r7.json` and `jitter_baseline.json` are now superseded by `battery/baseline_head.json`.
Keep `r7.json` for historical reference but **do not** use it as a gate — §4b.5 explains why
it cannot be reproduced.

## 7. Housekeeping

* Uncommitted in the working tree: the **Catmull-Rom resampler** for the pitch-sync warp
  (`main.cpp:772–786`), fixing an audible rattle on bright sung notes from linear
  interpolation. This looks correct and finished — commit or discard it deliberately rather
  than leaving it floating.
* Branch `oscillatorParameter` is **in sync with its remote** but **28 commits ahead of
  `main`** — the memory note saying "unpushed" is stale. Decide whether main should be caught up.
* `.DS_Store` and Xcode `xcuserstate` are tracked and dirty — add to `.gitignore`.
* `main.cpp` is 2,640 lines with 12 positional CLI args and ~40 tuned constants. Not urgent,
  but Stage 2 will push it past maintainable; plan a split into translation units with the
  settings block extracted into one documented struct.

---

## 8. What would change this assessment

* **Stage 0.2 shows oracle-fit is near-perfect on cymbals and dense mixes** → the model is
  adequate everywhere and the entire gap is estimation. Stage 2 shrinks; a learned estimator
  becomes much more interesting. *(I doubt this, but it is the cheapest way to find out.)*
* **Stage 1.1's classical ensemble leaves a big labelled gap on dense mixes** → the report's
  Stage 3 learned gate is justified for that class only.
* **Stage 2's S+T+N lands and the engine still can't shift cymbals convincingly** → that is
  the point at which a genuinely different representation (not a bolt-on) deserves a
  rethink — and the oracle ceiling from 0.2 would already tell us how much room is left.
* **A phase-aware neural decoder benchmark beats 6.5 dB LSD on this corpus** → revisit §3.2.
  The BigVGAN number is a mel-bottleneck ceiling, not a universal one.
