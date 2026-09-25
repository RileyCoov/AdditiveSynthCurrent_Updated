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

### 4d.6 The shifted-transient mechanism, pinned (2026-09-22)

Three measurements taken before planning the transient work. Together they overturn the
Stage-2 design that had been on the books since Sep 12.

**(a) The transient energy is NOT missing from the model — it is misplaced in time.**
Unity residual `x - x_hat` measured in the 60 ms after each onset:

| file | residual in 60 ms after onset | flatness input → residual |
|---|---|---|
| DrumLoop | **−20.7 dB** rel input (worst −10.2) | 0.031 → 0.081 |
| HappyMono | −19.3 dB (worst −12.9) | 0.166 → 0.209 |
| 1985 | −19.9 dB | 0.123 → 0.264 |

The sinusoidal model already accounts for ~99% of the energy at a drum hit. A "transient
component" that renders what the residual holds would be a −20 dB layer — it cannot fix a
+8.8 dB pre-onset error. **The planned S+T+N transient component is therefore the wrong
tool for this defect**, and the same measurement kills the cheaper variant (pasting the
unity residual unshifted at onsets): there is not enough there to matter.

Pasting the *original* unshifted is separately already rejected in code — see the R6(a)
comment at the residual mix: under shift it stamps original-pitch onsets over the shifted
tone. So neither paste is available; the fix has to be temporal.

**(b) The smear is one synthesis hop, not one analysis window.** Excess output level before
each hit, DrumLoop up5:

| ms before onset | −120 | −100 | −85 | −70 | −55 | −40 | −25 | −10 |
|---|---|---|---|---|---|---|---|---|
| Sep 21 | +1.5 | +1.0 | −1.4 | −1.3 | −0.4 | +1.2 | **+5.8** | **+8.8** |
| joint off | +2.3 | +1.8 | +0.9 | +2.7 | +1.9 | +0.8 | +2.4 | +6.4 |

Nothing before ~30 ms. `LONG_SIZE = 2048` at 48 kHz = 42.7 ms, hop 1024 = **21.3 ms** — the
smear is exactly the leading half of the synthesis frame that contains the attack. A
stationary sinusoid given the attack's amplitude fills that frame's whole span, including
the 21 ms before the hit. (Not the 4096 analysis window: that would have shown at −85 ms.)

**(c) The transient detector misses 40% of the hits.** `TRANS_DEBUG` (added here) vs the
measured onsets on DrumLoop:

```
detector: 0.299  0.619  0.960  1.259  1.451  1.600  1.920  2.261  2.560
onsets:   0.31   0.42   0.56  0.65   0.80   1.28   1.62   1.85   1.92   2.27
```

Recall 6/10 (missed 0.42, 0.56, 0.80, 1.85), precision 6/9 (0.960, 1.451, 2.560 fire with
no onset), and matched hits land up to 1.5 frames late (0.619 for a 0.65 hit). A missed
onset never reaches the short-window path at all, so it takes the full 21 ms smear. This is
the largest single lever on the defect and it is a detector fix, not a model change.

Note `SHORT_SIZE = 256` (5.3 ms), so where the detector *does* fire the frame is short
enough — but `processFrames` estimates each short frame's spectrum from a **long-window FFT
centred on it** (`m_transientLongSpecs`, added deliberately so the spectrum evolves through
the transient). Amplitude read from a 2048-sample window is smeared by construction even in
a 5 ms frame. That is the second lever.

### 4d.7 CORRECTION to 4d.6(c), and the transient rig (2026-09-22)

**4d.6(c) was wrong and is withdrawn.** It reported the transient detector at 6/10 recall
and 6/9 precision on DrumLoop. That scoring used an envelope "ground truth" that fired on
+1.7 dB and +0.7 dB steps -- not onsets. Rescored against a defensible reference (>=6 dB
rise in a 10 ms envelope, above -32 dBFS), and allowing for the fact that the detector
legitimately flags the FRAME whose 2048 window contains the hit (so it fires 9-31 ms early
by design):

```
hits:     0.31  0.65  0.97  1.14  1.28  (1.45)  1.62  1.95  2.27
detector: 0.30  0.62  0.96   --   1.26   1.45   1.60  1.92  2.26   (+ 2.56 at EOF)
```

Recall 7/8, precision 7/9 -- the detector is **fine**, and Step 1 (detection recall) is not
the lever it was billed as. One quiet hit at 1.14 is missed; that is a tuning matter, not a
rewrite. The lesson is the recurring one in this project (R3, R6, R7, §4b.5): an absolute
envelope metric invented a defect that direct inspection disproves. Ground truth first.

**The transient rig.** Busy percussion makes every envelope metric ambiguous, because
"before this hit" is also "after the previous hit" and slow release reads as pre-echo. The
rig removes the ambiguity: one hit at t = 0.75 s in digital silence (kick, snare, click,
pitched note). Any energy before the hit is pre-echo, full stop.

| rig file | unity | up5 before | up5 after | down5 before | down5 after |
|---|---|---|---|---|---|
| kick | **−240 dB** (silent) | −12.4 | **−36.1** | −9.1 | **−33.1** |
| pitched note | −240 | −10.4 | **−33.3** | −10.7 | **−33.5** |
| snare | −240 | −10.9 | −15.5 | −7.9 | −17.2 |
| click | −240 | −7.1 | −6.5 | −5.9 | −6.3 |

(dB below the hit's own peak; "after" = `1f4c810`.) Onset delay on kick/note falls from
~10 ms early to <1 ms. Unity shows the passthrough doing its job: exactly zero pre-echo,
which is why none of this was ever visible at unity.

**Reading.** The smear was the long-window amplitude, as §4d.6b argued -- but the mechanism
is specifically the *amplitude read*, not the synthesis frame length, and the fix is a
level correction rather than anything structural. What remains after it is the **click**:
a pure noise transient has no partials to carry the correction, so it is untouched. That
is the noise-model case, and it is now the only unaddressed part of the transient defect.

### 4d.8 Verification on REAL transients, and the verdict (2026-09-22)

Riley heard no difference either way on the A/B set, so "keep it" turned on whether the
synthetic-rig win transfers to real material with no cost. Test: real hits cut out of the
corpus and spliced into digital silence (0.5 ms fades) -- real material, unambiguous
ground truth.

| hit (source) | shift | pre-echo off → on | peak error off → on |
|---|---|---|---|
| kick (DrumLoop 1.615 s) | +5 / −5 | −16.0 → **−38.3** / −15.3 → **−37.3** | −2.5 → −2.3 / −3.9 → −3.0 |
| snare (DrumLoop 1.945 s) | +5 / −5 | −24.3 → **−41.2** / −26.2 → **−44.1** | −3.4 → −2.6 / −2.1 → −1.8 |
| Happy hit (0.330 s) | +5 / −5 | −25.5 → **−49.2** / −30.2 → **−53.8** | +0.1 → +0.5 / −0.7 → −0.8 |
| piano (0.180 s) | +5 / −5 | −16.9 → **−38.7** / −12.2 → **−33.5** | −0.0 → +0.1 / −0.2 → −0.4 |
| cymbal (48kCymbal 0 s) | +5 / −5 | −10.8 → −15.3 / −7.8 → −14.2 | −0.9 → −1.2 / −2.3 → −2.0 |

**16–24 dB less pre-echo on every pitched real hit, 4–6 dB on the cymbal, and the onset
peak is equal or slightly more accurate in 9 of 10 cases.** Onset delay flips from −5 ms
(early) to +4 ms, which matches what unity does (+3.9 ms) — the negative delay was the
pre-echo ramp crossing the threshold, not the attack.

**Two gating experiments, both rejected as dead weight:**
* *Frame-level onset gate* (correct only frames with a ≥6 dB rise of their own) — destroys
  the win: the frames that need correcting are the ones BEFORE the hit, which by
  definition have no rise.
* *Region-level onset gate* (correct only where the surrounding long window contains a
  ≥6 dB step) — bit-identical to no gate on the whole corpus. The "false positive"
  detections on Fairlight C3 and the Female *do* contain real level steps (the note's own
  quiet attack at −67 dBFS), so the correction was acting at genuine onsets all along. The
  earlier false-positive reading came from a 30 ms lookback that is too short for a soft
  attack. Removed rather than kept as a no-op.

**Verdict: keep.** Win confirmed on real material, unity unaffected (≤0.01 dB on all 15),
onset peaks equal or better, no audible difference reported in either direction. The one
negative metric (Fairlight `jitter_db` 0.46 → 0.72 at ±5) is an absolute trajectory metric
at a shifted condition, where the battery has no same-pitch reference and its own
documentation says the absolute forms are unreliable (§4b.5); the reference-relative
`jitter_dev_db` moves 0.095 → 0.103, and short-time level stability is unchanged
(10.69 → 10.74 dB std). What the corpus-wide difference under shift actually reflects is
phase divergence: any perturbation to a track's amplitude changes its propagated phase
from that frame on, so a legitimate correction at one quiet onset re-randomises the
waveform without changing what it sounds like.

**Why it is inaudible here:** at −16 dB under the hit and ~5 ms early, the pre-echo is
masked by the hit itself in busy material. It is a correctness fix with headroom value,
not an audible one on this corpus.

### 4d.9 Step 4 groundwork — the choir gap is per-partial FREQUENCY MOTION (2026-09-22)

**The gate: `battery/polyvoice_rig.py`.** Three steady voices, one vibrato voice, and three
voices with *independent* vibrato rates and phases (a major triad, the shape of
choir-burst). Engine today, unity, residual off:

| rig signal | amp_energy | shape_corr | srr_db |
|---|---|---|---|
| POLY_STEADY (3 steady voices) | 1.004 | 1.000 | **30.14** |
| ONE_VIB (1 voice, vibrato) | 1.000 | 1.000 | **38.09** |
| POLY_VIB (3 voices, independent vibrato) | 0.849 | 0.975 | **13.11** |

Polyphony alone is fine (30 dB). Single-voice vibrato is fine (38 dB) — and demonstrably
*because* of pitch-sync: forcing it off drops ONE_VIB to 12.91 dB, a **25 dB** swing, while
POLY_VIB is unmoved (13.11 → 13.11 auto/off, 14.53 forced). Forcing it on POLY_STEADY costs
17 dB (30.14 → 12.76), so the auto-gate is refusing correctly. **POLY_VIB at 13.11 dB
reproduces the real choir (15.5 dB) and is the gate for this step.**

**Where the 17 dB sits.** Fitting the rig's 36 *known-true* partials per frame, 2048 window,
hop 1024, and overlap-adding:

| basis | SRR |
|---|---|
| engine today | 13.11 |
| stationary cos/sin at the TRUE instantaneous frequencies | **24.39** |
| chirped: true frequency AND true df/dt | **46.18** |

So the gap decomposes into **~11 dB of frequency-estimation error** (13.11 → 24.39, what the
engine loses by mis-measuring moving partials) and **~22 dB of stationarity** (24.39 → 46.18,
the ceiling of any constant-frequency-per-frame model on this material).

**The chirp must be in BOTH the fit and the render:**

| fit / render | SRR |
|---|---|
| stationary / stationary | 24.39 |
| stationary / chirped | **1.44** |
| chirped / stationary | 20.01 |
| chirped / chirped | **46.18** |

Either half alone is *worse than neither* — the same "estimate on the basis you synthesise"
law that separated joint_mode 1 from joint_mode 2 (§4c.1), now quantified for chirps.

**It is robust to a sloppy slope.** With the chirp rate in error by 5 / 15 / 30 / 50%:
46.5 / 41.7 / 35.9 / 31.3 dB. Even a 50%-wrong slope beats stationary by 7 dB, so the slope
can come from a plain centred difference of the track's own frequency trajectory (available
post-tracking, exactly like the file-start credit's lookahead) rather than from a
higher-order estimator such as DDM.

**One cheap implementation route is already ruled out.** `joint_amp_phase` builds its normal
matrix from closed-form Hann DTFTs, which assume stationary atoms. Keeping that stationary
Gram and chirping only the correlation vector gives **−60.66 dB**: G and b must come from the
same basis. The chirped Gram has to be computed numerically for chirped pairs (O(N) each,
inside the existing band/neighbour limits), with the closed form retained for
stationary-stationary pairs so material without vibrato pays nothing.

### 4d.10 Step 4 RESULT: the choir gap is the ANALYSIS WINDOW, not the estimator (2026-09-22)

The chirped estimator of §4d.9 is built (slopes, chirped joint solve with a numeric Gram
for chirped pairs, chirped synthesis) and it **does not pay**: rig POLY_VIB 13.11 → 13.26 dB.
Chasing that down produced the real answer.

**1. It is not the slope estimate.** A centred difference over the track's trajectory
recovers the chirp rate to ~23% median relative error — well inside the ±30–50% the solve
tolerates (§4d.9).

**2. It is the frequency estimate, and no estimator fixes it.** Frequency error on the rig's
partials, 4096 window:

| estimator | rms | median |
|---|---|---|
| parabolic (what the engine uses) | 13.3 Hz | 4.0 Hz |
| reassignment / derivative method | 11.0 Hz | 5.2 Hz |
| 2-D (f, df) matched filter | 15.0 Hz | 3.3 Hz |
| joint coordinate descent, 3 rounds, interference removed | 11.6 Hz | 3.0 Hz |

Everything plateaus near 3 Hz. And the chirp needs better than 1 Hz to pay — with frequency
error of 0 / 0.5 / 1 / 2 / 5 Hz the chirped fit scores 46.2 / 36.6 / 30.8 / 25.2 / 17.1 dB
against a stationary 24.4 / 24.2 / 23.6 / 21.6 / 16.4. At the engine's error the two models
are indistinguishable, which is exactly what the engine measured.

**3. Why they plateau — the window is too long for the motion.** How far a partial's true
frequency travels *within one analysis window*, and what is left after fitting the best
straight line through it (i.e. the error a linear-FM atom cannot represent):

| window | ms | median sweep | in bins | residual after a linear fit |
|---|---|---|---|---|
| **4096 (the engine)** | 85.3 | **75.8 Hz** (max 266) | **6.5** | **8.77 Hz** |
| 2048 | 42.7 | 40.7 Hz | 1.7 | 2.09 Hz |
| 1024 | 21.3 | 21.2 Hz | 0.5 | 0.52 Hz |
| 512 | 10.7 | 10.8 Hz | 0.1 | — |

At 85 ms of 5.5 Hz vibrato the window spans ~47% of a vibrato cycle. No constant (0th
order) and no chirp (1st order) describes the partial there — 8.8 Hz of curvature is
unmodellable in principle, which is the floor every estimator above ran into.

**Conclusion.** The choir's ~13 dB is not reachable by better estimation on the current
analysis window. It needs either a **shorter window for moving partials** (where the chirp
model becomes valid — the residual falls to 0.5 Hz at 1024) at the cost of resolution and
worse inter-partial interference, or the motion removed first, which is what pitch-sync
does for one voice and would require **multi-f0 separation** for a choir. Both are larger
projects than Step 4 was scoped as, and multi-f0 separation in particular is a different
kind of engine.

This also retro-explains [[project_voice_reassignment]]: reassignment failed on vibrato in
Aug not because the method is wrong but because at 4096 there is no instantaneous frequency
to find — the note in that memory, "freq already accurate, window too long vs vibrato", was
right and is now quantified.

**Disposition.** The chirp path is kept, `chirp_mode = 0` (env `CHIRP`), threshold
300 Hz/s (measured: at 20 Hz/s estimator jitter chirps *steady* partials and costs the rig's
steady case 5.6 dB). Output with it off is byte-identical to before. It is the validated
half of a short-window/chirp pair and should be switched on only with that partner.

### 4d.11 Steps 5 and 6: the same wall as Step 4 — analysis resolution (2026-09-23)

**Fairlight's residual is not junk partials.** The hypothesis from §4d.4 (199 sub-40 dB
partials diluting the joint solve) is wrong: **98% of the residual sits AT partials**, 2%
between them. And the loud partials are rendered well — per-partial error over 193 frames
on C2: amplitude 0.3–1.0 dB mean, phase 0–1°. Ranking partials by residual contribution:

| freq | level | share of residual | local SRR |
|---|---|---|---|
| 64.5 Hz | 0 dB | **30.6%** | 34.2 dB |
| 1892.6 | −38.9 | 11.0% | 5.7 |
| 1699.2 | −23.9 | 8.5% | 17.5 |
| 893.6 | −45.1 | 8.0% | 4.3 |
| 155.3 / 105.5 | −51 / −53 | 1.9 / 1.3% | 0.4 / −0.4 |

Two populations: the fundamental, which is simply so loud that a 34 dB local error still
dominates absolutely, and a tail of −39 to −53 dB partials reproduced at 0–6 dB. The weak
ones sit in **clusters spaced ~17.6 Hz** (1640.6, 1658.2, 1675.8, 1699.2, 1728.5, 1740.2),
far inside the 47 Hz Hann main lobe at 4096. Unresolvable.

**But a longer window does not fix them.** An oracle stationary fit at increasing window
lengths (peaks detected and amp/phase solved at that window, OLA'd, scored):

| file | engine | N=4096 | N=8192 | N=16384 |
|---|---|---|---|---|
| Fairlight C2 | 26.8 | 25.2 | 21.2 | 21.9 |
| Fairlight C3 | 23.9 | 22.9 | 19.6 | 14.9 |
| Piano | 25.3 | 20.7 | **25.5** | 23.0 |
| SaintSaëns | 22.0 | 12.7 | 13.7 | **16.1** |
| choir | 15.5 | 12.2 | 9.2 | 8.8 |

Fairlight *loses* 3–8 dB with a longer window, because its partials carry ~25 dB of
amplitude modulation at 0.3–2.8 Hz plus that 17.6 Hz structure: resolve the partials better
and you violate constant-amplitude worse. It is the same trade Step 4 hit, in amplitude
rather than frequency. Raising the engine's own 16384 LF tier above its 200 Hz cutoff
confirms it in the engine: sweeping `LF_CUTOFF` to 400/800/1600 Hz gives at best +0.9 dB
(1985) and +0.7 (Happy) against −1.5 (Fairlight C2), −1.2 (Female) and −1.3 (SaintSaëns).

**What an adaptive window would actually buy.** Residual per band at each window, taking the
best window per band:

| file | 20–200 | 200–800 | 800–3k | 3k–12k | whole (4096) | whole (per-band best) |
|---|---|---|---|---|---|---|
| Fairlight C2 | 40.7 @2048 | 26.4 @4096 | 23.6 @2048 | 18.9 @4096 | 31.9 | **33.3** |
| SaintSaëns | 14.5 @8192 | 16.9 @4096 | 14.3 @2048 | 9.6 @1024 | 12.8 | **14.8** |
| choir | 15.8 @4096 | **19.2 @2048** (vs 14.6 @4096) | 11.8 @2048 | 8.3 @1024 | 12.2 | **15.6** |

The best window differs per band *and* per file, and no single choice wins: +1.4 dB
(Fairlight), +2.0 (SaintSaëns), +3.4 (choir). Note choir's 200–800 Hz band gaining **4.6 dB**
from a *shorter* window — §4d.10's diagnosis showing up independently.

**The unifying conclusion.** Every remaining per-class gap — choir 13 dB, Fairlight 11,
Piano/SaintSaëns 7–8 — is an analysis time-frequency limit, not an estimator defect. The
material modulates (in frequency, in amplitude, or both) within any window long enough to
resolve it. Multi-resolution analysis recovers **2–3 dB** of it; the remainder needs atoms
that model modulation *and* a way to estimate their parameters under interference, which is
a different analysis architecture, not a fix.

So the engine is close to the practical ceiling of a single-window STFT front end with
constant-amplitude, constant-frequency atoms per frame. Stage 0's oracle numbers (38–65 dB)
remain true and remain out of reach: they come from a per-frame free-frequency fit with no
tracking, no continuity and no decomposability — not something this engine can become.

## 5-REVISED. The plan (2026-09-22)

Supersedes §5 below, which is kept for its record of what was tried. Reordered by §4d.4's
ceiling table and §4d.6's mechanism findings. Each step is independently shippable,
ear-gated, and revertable behind an env knob.

### Step 1 — Transient detection recall — DEMOTED (see §4d.7)

Rescoring against a correct ground truth puts the detector at 7/8 recall, 7/9 precision.
One quiet hit missed on DrumLoop. Worth a threshold tweak eventually; not a lever.
Original text kept below for the candidate fixes, which still apply if it is revisited.

#### (original framing, superseded)

Recall 6/10 and precision 6/9 on DrumLoop (§4d.6c). Every missed hit is rendered by a long
frame and takes the full 21 ms pre-smear; every false positive spends short frames where
they are not needed. `transientNegotiationTactics` sums positive per-bin dB differences
across the whole spectrum, divides by bin count, compares to one global threshold, and
suppresses any detection within 2 frames of the previous one.

Candidate fixes, cheapest first: per-band flux (a kick's energy is <200 Hz; a full-spectrum
mean dilutes it), a threshold relative to a running median rather than a constant,
half-hop resolution for placement, and replacing the blanket 2-frame suppression with an
energy-rise test so consecutive real hits survive.

**Gate:** recall ≥ 9/10 and precision ≥ 8/10 on DrumLoop; no regression on the
tonal files' transient counts (a false positive on the Female costs a short-window region
in the middle of a vowel); shifted pre-onset excess drops.

### Step 2 — Short-frame amplitude from a short window — DONE and KEPT (`1f4c810`, §4d.8)

Where the detector does fire, the amplitude still comes from a 2048-sample FFT centred on
the 256-sample frame, so it is smeared by construction (§4d.6). Estimate amplitude for
transient short frames from a short-window magnitude (or from the frame's own time-domain
energy), keeping the long-window FFT for *frequency*, which is what it was added for.

**Gate:** shifted DrumLoop pre-onset mean < +2 dB (from +7.7), worst < +6 (from +19);
attack delay < 2 ms (from 4–6); peaks within ±3 dB; unity byte-identical (this path is
shift-relevant but runs at unity too — expect a small unity change and ear-check it).

### Step 3 — Envelope correction under shift, transient regions only (fallback)

Only if Steps 1–2 leave audible smear. The sinusoids have the right spectrum and the wrong
temporal envelope, so apply a corrective gain in transient regions: target = the input's
short-time envelope, time-translated (transients translate, they do not transpose). Gain
only — it cannot invent content and cannot bleed original pitch, which is what killed both
paste variants (§4d.6a). Risk is the usual one for gating a smeared signal: too aggressive
and the gate itself is audible.

### Step 4 — choir — BLOCKED, and demoted (see §4d.10)

The gap is the 85 ms analysis window, not the estimator: a vibrato partial sweeps 6.5 bins
across it and retains 8.8 Hz of curvature after the best linear fit, so no 0th- or 1st-order
model can represent it there. Reopening this means a shorter window for moving partials, or
multi-f0 separation so each voice can be warped the way pitch-sync warps one. Both are
bigger than this step was scoped as. The chirped solve and synthesis are built and gated
(`CHIRP=1`) awaiting the short-window partner.

#### (original framing, superseded)

> Groundwork done, §4d.9: gate built, gap decomposed (11 dB frequency + 22 dB stationarity),
> chirp shown to need fit AND render, slope shown to be robust to 30–50% error, and the
> cheap stationary-Gram shortcut disproved. What remains is the C++ implementation.

~13 dB below ceiling, residual still tonal. Many voices with independent vibrato, so
pitch-sync correctly refuses to engage (no single f0). Demodulate **each track** against its
own smoothed frequency trajectory before estimating amplitude and phase — the
generalisation of the warp that fixed the Female. Extend `vibrato_rig.py` to a two-voice
rig first so there is a gate before any engine code.

**Gate:** rig two-voice shape_corr; choir residual_srr 15.5 → >22; Female not regressed.

### Step 5 / Step 6 — CLOSED by §4d.11

Fairlight's gap is closely-spaced partials (~17.6 Hz clusters) carrying 25 dB of amplitude
modulation — resolving them needs a longer window, representing them needs a shorter one.
Piano/SaintSaëns are the same trade. Neither is an estimator problem. What remains available
is a per-band adaptive window worth 2–3 dB (see Step 7).

#### (original framing, superseded)

The oracle reaches 38 dB with 16 partials where the engine renders 213 for 27 dB. A
synthetic wavetable, so the truth is knowable: dump the engine's tracks against the file's
actual harmonic series and find out whether the loss is spectral motion within the frame
(shares a fix with Step 4), the 40 dB-down junk partials diluting the joint solve, or
something specific to wavetable crossfades.

### Step 6 — Piano / SaintSaëns (~7–8 dB): per-band high-resolution estimation

Closely-spaced partials below the Fourier limit (separation < 2.28·Fs/M). ESPRIT or
matrix-pencil per band, gated by model-order selection, applied only where a pair is
provably unresolved. Highest machinery-per-dB in the list; do it last.

### Step 7 — the noise/click transient under shift — DONE (§4d.12), awaiting ear

The one remaining defect that is *audible* rather than only measurable, and the only part of
the shifted-transient problem `1f4c810` could not touch. Rig pre-echo after that fix: kick
−36 dB, pitched note −33, snare −16, **click −6.5 (unchanged)**. A pure noise transient has
no partials for a level correction to act on; it is represented by the stochastic residual,
which under shift is a magnitude-deficit fill above `residual_hp_hz_shift`, with no temporal
structure of its own.

The fix is the noise half of S+T+N, scoped to onsets: fit the residual's per-band energy
envelope at short resolution through the transient and render band-limited noise that
follows it, unshifted (noise does not transpose). This is also the cymbal case — the one
class Stage 0 called structurally limited (engine 10.9 vs 13.1 at K256).

**Gate:** `battery/transient_rig.py` click pre-echo −6.5 → below −20 dB with the hit's own
peak within 2 dB; 48kCymbal/Out48k shifted ridge glide toward the input's 2208 Hz/step;
DrumLoop and Happy not regressed; tonal files byte-identical (no transient regions). Then ears.

### 4d.12 Step 7 DONE — the click's pre-echo was the residual, not the model (2026-09-23)

`1f4c810` fixed the tonal half of the shifted transient and left the click at −6.5 dB. The
cause turned out not to be the sinusoidal model at all. Rendering the rig click with the
stochastic fill switched off:

| | pre-echo |
|---|---|
| full engine, up5 | −6.5 dB |
| **residual fill off**, up5 | **−23.8 dB** |

So the tonal model was already clean (the amplitude fix did its job) and the pre-echo is the
**noise fill**: it is built from 1024-sample (21.3 ms) STFT frames, which spread a click's
energy ±10 ms — measured as starting 15 ms before the hit.

**Fix (`residual_env_cap`, default on):** noise carries no phase structure worth protecting,
so shape it in the time domain instead of chasing resolution. A short moving-RMS compares
the fill against the input; wherever the fill would carry more local energy than the input
does, it is scaled down. A cap, never a boost, so it can only remove energy the input does
not support.

Rig, pre-echo dB below the hit's own peak:

| rig | up5 before → after | down5 before → after | peak error |
|---|---|---|---|
| click | −6.5 → **−31.0** | −6.3 → **−31.5** | −7.4 → −5.8 (better) |
| snare | −15.5 → **−35.5** | −17.2 → **−33.2** | unchanged |
| kick | −36.1 → −36.2 | −33.1 → −33.2 | unchanged |
| note | −33.3 → −33.4 | −33.5 → −33.6 | unchanged |

Real hits spliced into silence:

| hit | up5 before → after | down5 before → after |
|---|---|---|
| **cymbal** | −15.3 → **−32.6** | −14.2 → **−32.7** |
| drum snare | −41.2 → −46.8 | −44.1 → −48.0 |
| drum kick | −38.3 → −41.5 | −37.3 → −41.6 |
| piano | −38.7 → −40.5 | −33.5 → −36.2 |

Onset peak error is unchanged on every one. **Unity is unchanged to 0.00 dB on all 15
files** — the cap only fires where the fill exceeds the input, which at unity it never does.

**Battery: zero metrics changed by more than 5%.** It has no shifted-condition pre-echo
metric, so it is blind to this in the same way it was blind to `1f4c810`; the rig and the
spliced-hit test are the instruments that see it.

**The cymbal chirp is NOT fixed by this** and was the wrong gate to have set: upper-band
centroid drift per frame is 139 Hz before and 139 after (input 197). The cap constrains the
fill in time, not the tonalised tracks' frequencies under shift, which is what the chirp is.
That remains open and is a tracking/representation issue, not a noise-envelope one.

### 4d.13 What else is below the ceiling — a full ablation of the non-analysis stages (2026-09-23)

§4d.11 concluded the remaining per-class gaps are analysis limits. That invites the obvious
question: is anything *outside* the analysis also costing fidelity? Ablating each stage at
unity, mean Δ residual SRR over 10 files:

| stage ablated | mean Δ | reading |
|---|---|---|
| **amplitude EMA removed** | **+2.29 dB** | **the only positive — a real cost** |
| joint solve off (legacy per-peak) | −11.28 | confirms §4c |
| MQ oscillator bank instead of OLA | −2.58 | OLA is the better synthesis path; MQ costs 5–6 dB on dense mixes |
| unity dedup on | −1.26 | correctly off by default |
| stochastic residual off | +0.07 | SRR-neutral at unity; it exists for perceptual fill |

**One lever, and it is the amplitude EMA.** The joint solve's output is passed through the
tracker's asymmetric amplitude EMA (added in `5bbbf37` as variance control after the raw
solve made amplitudes jumpy). It is a low-pass on the amplitude trajectory and cannot
distinguish estimator noise from real modulation: removing it gains Fairlight C2 +5.3 dB,
Piano +4.5, take-me-out +3.3, Happy +3.3.

**Median-of-3 does not transfer from frequency to amplitude** (−4.3 dB mean). `e5c2b41` used
exactly that filter on the frequency trajectory for exactly this reason — kill spikes, pass
ramps — but amplitude modulates far faster than frequency does: Fairlight carries ~17.6 Hz
AM against a 46.9 Hz frame rate, under three frames per cycle, so a median destroys the
modulation it is supposed to preserve. A level-gated EMA (raw for partials within 20 dB of
the frame's loudest) recovers +1.74 of the +2.0 but does not buy back the variance.

**The battery on raw vs EMA**, scoring `residual_flat_ratio` higher-is-better as `run.py`
does: **113 better / 46 worse**. The split matters —

| family | better / worse |
|---|---|
| residual_srr | **64 / 3** |
| residual_flat_ratio | 7 / 1 |
| env_p2p_dev (reference-relative) | 13 / 10 |
| jitter_dev (reference-relative) | 7 / 3 |
| env_p2p ABSOLUTE | 11 / 18 |
| jitter ABSOLUTE | 8 / 10 |

The reference-relative forms are break-even to better; only the absolute forms regress, and
those are the ones `8df4106` documents as unable to tell "the engine over-smoothed" from
"the input really varies". `shape_consistency` — the metric built for the saw wobble
complaint — does not move at all.

**Left at the EMA anyway, pending ears.** The cost concentrates on the voice (Female
`env_p2p_full` 7.93 → 8.30 at unity, 4.57 → 5.62 at +5) and that is precisely the "shaky
voice" class `e5c2b41` was written to fix. +2 dB against a possible return of a
previously-fixed artifact is an ear decision. `JOINT_SMOOTH=0`; A/B in
`~/Desktop/comparingOut/sep23_ampEMA_AB/`.

**Everything else outside the analysis is clean.** Synthesis (OLA), the residual, dedup and
the limiter cost nothing measurable at unity; MQ is worse and is already not the default.
So §4d.11 stands, with this one amendment: there is ~2 dB available in the amplitude
smoothing policy, and it is a *policy* choice rather than an analysis limit.

### 4d.14 Robustness and performance survey (2026-09-23)

Before planning further fidelity work, a check for defects that fidelity metrics would never
show.

**Performance is fine, and my earlier "10-minute render" flag was wrong.** Release build:
0.3–1.3× realtime across the corpus, shifted or not. The 10-minute case was the *Debug*
binary that `render.sh` drives through Xcode DerivedData. No perf work is warranted.

**Extreme shifts hold up.** ±12 and ±24 semitones on five files: no non-finite samples,
peaks bounded by the limiter, output level within ~2 dB of input everywhere. We had only
ever tested ±7.

**Edge-case inputs hold up.** Silence, DC offset, white noise, a clipped sine, −84 dBFS
material, a 50 ms file, a file whose onset is at 2.5 s, DC+tone, and a 20 Hz→20 kHz sweep:
no crashes, no non-finite output, nothing above full scale.

Two notes from that sweep, neither a bug:
* A pure DC input renders silent under shift (DC is below `lf_min_hz`, and a shifted DC is
  still DC). Removing it is the right behaviour.
* **A full-range sweep loses 5.9 dB** (−7.5 dB in, −13.4 out at unity). A 20 Hz→20 kHz sweep
  in 2 s moves at ~10 kHz/s, which is the §4d.10 window limit at an extreme rate. It is the
  same wall, confirmed from a third direction.

**Conclusion: the engine is robust.** There is no hidden defect class here. What is left is
what the analysis cannot represent, and what the corpus cannot tell us.

## 5b. BLIND LISTENING TEST RESULTS (2026-09-24)

15 clips, each original vs engine resynthesis at unity, randomised A/B, RMS level-matched,
sent out blind. **11 "can't tell". 4 identified — and all 4 correctly**, with a specific
symptom named for each. Set: `~/Desktop/comparingOut/sep23_listening/`.

| clip | source | verdict | symptom named by the listener |
|---|---|---|---|
| 1 | Female sung line | correct | "a double hit in the first utterance… a clear tell and makes it obvious of a bug" |
| 9 | Fairlight C3 | correct | "A has some sort of onset that B is not including before the note starts" |
| 10 | Fairlight C2 | correct | "content missing… stutters in the original… too smooth, like we sanded out the fine details" |
| 11 | choir | correct (hedged) | "the choir maybe isn't as full" |

### 5b.1 The meta-finding: SRR does not predict audibility

| group | mean SRR | mean p10 local SRR | mean worst 20 ms block | mean spectral flatness |
|---|---|---|---|---|
| **detected (4)** | **21.3 dB** | 14.6 | **1.9** | **0.0005** |
| not detected (11) | **21.5 dB** | 21.2 | 8.9 | 0.0061 |

Average SRR is **identical** across the two groups. Fairlight C2 at 26.8 dB was identified;
both cymbals at 8.8 and 10.2 dB were not. What separates them:

* **Exposure.** The detected files are 12× more tonal (flatness 0.0005 vs 0.0061). On sparse
  tonal material an error is naked; on cymbals and drums it is masked.
* **Worst case, not average.** Detected files average a worst-block local SRR of 1.9 dB
  against 8.9 for the rest.

**The battery's whole-file `residual_srr_db` is therefore the wrong gate for perceived
quality**, which is why four consecutive rounds of measurable improvement were inaudible.
A perceptually aligned gate would be worst-case (p10 / min) local SRR, reported separately
for tonal material. That is a battery change, and it is cheap.

### 5b.2 Clip 1, Female — a REGRESSION I introduced, and the clearest defect found

First 10 ms blocks, dBFS:

| t (ms) | 0 | 10 | 20 | 30 | 40 |
|---|---|---|---|---|---|
| input | −83.7 | −45.7 | −29.5 | −18.3 | −16.7 |
| engine (file-start credit ON) | **−28.3** | **−27.4** | −26.7 | −20.8 | −18.6 |
| engine (credit OFF) | −87.7 | −55.6 | −34.1 | −21.1 | −18.1 |

The engine emits **−28 dBFS where the input is silent**, 30 ms before the real onset: a
55 dB error, and audibly a second, earlier hit. That is the pre-echo flagged as a known
side effect when the file-start credit landed (`bbc8d6a`) — now confirmed audible, and the
single worst defect in the test. `Female` is also the only file whose worst 20 ms block has
a *negative* local SRR (−20.5 dB), at t=0.

Scope: 4 of 15 files carry it, in every case one whose onset is delayed —
Female **+20.9 dB**, 1985 +8.0, take-me-out +7.1, SaintSaëns +6.5. Only Female was heard,
because it is an exposed solo voice while the others are dense.

### 5b.3 Clips 9 and 10, both Fairlights — spectral change without level change

Half of each Fairlight's residual lives in **6% of the time** (C2: 12 of 198 20 ms blocks;
C3: 8 of 141), at discrete moments — C2 at 0.72, 2.00, 2.52 s; C3 at 0.52, 1.42 s. At those
moments:

* the **level is flat** (±1 dB across the event), and
* the **spectrum changes 6–10 dB per bin**.

That is a wavetable switch: the waveform's shape changes abruptly at constant loudness.
`transientNegotiationTactics` sums *positive* per-bin dB differences against one global
threshold — it is an onset/energy-rise detector — so it fires **once** on C2 (at 0.064 s)
and three times on C3, missing all of these. The events are therefore analysed with 85 ms
windows and smeared into gradual crossfades. "Too smooth, like we sanded out the fine
details" is a literal description of that.

Ruled out along the way, all measured: per-partial AM is preserved (correlation
0.98–0.999, depth within 1.5 dB); whole-file band levels are within 1 dB; envelope
modulation from 1–400 Hz is within 1 dB; beating inside the close-spaced 17.6 Hz clusters
survives (p2p within 1–2 dB, rates identical). The defect is **not** smoothing, and not the
amplitude EMA — it is event detection.

Clip 9's "missing onset before the note" has a second component: through the quiet pre-note
passage C3's engine output runs **4–6 dB down above 4 kHz** (input −33.7 dB, engine −38.7)
— the noise floor and room tone ahead of the note are under-filled.

### 5b.4 Clip 11, choir — mis-placed partials, not missing ones

Partial *count* is right (input 80 above −60 dB, engine 82). But **10 of the 80 input
partials have no engine partial within 10 Hz**, and they carry **9.2% of the partial
energy** — one of them only 10 dB below the loudest. Energy rendered at the wrong frequency
both removes the right partial and adds a wrong one, which thins the sound: "not as full".

This is §4d.10 heard rather than measured — frequency estimation under polyphonic vibrato,
blocked on the 85 ms analysis window. It is the one detected defect with no cheap fix.

### 5b.5 What the four findings point to, in order

1. **Gate the file-start credit on input energy** (5b.2). A regression, the loudest defect,
   and the cheapest fix: the credit should not release a frame-0 birth into a span where the
   input has no energy. The residual envelope cap (§4d.12) is the same idea already proven
   on the noise path — an output may not carry energy the input does not support. Expect it
   to fix Female outright and improve three other files.
2. **Detect spectral-change events, not just onsets** (5b.3). A symmetric spectral-distance
   measure alongside the existing positive-flux sum, feeding the same short-window path.
   Addresses both Fairlights, and plausibly any wavetable, granular or heavily-modulated
   source — a class the corpus barely represents.
3. **Change the battery's gate to worst-case local SRR** (5b.1), reported separately for
   tonal material. Without it we cannot see the defects that listeners actually hear.
4. **Under-filled noise floor in quiet passages** (5b.3) — smaller, and possibly a
   consequence of 1 and 2.
5. **Choir frequency accuracy** — unchanged from §4d.10: architectural, and the only one of
   the four that Phase D would be needed for.

Nothing here argues for the Phase D rebuild. Three of the four audible defects are ordinary
bugs in event handling, and none of them is the analysis ceiling.

## 5c. THE UP-SHIFT REPORT — I could not find an objective defect (2026-09-24)

Round-2 feedback: unity is now indistinguishable, **down-shift is good on everything**, and
the remaining problem is **up-shift** — "slight gaps in the audio" on the drum loop, and on
the piano "amplitude modulation or harmonic distortion".

Pitch shift here preserves timing, so the input's own envelope IS a valid reference for
shifted output. That gives real instruments for the first time. I ran eight of them. **Every
one says up-shift is equal to or better than down-shift.**

| test | result |
|---|---|
| output limiter pinning | not engaged on any level-matched render |
| Nyquist drop / anti-alias fade (up-only by construction) | DrumLoop has **0.00%** of its energy above 14.4 kHz, so nothing reaches the fade |
| envelope fidelity, 5 ms blocks | DrumLoop dips >6 dB on **9.5%** of blocks up vs **15.3%** down |
| envelope fidelity, **transposed** bands | DrumLoop 10–15 kHz dropout **0.0%** up vs **14.3%** down |
| frame-rate (46.9 Hz) AM | no up-specific excess on any file |
| roughness — excess 20–150 Hz modulation vs input | up is **better by 1–17 dB** on all five files tested |
| unshifted-residual share of output | −35.8 dB on DrumLoop, and symmetric up/down |
| per-partial shift-ratio accuracy | +0.3 cents error, 0.8–1.5 cents scatter; identical both ways |

### 5c.1 Two methodological traps, recorded so they are not repeated

**Band comparisons must transpose.** Comparing an output band [lo, hi] against the input's
[lo, hi] is wrong for shifted audio: the content has moved. My first pass that way reported
"up-only HF dropout, 13.2% on DrumLoop" — the correct comparison against
[lo×ratio, hi×ratio] turns that into **0.0% up and 14.3% down**, i.e. the opposite sign.
Any future shifted-audio measurement must transpose the analysis band.

**A 7-semitone shift is confounded for harmonic material.** The ratio is ~3:2, so a shifted
harmonic lands on another original harmonic (2 × 1.498 ≈ 3). Any "is there energy at the
original pitch" test therefore cannot separate a ghost from the interval itself at ±5
semitones — the piano's apparent 4 dB of original-pitch excess under up-shift is not
evidence of anything. Use a non-just interval (e.g. +6 or +8 semitones) for ghost tests.

### 5c.2 Best interpretation: the artifacts are not worse, they are more exposed

Given eight unweighted measures that all favour up-shift and two independent listeners who
hear the reverse, the likeliest explanation is that our instruments are unweighted and
hearing is not:

1. **Up-shift moves the engine's error into 2–5 kHz**, where the equal-loudness contour
   peaks and the ear is most sensitive. An error that is objectively identical is
   subjectively louder there. Down-shift moves the same error away from that region.
2. **Modulation rates scale with the shift.** A residual error that beats at 30 Hz at unity
   beats at 45 Hz up (squarely in the roughness band, heard as distortion) and 20 Hz down
   (heard as flutter, far more forgiving). This predicts the piano report specifically —
   "amplitude modulation or harmonic distortion" is what roughness sounds like — even though
   the *measured* roughness excess is lower up than down.
3. **The noise/tonal relationship inverts.** The residual stays at the original pitch by
   design (`4d85b5f`: noise does not transpose). Down-shift therefore leaves the air *above*
   the harmonics, which is how natural sounds are built; up-shift moves the harmonics up
   *past* the noise, which is not. The residual is only −36 dB on DrumLoop but **−18 dB on
   the cymbals**, and a drum loop's top end *is* cymbals — the one place this would be heard.

Note what this implies: there may be no bug to fix. "It can be fixed" is a reasonable
instinct, but if the cause is exposure rather than magnitude then the remedy is *less total
error* — the analysis rebuild — not a targeted patch. I do not yet have evidence either way,
and that is the gap the plan below closes.

## 6b. PLAN — build the missing instrument, then test the three hypotheses (2026-09-24)

### B1. A perceptually weighted residual metric (days, no engine change) ← FIRST

Every metric in this project is unweighted, which is precisely why four rounds of
improvement were inaudible (§5b.1) and why up-shift now measures better than it sounds.
Weight the residual before scoring it: group into Bark bands, apply an equal-loudness
weighting, and compute SRR per band against the input's masked threshold rather than its raw
level. Report it alongside `residual_srr_worst`.

**Gate:** the new metric must rank the round-1 clips in the order the listener did (the four
identified worse than the eleven not), and must rank up-shift worse than down-shift on
DrumLoop and Piano. If it cannot reproduce judgements we already have, it is not the
instrument. That is a real falsification test, and it is worth more than any fix.

### B2. A/B the three up-shift hypotheses (cheap, decisive, no reference needed)

Render each file three ways under up-shift and send blind:
* **(a) residual transposed** with the tonal content instead of left at the original pitch;
* **(b) residual removed** entirely;
* **(c) current** behaviour.
If (a) or (b) fixes the drum loop by ear, hypothesis 3 is the cause and the fix is a knob:
transpose the noise, or transpose the part of it that belongs to the source rather than the
room. If neither changes anything, hypothesis 3 is dead and exposure (1 and 2) is the answer.

### B3. The round-trip test — a true reference for shifted audio, finally

Shift up 5 semitones, then shift the result back down 5, and blind-A/B against the original.
Two shifts compound the error, so it is a harder test than one, but for the first time it has
a **real reference**. It also directly measures the thing the report is about: whether the
up-shift path is worse than the down-shift path, since a round trip through both should
return to where it started.

Needs only a script that renders twice; no engine change.

### B4. Only then decide

If B1 gives a weighted metric that agrees with ears, and B2 finds a cause, fix it. If B2
finds nothing and B1 says up-shift error is merely better-exposed, the honest report is that
up-shift is as good as this architecture makes it, and further gain needs §4d.11's rebuild.

## 5d. THE PITCH-SHIFT A/B RESULT — my perceptual metric was wrong (2026-09-25)

The first blind pitch-shift test (`sep24_shift_AB`, current vs `JOINT_SHIFT=1`) came back
**unanimously for the shipped engine** on all three files the metric predicted would improve
most. The listener's words: the `JOINT_SHIFT=1` version "seems to have tremolo or a pitch
that can't hold itself and warbles too much, whereas the current is cleaner sounding."

That is a clean, prospective falsification. It reverses §5c/6b's conclusion.

### 5d.1 The metric that predicted correctly was the one I dismissed

`env_p2p` is envelope peak-to-peak — **tremolo**. `jitter_db` is trajectory jitter — **a
pitch that can't hold itself and warbles**. Both are the *absolute*, output-only forms I had
written off as unreliable at shifted conditions. On the three judged files:

| file | semi | metric | current | JOINT_SHIFT=1 |
|---|---|---|---|---|
| DrumLoop | −5 | env_p2p_full | 8.84 | **10.76** |
| DrumLoop | −5 | jitter_db | 1.64 | **1.84** |
| HappyMono | +5 | env_p2p_full | 4.03 | **4.47** |
| HappyMono | +5 | jitter_db | 1.66 | **1.78** |
| Female | +5 | env_p2p_full | 4.57 | **5.31** |

`env_p2p_full` is worse in all six file/direction combinations and `jitter_db` in four of
six. The battery's verdict — 15 better / 32 worse, entirely in those two families — was
**right**, and my dismissal of it was wrong.

### 5d.2 Why the perceptual metric failed, precisely

`shifted_nmr` compares Bark spectra **frame by frame** against a transposed target with a
single global gain removed. It is therefore **blind to frame-to-frame instability**: a signal
that wobbles in amplitude and frequency can match the target spectrum well in every
individual frame while warbling audibly. The artifacts that dominate shifted audio are
*temporal*, and the metric is *per-frame spectral*. Masking and absolute threshold do not
help with that — they were the right additions to the wrong measure.

**The deeper lesson is about validation, not signal processing.** Gates 1 and 2 were
retrospective: I tuned the measure until it reproduced two judgements I already had. It then
failed its first genuine prediction. A metric fitted to known answers is not validated by
them. Any future metric must be gated on a **prospective** test — predict, then listen.

### 5d.3 What the joint solve is actually doing under shift

My mechanism (amplitudes fitted for relative phases that synthesis discards) may be real, but
it is swamped by the opposite effect: **the joint solve is variance reduction as well as bias
reduction.** Its amplitudes are far more stable frame-to-frame than per-peak reads, and under
shift — where propagated phase integrates every error — stability matters more than
per-frame spectral accuracy. Removing it returns amplitudes to the noisier per-peak estimate,
and propagation turns that noise into tremolo and warble.

This also explains the earlier ablation result that I under-weighted: switching the amplitude
EMA off *raised* shifted NMR (+0.4 dB) even as it gained 2.3 dB of unweighted SRR at unity.
Under shift, smoothing helps. That is now confirmed by ear in the strongest way available.

**`JOINT_SHIFT` stays at 0 and should be treated as a closed question**, kept only as a
documented negative result.

### 5d.4 Where this leaves the up-shift complaint

Still real — two listeners independently — and still unlocalized. What is now excluded:

* noise placement (§5c, residual is −18..−36 dB, three placements within 0.1 dB);
* removing the joint solve (this test — makes it audibly worse);
* gross added modulation in the shipped build: measured against the input's own envelope
  (timing is preserved, so excess is modulation the engine added), the current engine's
  shifted envelope peak-to-peak is **at or below the input's** on nearly every file and band
  — DrumLoop +0.3/−0.7 dB, Female −0.6/+0.5, HappyMono −1.9/−2.8, cymbal +0.4/−2.4. There is
  no gross tremolo to remove;
* Nyquist/anti-alias, limiter, pitch accuracy, roughness (§5c).

So the shipped engine is not obviously doing anything wrong under up-shift, and I have no
instrument that localizes what the listeners hear. That is the bottleneck now — not a missing
fix, a missing observation.

## 6c. PLAN — get better observations, not better metrics (2026-09-25)

### C1. Localized listening (highest information, needs testers not code)

Every defect this project has actually fixed was fixed because a listener named a *place*:
"the first utterance", "before the note starts", "the stutters in the original". The
up-shift complaint has no place attached — "slight gaps" and "amplitude modulation" could be
anywhere. Cut the up-shifted drum loop and piano into 1–2 second excerpts, numbered, and ask
which excerpts contain the problem. With timestamps I can analyze those moments directly, the
way §5b.3's wavetable switches were found from "too smooth".

### C2. The external anchor (Audacity) — settles the question of scale

`sep24_shift_AB/add_audacity_here/`. If the engine's up-shift is *preferred* over a
conventional phase-vocoder shift, the complaint is about pitch shifting in general and the
engine is at a good place. If Audacity wins, there is a concrete target and a difference to
analyze. Either answer is decisive and it costs minutes. My own phase-vocoder reference was
not accurate enough to serve (a −5 shift came out at a 0.98 frequency ratio), so this needs
real Audacity renders.

### C3. Metric policy, corrected

* **Unity**: reference-relative forms (`env_p2p_dev`, `jitter_dev`) plus
  `residual_srr_worst`. Unchanged.
* **Shifted**: the **absolute** `env_p2p_full` and `jitter_db` are the gate. Their known flaw
  — "less variation always scores better" — is a flaw at unity, where the input's own
  variation is the truth. Under shift, propagated phase can only *add* modulation, so less
  really is better. They just predicted a listening result that a purpose-built perceptual
  model got backwards.
* `shifted_nmr` is retained in `battery/probes/` as a documented failure, not a gate.

### C4. Do not

* Build another perceptual metric before C1/C2. The gap is observations, not models.
* Revisit `JOINT_SHIFT`, noise placement, or the Nyquist path — all closed by measurement or
  by ear.

## 5e. WHAT AUDACITY ACTUALLY DOES, AND WHY OUR SHAPE FIX NEVER FIRES (2026-09-25)

### 5e.1 Correcting the premise: ChangePitch.cpp contains no pitch shifter

The supplied `ChangePitch.cpp` is GUI and parameter arithmetic — note/octave/cents/percent
conversion, `DeduceFrequencies()` for the dialog's "Estimated Start Pitch", and validators.
The only DSP lines are the delegation:

```
mSoundTouch->setPitchSemiTones((float)(m_dSemitonesChange));
return EffectSoundTouch::Process();
```

So the algorithm is **SoundTouch**: a **WSOLA time-stretch** to 2^(n/12)× the length followed
by **resampling** back. (The checkbox switches to SBSMS, a subband sinusoidal model that also
pitch-shifts by resampling *after* sinusoidal time-stretching.)

This matters more than it looks. WSOLA never models the sound: it splices waveform segments
chosen for local similarity. It therefore has **no partials, no tracks, no phase propagation
and no parameter estimation** — and so none of our artifact classes. That is why it is robust
on a drum loop. Its own costs are different (transient doubling, and it transposes the
spectral envelope with the harmonics, which is the "chipmunk" colour).

So Audacity is not a template for this engine. It is a demonstration that a method with no
model has no model errors — useful as an A/B anchor, not as a design to copy.

### 5e.2 The finding: our shape-invariant machinery almost never engages

The harmonic lock (`1e85894`, `33ef83b`) exists precisely to stop relative-phase drift under
shift. `HLOCK_STATS` (added here) counts how many track-frames actually render under it:

| file | track-frames locked, up5 |
|---|---|
| 300hzSaw (what it was built and validated on) | **69.1%** |
| HappyMono | 14.6% |
| 440sawtooth | 8.8% |
| PianoSampleMono | 5.4% |
| Female sung line | 2.7% |
| **DrumLoopShort** | **2.5%** |

On real material it fires on 2–6% of partials. **97.5% of the drum loop's partials are
free-running per-track phase integration** — exactly the mechanism that produces
crest-factor drift and warble. This explains the §5c ablation showing harmonic-lock-off as
bit-identical on every real file, and it explains why the saw was fixed (shape 0.79 → 0.999)
while the complaint about real material survived. The gate is the cause: a stack needs a root
≥50 Hz within 40 dB of the loudest, k=2 or 3 present, ≥3 members, ≥30% of frame energy, *and*
a root that has passed the steady-tone frequency lock. Real sources rarely satisfy all six.

### 5e.3 The drift is measurable and up-asymmetric

`waveform_shape_consistency` (the metric built for the saw wobble) on shifted renders, where
the source is periodic enough for it to mean something:

| file | input | unity | **up5** | down5 |
|---|---|---|---|---|
| 300hzSaw | 0.999 | 0.999 | **0.984** | 0.993 |
| 440sawtooth | 0.996 | 0.996 | **0.973** | 0.987 |
| Female sung line | 0.336 | 0.312 | **0.093** | 0.244 |

Up-shift degrades waveform shape more than down-shift on all three, and on the Female line it
collapses to a quarter of the input's value. That is the crest-drift signature, it is
up-asymmetric as reported, and it sits on files where the lock covers 2.7–8.8% of partials.

(Piano 0.745–0.751 and the mixes 0.04–0.16 are outside this metric's competence — it needs a
single periodic f0 — so they neither support nor contradict.)

## 6d. PLAN — generalize shape invariance, then handle percussion separately (2026-09-25)

Grounded in the supplied reference's P1–P4 and the measurements above. Ordered by evidence,
not by expected size.

### D1 (build first) — relative-phase-shift (RPS) synthesis, replacing the gated lock

Current: capture θ_k once (`L.captured`) and freeze it, only for stacks passing six tests.
RPS (Saratxaga 2009): **θ_k(t) = φ_k − k·φ_1, re-measured every analysis frame and
interpolated**, with synthesis φ_k' = k·φ_1' + θ_k where only φ_1' integrates β·ω_1.

Why this is the right first move: it removes the steady-root requirement, which is the main
reason the lock never fires — we no longer need the root's frequency to be stable, because the
relationship is re-measured rather than assumed. It also subsumes the frozen-offset design,
and per the reference it should let us *remove* most of the frequency-EMA hacks, since jitter
in upper partials stops accumulating as relative-phase error when only φ_1 is integrated.

**Prospective gate — declared before building, per §5d.2's lesson:**
1. lock coverage on Female/Piano/Happy rises from 2.7/5.4/14.6% to >50%;
2. `waveform_shape_consistency` at up5 recovers toward the input (Female >0.25, saws >0.99);
3. `env_p2p_full` and `jitter_db` at ±5 **no worse** than now — this is the pair that
   correctly predicted the last listening result, and RPS must not trade warble for shape;
4. then blind A/B. If (3) fails, stop — that failure mode is exactly what testers rejected.

### D2 — a transient/percussive layer, because RPS cannot cover drums

Shape invariance is only defined for quasi-harmonic sources. A drum loop has no k·φ_1
relationship, so D1 will leave its 97.5% untouched. Percussion needs the other route
(Verma & Meng; Levine & Smith; Fierro & Välimäki 2023): detect transients, carry that layer
through with envelope filtering rather than resynthesising it from sinusoids, and reset
partial phases at onsets instead of inheriting rebirth continuity across them.

This is the one that addresses the drum loop specifically, and the drum loop is what the
listener singled out. It is larger than D1 and depends on nothing in it, so it can run in
parallel.

### D3 — spectral-envelope amplitude resampling, with a continuous knob

`A_k' = E(β f_k)` instead of `A_k`, with `A_k' = A_k^{1−γ}·E(β f_k)^γ` so γ is a dial.

Two honest caveats. **It was prototyped and rejected by ear in July** (see
[[project_neural_direction]]) — on an engine 10 dB worse with no joint solve, so worth
retrying, but not a fresh idea here. And **the symptom does not match**: the reference calls
this the most likely cause of up-shift sounding worse, but its artifact is *chipmunk timbre*,
whereas the reports are "slight gaps" and "amplitude modulation". Envelope preservation may
well improve naturalness without touching what was complained about. Test it, but as a
timbre improvement with its own A/B, not as the fix for this defect.

### D4 — a stretch-then-resample reference mode inside our own tooling

The reference's recommendation 4, and cheap: WSOLA (or our own oscillator bank time-stretched
by β, then resampled by 1/β) purely as a comparison render. It gives us the Audacity-class
baseline on demand without external tools, which §6c/C2 currently depends on Riley for. Not a
product path — it inherits the formant shift — but a permanent A/B anchor.

### Explicitly not doing

* Adopting WSOLA as the engine (§5e.1 — it works by having no model, which forfeits
  everything parametric this project exists for).
* Another perceptual metric before D1's prospective gate is tested (§5d.2).
* `JOINT_SHIFT`, noise placement, Nyquist — closed.

## 5f. D1 (RPS synthesis) FAILS ITS PRE-DECLARED GATE (2026-09-25)

Built as specified in §6d/D1: theta_k = phi_k - k*phi_root re-measured from the analysis
phases every frame, wrapped-EMA smoothed, members re-anchored per frame to
phi_k' = k*phi_root' + theta_k, keeping each member's own frequency (no harmonic
quantisation), with a looser tolerance since RPS absorbs deviation rather than accumulating
it. `RPS`, `RPS_TOL`, `RPS_MAXK`; default 0 and shifted output byte-identical when off.

**Gate 1 — coverage — PASSES, dramatically.** Locked track-frames at up5:

| file | frozen | RPS (2% tol) |
|---|---|---|
| DrumLoop | 2.5% | **96.6%** |
| Piano | 5.4% | **94.8%** |
| Female | 2.7% | **95.6%** |
| HappyMono | 14.6% | **99.1%** |

So the diagnosis in §5e.2 was right: the tolerance, not the steady-root test, was what kept
the lock from firing, and RPS does let it fire.

**Gates 2 and 3 — FAIL.** At +5 semitones, across tolerances 2% / 0.5% / 0.2% and with k
capped at 3 / 6 / uncapped:

| file | metric | input | frozen | RPS (best of six settings) |
|---|---|---|---|---|
| Female | shape_consistency | 0.336 | 0.093 | 0.167 (target >0.25) |
| Female | jitter_db | 2.135 | **1.323** | 1.436 — **worse at every setting** |
| Piano | shape_consistency | 0.751 | **0.745** | 0.709 — worse at every setting |
| 440saw | shape_consistency | 0.996 | 0.973 | 0.974 (target >0.99) |
| HappyMono | jitter_db | 1.361 | **1.568** | 1.679 — worse at most settings |

Gate 3 said "`env_p2p_full` and `jitter_db` at ±5 no worse", and "if (3) fails, stop — that
failure mode is exactly what testers rejected". It fails on the Female at all six settings
and on HappyMono at most. **Stopped, per the rule.** `harmonic_rps_mode` stays 0.

### 5f.1 Why it probably failed — and it is not the obvious reason

My first hypothesis was k-amplification: RPS anchors a member at k*phi_root, so error in the
root's *propagated* phase is multiplied by k. `RPS_MAXK` tests it directly and **capping k
does not help** (Female shape 0.030 at k<=3 vs 0.038 uncapped; jitter 1.449 vs 1.436). So
that is not the whole story.

The likelier explanation is the opposite of the intended one. Under shift the root's phase is
*propagated*, so it already carries accumulated error. Independent per-track integration
gives every partial its own **uncorrelated** error, and uncorrelated errors partially cancel
in the sum. Anchoring every member to a shared root makes those errors **correlated**, so
they add coherently and appear as waveform-level modulation — which is exactly what
`jitter_db` and `env_p2p` report, and exactly what a listener called "tremolo" and "a pitch
that can't hold itself" in §5d. Shape invariance assumes an accurate common reference; with a
drifting one, sharing it concentrates error instead of removing it.

This is consistent with the one place the frozen lock does work: the 300 Hz saw, where the
root is a stable synthetic tone whose frequency is estimated almost exactly (69% coverage,
shape 0.79 -> 0.999). Shape invariance is not limited by harmonicity here — it is limited by
**reference accuracy**, and §4d.10 already measured that real material's frequency estimates
plateau at ~3 Hz error.

### 5f.2 What this closes and what it leaves

* **Closed:** RPS/shape-invariant synthesis as a fix for shifted real material, unless and
  until the root frequency estimate improves — which §4d.10 established is an analysis-window
  problem, i.e. architectural. Kept behind `RPS=1` as a documented negative with working code,
  because it becomes viable the moment reference accuracy does.
* **Unaffected:** D2, the transient/percussive layer. It never depended on D1, it is the one
  aimed at the drum loop the listener singled out, and it does not need a harmonic reference
  at all — which, after this result, is a point in its favour rather than a limitation.
* **Re-weighted:** D3 (envelope amplitude resampling) now looks relatively more attractive,
  since it is an amplitude-domain change that does not touch phase and so cannot produce this
  failure mode. Its caveats from §6d stand (already rejected by ear once; chipmunk timbre is
  not the reported symptom).

**Method note.** Declaring gate 3 before building is what stopped this becoming a second
rejected listening test. The cost was one day; the previous undeclared version cost a
listening round and a wrong conclusion in the doc.

## 6. THE PLAN AFTER SEPTEMBER (2026-09-23)

The four changes landed since the re-baseline — file-start credit, transient short-frame
amplitude, the residual envelope cap, and the amplitude-EMA question — are **all inaudible
to the listener on this corpus**, while being clearly correct on measurement (+16 to +24 dB
of pre-echo, +2 dB of SRR). That is the single most important fact for planning: *residual
SRR is no longer the binding constraint on perceived quality*, and we currently have no
instrument that measures what is.

Everything below follows from that.

### Phase A — build the instrument that can still hear a difference (do first)

Two halves, both cheap, and neither is engine work.

**A1. Widen the corpus.** Every conclusion in this document rests on 15 files that the
engine has been tuned against for three months. New material is where an audible defect
will actually be found. Target ~20 more: speech (male and female), solo strings and winds,
harpsichord/glockenspiel, acoustic guitar, dense electronic, a full mix with vocals, a
sparse ambient pad, and something deliberately ugly (distorted guitar, heavy compression).
Run the existing battery over it and look for classes that score far below their peers.

**A2. A listening protocol that can resolve small differences.** Informal A/B has now
failed four times in a row to separate renders that differ by 16–24 dB on a targeted
metric — which is the expected outcome for masked artifacts, not a failure of the listener.
A short ABX harness (same clip, randomised A/B/X, forced choice, ~10 trials) over the
contested pairs would answer "is this audible at all" definitively, and would let every
future change be gated on a number that means something perceptually.

**Gate for Phase A:** a list of classes where the engine is audibly (ABX-confirmed) worse
than the input, ranked. If that list is empty on a 35-file corpus, the fidelity project is
finished and the work is Phase C.

### Phase B — fix what Phase A finds

Deliberately unspecified. The defects worth fixing are the ones a listener can identify, and
we do not yet know what they are. Two candidates are already on the books and should be
re-checked under the Phase A protocol before any more effort:
* **the cymbal chirp under shift** — the oldest open item, believed audible, never resolved;
  the residual envelope cap (§4d.12) did not touch it because it is the tonalised tracks'
  frequencies, not the noise envelope.
* **the amplitude EMA** (§4d.13) — +2 dB, inaudible in informal listening, one env var
  either way. Decide it with ABX rather than argument.

### Phase C — capability, which is what the parametric model is actually for

The engine's premise is that audio becomes *parameters*, not just that it reconstructs well.
Pitch shift is the only thing built on that premise so far. The decomposability constraint
has been paid for repeatedly; this is where it pays back.

* **Time-stretch.** Parameters are already time-indexed with explicit frequencies, amplitudes
  and phases; stretching is re-timing the synthesis grid, and it is the one operation a
  parametric model does better than a phase vocoder. Largest capability-per-effort on the
  list.
* **Formant-preserving shift.** Prototyped and rejected by ear in July, on an engine that was
  10 dB worse and had no joint solve. Worth re-testing now, behind a switch.
* **Component-level editing** — isolate/mute/gain a partial or a harmonic stack, which the
  harmonic-lock grouping already identifies. The residual meter makes it honest about what
  is not captured.

### Phase D — the analysis rebuild, only if Phase A demands it

§§4d.10–4d.12 established that the remaining per-class gaps (choir 13 dB, Fairlight 11,
Piano/SaintSaëns 7–8) are all one thing: a single-window STFT with constant-amplitude,
constant-frequency atoms cannot represent material that modulates within any window long
enough to resolve it. Closing that needs multi-resolution analysis with modulation-aware
atoms, and for choir specifically, multi-f0 separation so each voice can be warped the way
pitch-sync warps one.

This is a rebuild of the analysis stage, not a fix to it. **Do it only if Phase A shows
fidelity still limits what a listener hears** — otherwise it is a large investment in a
number nobody can perceive.

### Not worth doing

* More estimator work at the current window (§4d.10, §4d.11 — measured, three times).
* Per-band adaptive window as a standalone project: worth 2–3 dB (§4d.11), which on current
  evidence is inaudible. Fold it into Phase D if Phase D happens.
* Anything learned in the signal path (§3, unchanged).
* Performance work (§4d.14).

### Step 8 — per-band adaptive analysis window

Worth a measured 2–3 dB (§4d.11) and it is the only broad lever left. Choose the analysis
window per frequency band from a local stationarity test — the machinery already exists in
miniature as the LF tier's `lf_max_flux` gate, and `LF_CUTOFF`/`LF_MAX_FLUX` are now
env-exposed for experiments. Expect it to be inaudible on this corpus, like the transient
amplitude fix; justify it as headroom, not as a listening win.

**Gate:** no class regressed; choir 200–800 Hz band +4 dB; whole-file +1.5 dB or better on
choir/SaintSaëns; ears neutral.

### Beyond Step 8 — what is left is architectural

§4d.11 measured the remaining per-class gaps as analysis time-frequency limits. Closing them
means atoms that model amplitude and frequency modulation, and a way to estimate their
parameters under inter-partial interference — a different front end, not a fix to this one.
Worth doing only if someone wants to rebuild the analysis stage; the current one is close to
its practical ceiling.

### Explicitly dropped or demoted

* **S+T+N transient component** — the residual at onsets is −20 dB (§4d.6a): the energy is
  already modelled, so a transient *component* is not what this defect needs. Revisit only
  if Steps 1–3 fail.
* **Noise component for cymbals** — still the one structural case, but the gap is now ~2 dB
  at matched budget, behind everything above.
* **Anything learned in the signal path** — §3, unchanged and strengthened.
* **Unshifted paste of the original or the residual at onsets** — §4d.6a.

## 5. The plan (superseded — kept for the record)

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
