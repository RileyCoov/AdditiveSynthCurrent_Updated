# Root-Cause Diagnosis of the Additive Pitch-Shift Engine (commit e40c681): a Transient-Handoff Gain Hole and Phase-Incoherent Duplicate Partials

Yes, your engine has a real DSP error. In fact it has two separate engine defects, and neither one is a "birdie". **Defect T** is a transient-segment handoff locked to the frame grid. Under shift it cuts a 16.7–19.3 ms broadband hole (−5 to −6.5 dB) and puts an HF pre-echo of about +11 dB in front of every detected hit [M]. **Defect P** is a set of spurious partials that are the same in both shift directions and sit +14–19 dB above the source. They are already present in the analysis/tracking/least-squares parameter set and only become audible once phases are propagated instead of re-measured [M][I]. 1985 does not have a different bug. It exposes Defect T about 3.4–3.8× more often than the other files, and its Defect-P tones fall in the sustained midrange at the loudest level in the frame [M].

**TL;DR**
- **What is broken.** (1) Shifting renders each transient segment with long and short layers whose weights do not add to 1, and with a pre-onset long frame whose analysis window already contains the attack. The result is a −5…−6.5 dB hole over ≈S+1480…S+2460, +11.3/+11.4 dB of 6–16 kHz pre-echo, and attacks that arrive 10–13 ms late [M][I]. (2) Partials that nearly coincide (coasting tracks, anti-phase joint-LS pairs, tier duplicates) cancel when phases are re-measured. Under independent phase propagation they add incoherently, giving smooth +14–19 dB tonal swells. 17 of the top 25 up-shift events have a matching down-shift event (p≈0.002) [M][I].
- **Why 1985.** The cause is exposure, not a unique failure. 1985 has the densest hits (3.70 segments/s) with steady bass, guitar and vocal underneath, which gives 1.96–2.17 audible holes/s versus 0.57/s on HappyMono. Its spurious-partial rate (5.4/s prominent, up-shift) is no higher than HappyMono's (6.2/s) [M].
- **What you are mishearing.** The tones are real but come from two different populations. The "periodic" part (2.88–3.00, 4.08–4.29, 6.31–6.88 s) is mostly the per-hit dropout plus pre-echo, which is not tonal. You missed the two strongest tones (1.50–1.58 s at 1,205→1,805/809 Hz, and 1.60–1.67 s at 925→1,383/621 Hz). Unity cannot validate the representation, because 26.6% of it is bit-exact passthrough and the rest has its phases re-measured every frame [M][D].

**Provenance legend:** [M] = measured on the 60 WAVs this session. [D] = your Oct 1 2026 engine description or notes (main.cpp was **not** in the upload, so every code-level statement is a hypothesis to check). [L] = literature. [I] = my inference or derivation.

---

## 1. Executive diagnosis

**There is a bug, and it is a combination (class G):**

| Defect | Class | Status | Core evidence |
|---|---|---|---|
| T: transient-segment handoff (gain hole + pre-echo + late attack) | A reconstruction + F transient handling | **Confirmed engine artifact**; mechanism **very likely** | Locked to the frame grid (minimum at S+~2070, IQR ±130 samples, slope ≈0 against detected onset); identical up/down; −0.2 dB at unity because passthrough replaces S…S+4096 [M] |
| P: direction-symmetric spurious partials | B analysis + C tracking, revealed by E (phase propagation) | **Confirmed artifact**; near-coincident-atom mechanism **very likely**, not proven | +14–19 dB over the mapped source, smooth 100–200 ms bell, absent at unity (within ~1 dB), 47/79 prominent events twinned [M]; your own result that forcing propagated phase at ratio 1.0 closes 79% of the gap on 1985 [D] |
| Bulk HF "musical noise" field | B/D (noise rendered as sinusoids) | **Confirmed**; partly a normal consequence of the design | 127/s up, 104/s down on 1985, 87–94% above 4 kHz, median −30/−33 dB re frame max [M] |
| Output ceiling → global AM | D/post | **Confirmed** (very likely a limiter) | Samples pinned at exactly 0.95050048828125 = 31146/32768; identical AM on every harmonic of the 440saw only [M] |
| 200 Hz LF-tier step | B/D | **Very likely** an engine artifact | HappyMono 180–200 Hz +6.0 dB vs 200–220 Hz +0.1 dB [M] |
| HF long-term excess | D/residual | **Confirmed** measurement; cause **possible** | 4–15 kHz +1.9…+2.7 dB on dense mixes, both directions [M] |

**Is unity correct?** Unity is correct as a *render* and not as a *representation*. Measured where the model is actually used, 1985 unity SER is 14.6 dB, i.e. error power 10^(−1.46) ≈ 3.5% of signal power [M][I]. It adds 0 tonal events on 1985 and has no hop-grid ripple (≤0.34 dB pk-pk) [M]. But unity is built to hide both defects. It re-derives phase from analysis every frame, and it crossfades the original input into every transient segment [D]. So unity sounding fine says nothing about whether the parameter set is sound [I].

---

## 2. Original vs. Unity

| Check | Finding | Tag |
|---|---|---|
| Passthrough fraction | 26.6% of active samples on 1985 equal the input within 1 LSB, in 41 runs of ~2,750–3,100 samples, one per drum hit (≈0.256 s spacing). DrumLoopShort 25.5%, HappyMono 17.9%, take-me-out 0.2% | [M] |
| Segment geometry | Each passthrough window is exactly 4 hops (4096 samples), S=k·1024 → S+4096 (34/35 on 1985; 37/38 HappyMono; 8/9 DrumLoopShort). Model→original crossfade completes ~600 samples after S; original→model ramp takes the last ~800 samples | [M] |
| SER (model-only regions) | 1985 14.6 dB (16.4 overall), flat across bands 11–17.5 dB; DrumLoopShort 9.6; HappyMono 15.4; take-me-out 16.4; steady tones 28–43 dB | [M] |
| Added tonal components | 0 events on 1985, 0–3 on other dense files | [M] |
| Hop-grid AM/ripple | ≤0.34 dB pk-pk folded mod 1024 (dense); ≤0.01 dB steady tones | [M] |
| Phase | 300 Hz sine: exact frequency, PM sidebands at 300±46.875 Hz at −69 dB. 46.875 Hz = 48000/1024 = frame rate, i.e. the phase re-anchoring at each frame leaves a fingerprint; inaudible | [M][I] |
| Clicks/discontinuities | None added at any ratio | [M] |
| Harmonic amplitude at Defect-P sites | Unity tracks the original within ~1 dB | [M] |
| Transients | Not assessable on 1985/HappyMono/DrumLoopShort: S…S+4096 *is* the input. take-me-out has no passthrough segments, so its transient segments could not be located | [M] |
| LF loss, noise-floor modulation at unity | No band-specific deficit visible in the flat band SER. Noise-floor modulation was not separately measured at unity | [M] |

**Why phase re-measurement masks the parameter-set defect [I].** Suppose a frame contains two atoms at f₁ ≈ f₂ (|f₁−f₂| < 1 bin = 11.72 Hz). At unity each is rendered with its *measured* phase, so the rendered complex amplitude is

  A_unity = a₁e^{jφ₁} + a₂e^{jφ₂},  |A_unity|² = a₁² + a₂² + 2a₁a₂cos Δφ,

which the joint LS solve has fitted to the observed |X| by construction. Any decomposition whose complex sum matches the frame reproduces the frame, including a near-cancelling anti-phase pair (cos Δφ ≈ −1). This is why your own [D] result, that forcing the propagated-phase path at ratio 1.0 closes 79%/45%/24% of the unity→shift spike gap on 1985/HappyMono/DrumLoopShort, is the decisive observation. Unity validates the *sum*, not the *atoms*.

**Answer to Q2:** at ratio 1.0 the engine has a **latent** reconstruction problem: wrong atoms whose errors cancel, plus a transient handoff that is never exercised. It produces no **audible** reconstruction error at unity, because the masking mechanisms are working as designed [M][D][I]. The −69 dB frame-rate PM is real but irrelevant.

---

## 3. Up a Fifth vs. Down a Fifth

| Metric (1985 unless noted) | Up | Down | Unity | Tag |
|---|---|---|---|---|
| Dropout, bass 40–160 Hz, [S+1536,S+2304] self-referenced | −5.7 dB | −5.8 dB | −0.2 dB | [M] |
| Dropout, mids 500–2k | −4.4 dB | −5.0 dB | — | [M] |
| HF pre-echo 6–16 kHz, [S+512,S+1536] vs original | +11.3 dB | +11.4 dB | — | [M] |
| Between-hit HF baseline excess | +2.0 dB | +2.5 dB | — | [M] |
| HF onset lateness (median, IQR ±15 ms) | +13.0 ms | +10.5 ms | — | [M] |
| Sharp HF impulses (original 11.1/s) | 0.8/s (shifted) | 0.8/s (shifted) | — | [M] |
| Prominent spurious partials | 50 (5.4/s, median src 3,044 Hz) | 29 (3.2/s, 2,564 Hz) | 0 | [M] |
| Bulk spurious events | 1,166 (127/s) | 961 (104/s) | — | [M] |
| 440saw AM (all harmonics alike) | 0.8–1.0 dB pk-pk, ~3.3 Hz | 1.5–2.4 dB pk-pk, ~21–30 Hz | — | [M] |
| Pure-sine accuracy | 449.571 vs 449.490 Hz (+0.3 cent) | 200.261 vs 200.227 Hz (+0.3 cent) | exact | [M] |
| 4–15 kHz long-term level | +1.9…+2.7 dB | +1.9…+2.7 dB | — | [M] |

**Interpretation.**
- **Defect T does not depend on ratio** (±0.1 dB pre-echo, ±0.1 dB bass dip up vs down) [M]. So it lives in segment scheduling and OLA weighting, not in frequency mapping, oscillator phase increments or resampling [I].
- **Defect P lives in source coordinates.** 17 of the top-25 up events have a down twin at the same *source* frequency (±3%) and overlapping time, against a circular-shift null of mean 4.5 and 99th percentile 9 (p≈0.002) [M]. An error created by transposition (mapping, aliasing, phase-increment scaling) would not reappear at the same source frequency for ratios 1.5 and 0.667 [I]. The parameter set shared by both renders is the origin.
- **Asymmetries that are real.** Up-shift has 1.7–2× more prominent events on all three dense mixes [M]. Candidate reasons [I]: (a) the residual is transposed on down-shift only [D], so the down output carries a transposed noise bed that may partially mask stubs; (b) mapping 2.5–3 kHz content up to 3.8–4.6 kHz moves it out of the dense vocal/guitar region where masking is strongest. The shipped up-shift-only duplicate suppression [D] should *reduce* up events, so it cannot explain the excess. Status: **unable to determine**; see E7. The AM rate difference (3.3 Hz vs 21–30 Hz) belongs to the ceiling/limiter (§6.6), not to the spectral engine [M][I].
- **Mechanisms on your list that the data rule out or put out of scope:** phase wrapping and oscillator resets (no clicks; the only phase artifact is −69 dB PM) [M]; hop-rate OLA errors on steady oscillators (≤ −100 dB hop-rate AM, ≤0.08 dB ripple under shift) [M]; bin quantisation and pitch-scaling error (+0.3 cent) [M]. Aliasing was **not measured**, so no claim is made.

---

## 4. Artifact inventory table

Samples are at 48 kHz. For the per-segment rows, S is the segment start listed in M1, snapped to k·1024.

| Time (s) | Samples | Frequency | Duration | Output | Artifact type | Confidence | Likely cause |
|---|---|---|---|---|---|---|---|
| 1.24–1.30 (±40 ms audible) | 59.6k–62.2k | src ≈530 Hz → 797 up / 352 down | ~100–200 ms bell | Up + Down | Spurious partial, +15–16 dB, within ~1 dB of frame max | Confirmed artifact [M] | Defect P: near-coincident atoms incoherent under propagation [I] |
| 1.50–1.58 | 72.2k–75.8k | ≈1,205 → 1,805 / 809 Hz | ~100–200 ms | Up + Down | Spurious partial, +14–18.5 dB, *loudest in frame* | Confirmed [M] | Defect P [I] |
| 1.60–1.67 | 76.8k–80.1k | ≈925 → 1,383 / 621 Hz | ~100–200 ms | Up + Down | Spurious partial, +18–19 dB | Confirmed [M] | Defect P [I] |
| 2.04–2.19 | 97.8k–105.2k | ≈2,777 → 4,160 / 1,852 Hz; stubs near apex and ≈2,580 Hz | ~150 ms | Up + Down | Glide replaced by steady stubs, +16–20 dB | Confirmed [M] | 25 Hz gate vs ≈75 Hz/hop glide → coasting at constant frequency [D][I] |
| 2.83–2.89 | 135.8k–138.7k | 4,559 Hz (src 3,042) | ~60 ms | Up only | Spurious partial, −6 dB re frame max | Confirmed [M] | Defect P or up-only path; undetermined |
| 4.37–4.69 | 209.7k–225.3k | ≈3,441 → 5,156 / 2,297 Hz | recurring ~300 ms | Up + Down | Spurious partial, +6–12 dB, up to frame max | Confirmed [M] | Defect P (+ coasting persistence) [I] |
| 5.45–5.68 | 261.6k–272.6k | ≈3,500 → 5,121–5,250 / 2,332 Hz | ~230 ms | Up + Down | Spurious partial | Confirmed [M] | Defect P [I] |
| Every segment S (35 on 1985; e.g. S=141,312 / 2.944 s) | S+512…S+1536 (e.g. 141.8k–142.8k) | 6–16 kHz (+mids 500–2k +1.5 dB) | ≈21 ms | Up + Down | Pre-echo, +11.3/+11.4 dB | Confirmed [M] | Last long frame analyses past the onset [I] |
| Every segment S (e.g. 2.944 s → 2.975–2.995 s) | S+1480…S+2460 (e.g. 142.8k–143.8k) | Broadband: 40–160, 500–2k, 2–6 kHz | 16.7–19.3 ms | Up + Down | Gain hole −5.1…−6.5 dB, minimum at S+~2070 | Confirmed artifact; mechanism very likely | Non-complementary long/short weights [I] |
| Every segment S | after S+~2400; tails ≥80 ms | HF | 10–13 ms late onset; +4–5 dB tails | Up + Down | Late re-attack, transient smearing; HF impulses 11.1/s → 0.8/s | Confirmed [M] | Short layer starts at grid onset, not true onset [I] |
| Whole file, 1985 | — | 87–94% at src >4 kHz | short events | Up 127/s, Down 104/s | Bulk musical-noise field, median −30/−33 dB | Confirmed; partly normal | HF noise rendered as sinusoids [I] |
| 5.81, 7.39 (and 6.31–6.88 HF events) | — | src 7–8.5 kHz | short | Up/Down | Weak HF events, −9 to −34 dB re frame max | Confirmed (weak) [M] | Bulk musical-noise class |
| HappyMono / take-me-out, whole file | — | 80–200 Hz src band; step at 200 Hz | stationary | Shifted | +2…+6 dB LF excess (180–200 Hz +6.0, 200–220 +0.1) | Very likely engine | LF-tier/long-tier duplication [I] |
| 440sawtooth up/down, whole file | 107/107, 101/107 25-ms blocks at ceiling | all harmonics | continuous | Up/Down | Global AM 0.8–1.0 / 1.5–2.4 dB pk-pk; −0.8…−1.0 dB level loss | Confirmed | Output ceiling/limiter [M][I] |
| 440sawtooth source | 132,300 samples | f0 = 478.9 Hz | — | Source | Mislabelled 44.1 kHz file played at 48 kHz | Confirmed source problem [M] | Wrong header/import |

---

## 5. Why 1985 is different

The difference is exposure, not complexity. These are the signal characteristics that matter [M]:

1. **Hit density.** 3.70 transient segments/s, against 2.63 for HappyMono and 3.10 for DrumLoopShort, at 8th-note spacing ≈0.256 s. Every one of these is a Defect-T event under shift.
2. **Steady sustained material under the hits.** The hole is only audible as a *dropout* if the band was steady before the hit. 1985 has 18/34 (bass) and 20/34 (mids) segments with the original band level in the dropout window within ±2 dB of the pre-hit level. That gives **1.96 and 2.17 exposed holes/s**, against 0.57/s on HappyMono and **0/8** on DrumLoopShort, where the hole lands on a rising attack and sounds like a weak or "hollow" attack instead. The dip is the same depth everywhere (−5.0 to −6.7 dB) [M]. The ratio is therefore ≈3.4–3.8× more exposure, not a different mechanism.
3. **Defect-P tones fall in the exposed midrange.** The prominent-event *rate* is not unique: HappyMono 6.2/s and take-me-out 5.6/s against 1985 5.4/s (up) [M]. But 1985's strongest events sit at 530, 925 and 1,205 Hz source frequency, at or within ~1 dB of the frame's loudest component [M]. That band is unmasked and pitch-salient [I].
4. **Fast vocal glides.** At 2.04–2.19 s the source glides ≈2,400→2,750→2,400 Hz in ≈0.1 s, ≈75 Hz per 1024-hop, three times the 25 Hz gate [M]. This produces a tracking failure that is not a transient.
5. **A probable reason the LF-tier defect does *not* show on 1985 [I].** The 16384-point LF FFT spans 0.341 s, longer than the 0.256 s hit spacing, so every LF window contains a hit. If the "bass is stationary" test is energy-based, it likely fails on 1985 and the LF tier stays off. This would explain why the 200 Hz step appears on HappyMono and take-me-out but not 1985. E4 tests this indirectly.

**README correction.** "Bass dropping out ~19 ms at each drum hit" is the right duration but the wrong description. The hole is broadband (bass −5.7/−5.8, mids −4.4/−5.0 dB), it is not specific to down-shift, and two of the four README times (1.350, 4.174 s) fall just *before* segment starts, not inside them [M].

---

## 6. Code-level diagnosis

**Caveat:** main.cpp was not uploaded. Each item gives the code region implied by [D], the intended math, my hypothesis about the actual behaviour [I], why it would produce the measured artifact [M], which renders it affects, why it is louder on 1985, and a test.

### 6.1 Transient handoff: long-frame suppression without a complementary window (Defect T, gain hole)
- **Where:** the transient-segment scheduler that suppresses long frames near an onset, and the OLA normaliser ("OLA normalised by window sum") [D].
- **Intended:** y(n) = Σ_m w_L(n−mH)·s_m(n) + Σ_k w_S(n−kH_S)·t_k(n), with **W(n) ≡ Σ_m w_L + Σ_k w_S = 1 for all n**. The two layers render the *same* waveform coherently, so the right condition is amplitude-complementary. The power-complementary condition w²(n)+w²(n+N)=1 is the Princen–Bradley/TDAC requirement for windows applied at both analysis and synthesis [L: Edler 1989 as summarised in US 8,095,359].\[1\]\[2\] It applies only if the layers are incoherent. For your symmetric Hann h[n]=0.5−0.5cos(2πn/(N−1)), N=2048, H=1024: h[n]+h[n+1024] = 1 − 0.5[cos(2πn/2047) − cos(2πn/2047 + π/2047)] ≈ 1 − O(7.7·10⁻⁴), a ripple of ≈0.007 dB. That fits the measured ≤0.01 dB on steady tones [I][M].
- **Hypothesised actual:** the long frame centred at S+2048 is dropped. The frame centred at S+1024 (span S…S+2048) and the frame at S+3072 (span S+2048…S+4096) remain. Then
  W_L(n) = ½ + ½cos(π(n−S−1024)/1024) on [S+1024, S+2048), and ½ − ½cos(π(n−S−2048)/1024) on [S+2048, S+3072),
  which falls through 0.5 at **S+1536**, reaches 0 at **S+2048**, and climbs back through 0.5 at S+2560 [I]. Measured: the dip starts at S+1480–1590, has its minimum at S+~2070, and recovers at S+2380–2460 [M]. Recovery comes earlier than S+2560 because the short layer starts at ≈S+2048–2300 [M]. The depth is −5…−6.5 dB rather than −∞ because the band-power windows smooth it and the residual and short layer partially fill the hole [I].
- **Normaliser clue [I]:** if the code divided by the *actual* per-sample Σw of the rendered frames, the S+1024 frame's tail would be amplified instead of producing a hole. The observed hole therefore implies division by a fixed/nominal Σ (≡1), or a separately normalised short layer.
- **Affects:** up and down equally. Unity is untouched because S…S+4096 is replaced by input [M].
- **Why 1985:** §5, items 1–2.
- **Test:** E2.

### 6.2 Pre-echo: the last pre-onset long frame analyses through the attack
- **Where:** the analysis for the long frame centred at S+1024. Its 4096 window spans S−1024…S+3072 [D][M].
- **Intended:** parameters for a frame used before the onset should come from pre-onset samples only. This is the purpose of the START window in block-switching coders. In AAC (2048/256) the LONG_START window is the long rising half over 0–1023, ones over 1024–1471, the short falling slope over 1472–1599, and **zeros over 1600–2047** [L: ISO/IEC 14496-3 framing as reproduced in Fraunhofer patents US 9,734,833 / US 8,457,975].\[3\]\[4\] Edler introduced these transition windows [L: Edler 1989].\[5\]
- **Actual:** the true attack (half-rise of 2–16 kHz at S+1770–1800) lies inside the analysis window [M]. A stationary-sinusoid fit spreads the attack's broadband energy across the whole synthesis frame S…S+2048. This is classic pre-echo: error that "spreads out in time over the complete block and therefore will also precede the time-domain transient" [L: Erne, AES Convention Paper 5489, 2001, §2]. Measured: +11.3/+11.4 dB at 6–16 kHz in [S+512,S+1536], against a +2.0/+2.5 dB baseline; DrumLoopShort +22.3/+20.0 dB; HappyMono +7.7/+8.3 dB [M].
- **Answer to Q12:** yes, the engine does interpret transient spectral content as stationary partials, but only in this one frame (≈21 ms of broadband HF), not as a tonal birdie. Prominent tones are *under*-represented in segments: 19–22% of events against 33% time coverage [M].

### 6.3 Onset quantised to the frame grid → late short layer
- **Where:** the transient detector and the short-frame start index [D].
- **Actual:** the dip's position is locked to the grid (IQR ±130 samples, slope ≈0 against detected onset) [M]. The short layer appears at ≈S+2048–2300 while the attack is at S+1770–1910 [M]. That produces HF onsets +13.0/+10.5 ms late, late re-attacks after S+~2400, and sharp-impulse density falling from 11.1/s to 0.8/s [M]. This is transient smearing in the phase-vocoder sense [L: Röbel 2003; Duxbury, Davies & Sandler 2002].\[6\]\[7\]

### 6.4 Near-coincident atoms → incoherent power excess under propagation (Defect P)
- **Where:** (i) coasting tracks that "re-read magnitude/phase at the old bin" and die only 70 dB below their own maximum (with a 3-frame cap below 200 Hz only); (ii) the joint Hann-weighted LS solve (±8 bins, Tikhonov 1e−3); (iii) cross-tier duplicates [D].
- **Intended math (coherent):** the frame contribution is a₁e^{jφ₁}+a₂e^{jφ₂} with |·|² = a₁²+a₂²+2a₁a₂cos Δφ.
- **Under shift [I]:** each track propagates its own phase, φᵢ(t)=φᵢ(t₀)+2π r fᵢ (t−t₀). Δφ then drifts at 2πr(f₁−f₂) and is no longer the measured value, so the expected power becomes a₁²+a₂². The excess over a measured anti-phase pair is
  G_avg = (a₁²+a₂²)/(a₁−a₂)², G_peak = ((a₁+a₂)/(a₁−a₂))².
  With a₂/a₁ = 0.78: G_avg = 1.61/0.0484 ≈ 33 (**+15 dB**), G_peak ≈ 65 (**+18 dB**). That matches the measured +14–18.5 dB at 1.50–1.58 s and +15–16 dB at 1.24–1.30 s [I][M]. Your toy result (+7.7 dB level error at 0.3-bin spacing) has the same sign and mechanism [D].
- **Why the LS solve produces such pairs [I/L]:** for a Hann-weighted inner product, two unit atoms spaced δ bins have normalised correlation ρ(δ) ≈ sinc(δ)/(1−δ²). At δ=0.3, ρ ≈ 0.94, so the Gram condition number is (1+ρ)/(1−ρ) ≈ 34. A Tikhonov λ=10⁻³ is far below 1−ρ ≈ 0.06 and does almost nothing. For n closely spaced signals, the CR bound on frequency estimates "is proportional to (δω)^−2(n−1)/N·SNR. Therefore, the bound increases rapidly as the signal separation is reduced" [L: H.B. Lee, IEEE Trans. Signal Processing 40(6):1507–1517, 1992], and the bound "is a strong function of the phase difference" [L: Dilaveroglu, IEEE TSP 46(2):463–474, 1998 and Signal Processing 1998]. Sub-bin pairs therefore have poorly determined amplitude splits, and the anti-phase split is one the data cannot rule out.
- **Fit to the temporal shape [I]:** a smooth 100–200 ms bell is 4.7–9.4 hops of 21.3 ms. That is consistent with atom amplitudes interpolated across a few frames while a coasting track overlaps a moving partial. It is **not** consistent with birth/death discontinuities, which would give clicks or onsets with sharp edges, none of which were measured [M].
- **Answer to Q10 ("false oscillator"):** the onset and offset are not sudden, and there is no abrupt f/amp/phase jump [M]. The offending component is not an invented frequency. It is a weak real source component (the source has "only a weak or brief component" there) amplified by +14–19 dB through incoherent summation [M][I]. Estimated frequencies are given in §4.
- **Answer to Q11 (phase):** a phase error on a *single* correctly estimated partial changes its power by 0 dB. A wrap or reset gives a click or PM sidebands, which is broadband and short. A smooth narrowband +15 dB burst needs ≥2 atoms whose *relative* phase is wrong [I]. So "phase error" is correct only in the sense of the inter-atom relative phase, i.e. a loss of vertical coherence: as summarised for the phase vocoder, "the vertical phase coherence between bins is not maintained" [L: Laroche & Dolson, IEEE Trans. Speech and Audio Processing 7(3), 1999]. The oscillator phase arithmetic is not at fault.

### 6.5 Fixed 25 Hz gate + coasting on glides
- **Gate condition [I]:** for harmonic n with fundamental slope ḟ₀, the per-hop jump is |Δfₙ| = n·ḟ₀·H/fs. The track breaks when this exceeds 25 Hz, i.e. for n > 25·fs/(H·ḟ₀). At 2.04–2.19 s, Δf ≈ 75 Hz/hop > 25 [M]. The broken track coasts at a constant frequency while re-reading the old bin, and the new peak is born as a separate track. The outputs show steady stubs at the apex and near 2,580 Hz in both directions [M]. Illustrative only, with assumed values: ±50-cent vibrato at 5.5 Hz on f0=220 Hz gives a peak ḟ₀ ≈ 223 Hz/s ≈ 4.8 Hz/hop, so harmonics n ≥ 6 break the gate at the vibrato extremes [I].
- **Prior art:** linear-prediction tracking [L: Lagrange, Marchand & Rault 2007]; HMM tracking [L: Depalle, García & Rodet 1993]; McAulay–Quatieri birth/death matching [L: 1986]; frequency-slope estimation ψ̂₀=(ω₊−ω₋)·Fs from one-sample-shifted spectra [L: Marchand, DAFx-12].\[8\]\[9\]\[10\]\[11\]

### 6.6 Output ceiling / limiter
- **Evidence:** many shifted outputs are pinned at exactly 0.95050048828125 (31146/32768). On the 440saw the peak sits at that ceiling in 107/107 (up) and 101/107 (down) blocks, every harmonic shows the *same* AM (≈−36 to −38 dB, dominant 3.3 Hz up / 21–30 Hz down), and the level is uniformly −0.8…−1.0 dB low. The 300 Hz saw never reaches the ceiling and has 0.0 dB AM [M]. So the AM is global gain modulation, not a spectral-engine effect [I].
- **Why shift triggers it:** the shifted 300 Hz saw's peak wanders between 0.64 and 0.86 against an original 0.50 [M]. Its inter-harmonic phase relations drift, so the crest factor rises by up to 20log₁₀(0.86/0.50) ≈ +4.7 dB [I]. At a fixed power spectrum, the peak factor is set by the component phases [L: Schroeder 1970].\[12\]

### 6.7 200 Hz LF-tier boundary
- **Evidence:** HappyMono 180–200 Hz +6.0 dB vs 200–220 Hz +0.1 dB; 80–200 Hz +2…+6 dB on HappyMono and take-me-out [M].
- **Inference:** +6.02 dB = 20log₁₀2 is exactly what a coherent duplicate of an equal-amplitude atom gives; an incoherent duplicate gives +3.01 dB [I]. The 16384-pt tier's peaks below 200 Hz are probably *added to*, rather than *replacing*, long-tier peaks near the boundary, and frozen-offset harmonic lock may keep them in phase [I][D]. The shipped dedup is "full-band on up-shift only" [D], and cross-tier pairs may escape it.

### 6.8 Residual and HF partial admission
- **Evidence:** bulk events at 127/s (up) on 1985, 87–94% above 4 kHz, median −30/−33 dB [M]; HF level +1.9…+2.7 dB on dense mixes [M].
- **Inference:** the admission threshold (0.00025×file-max with +12 dB sensitivity ramped 1.5→8 kHz [D]) is designed to turn noise peaks into short-lived sinusoids. That is exactly the case that motivated SMS's stochastic component [L: Serra & Smith 1990].\[13\]\[14\]\[15\]\[16\] The residual mask max(0, |X_in| − |X_model|) [D] is half-wave-rectified magnitude spectral subtraction, the canonical musical-noise generator.\[17\]\[18\] Musical noise is "time varying tones" from isolated surviving peaks [L: Boll 1979; Berouti, Schwartz & Makhoul 1979; Vary 1985; Cappé 1994].\[19\]\[20\]\[21\]\[22\] The +2 dB HF excess is plausibly propagated-phase HF partials adding incoherently to a residual computed against coherent model magnitudes [I]. Status: **possible**.

---

## 7. Other errors you should know about

1. **Transient smearing and late attacks on every hit** under shift: HF onsets +13.0/+10.5 ms late, impulse density 11.1/s → 0.8/s, HF decay tails +4–5 dB for ≥80 ms [M]. You will hear this as a softened, "flammed" or slightly behind-the-beat snare/hat, not as a tone.
2. **Limiter AM** on loud, high-crest material (440saw: 0.8–1.0 dB up, 1.5–2.4 dB down, and a −0.8…−1.0 dB level loss) [M]. On dense mixes, any output pinned at 31146/32768 is suspect.
3. **LF boom** of +2…+6 dB in 80–200 Hz with a hard step at 200 Hz (HappyMono, take-me-out). Female_Sung_Line 40–200 Hz is +4…+9 dB, at low absolute level [M].
4. **HF brightening**, +1.9…+2.7 dB at 4–15 kHz on dense mixes [M].
5. **Loss of vertical phase coherence** ("phasiness") [L: Laroche & Dolson, IEEE Trans. Speech and Audio Processing 7(3), 1999]: the shifted saw's peak wanders 0.64–0.86 [M]. It is audible mainly as timbral roughness on bright stacked harmonics and as the driver of item 2 [I].
6. **Musical-noise field**, 104–127 events/s at −30 dB, mostly above 4 kHz [M]. It is heard as HF "sparkle/swirl", not as discrete tones.
7. **Inaudible but real:** unity frame-rate PM at −69 dB; +0.3 cent positive frequency bias [M].
8. **Source problem:** "440sawtooth" is 478.9 Hz = 440×48000/44100, i.e. a 44.1 kHz file relabelled 48 kHz [M]. Any pitch-accuracy test against 440 Hz on that file is invalid.

Not observed: added clicks, hop-grid ripple, or hop-rate AM on steady tones [M]. Not measured: aliasing, chorus/phasing indices, decay-time changes beyond the HF tails.

---

## 8. Most likely root cause

- **A. Is there a bug?** Yes. Defect T is a reconstruction/transient-handling bug, confirmed in the outputs, with its mechanism very likely. Defect P is a parameter-set defect revealed by phase propagation. It is confirmed in the outputs and its mechanism is very likely [M][I].
- **B. Where?** T: segment scheduler + OLA normalisation + last-long-frame analysis window + onset quantisation (§6.1–6.3). P: coasting re-read, joint-LS conditioning, cross-tier duplicates, fixed gate (§6.4–6.5) [I].
- **C. Why 1985?** Exposure: 3.70 segments/s × steady sustained bed → 1.96–2.17 holes/s, and Defect-P tones that land at 530–1,205 Hz at frame-max level [M].
- **D. Why less audible elsewhere?** HappyMono has 0.57 holes/s. DrumLoopShort's holes fall on rising attacks (0/8 exposed) and are heard as hollow attacks. Spurious tones on HappyMono are as frequent but sit in less exposed contexts [M][I].
- **E. Definitely present elsewhere?** Yes. The dropout appears on HappyMono (−4.8/−3.6 bass, −5.8/−5.4 mids), DrumLoopShort (−6.7/−4.8, −6.8/−5.4) and Fairlight C3 (≈−5 to −8 dB). Pre-echo appears on DrumLoopShort (+22.3/+20.0 dB) and HappyMono. Prominent spurious partials occur at 6.2/s on HappyMono and 5.6/s on take-me-out [M].
- **F. Which "musical noise" is something else?** The 2.88–3.00, 4.08–4.29 and 6.31–6.88 s regions are dominated by the segment dropout/pre-echo, which is not musical noise and not tonal. The prominent tones are phase-incoherent duplicate partials, not musical noise in the speech-enhancement sense. Only the bulk HF field is true musical noise [M][I][L].
- **G. Single most likely root cause of what you *called* the artifact:** near-coincident atoms whose complex sum is correct only under measured phases (Defect P). For the *periodic* component, the single cause is the non-complementary long/short handoff (Defect T). Both are rated very likely, neither proven, until E1–E3 are run.
- **H. Next experiments:** E1, E2, E3 first (§9).

### Hypothesis ranking

| # | Hypothesis | Evidence for | Evidence against | Confirming observation | Falsifying experiment |
|---|---|---|---|---|---|
| 1 | T: W_L+W_S<1 at handoff (dropped centre frame, nominal normaliser) | Dip spans S+1480…S+2460 with minimum S+~2070, matching the Hann model's 0.5 points at S+1536 and zero at S+2048; frame-locked; ratio-independent [M][I] | Depth only −5…−6.5 dB (explained by smoothing/fill, not proven) | E2 dump shows min W ≤0.3 near S+2048 | E2 shows W≡1±0.01 → look at short-layer amplitude estimation instead |
| 2 | T': pre-echo from last long frame analysing past onset | +11.3/+11.4 dB in [S+512,S+1536]; window S−1024…S+3072 contains the S+1770 attack [M] | — | Truncating that analysis window at n0−64 drops pre-echo to ≤+3 dB | Pre-echo unchanged after truncation |
| 3 | P: near-coincident atoms incoherent under propagation | Direction symmetry p≈0.002; absent at unity; +15/+18 dB fits a₂/a₁≈0.78 anti-phase; [D] 79% closure at ratio 1.0 | Not yet shown at atom level | E3 finds ≥2 atoms within <1 bin, complex sum ≈15 dB below a₁+a₂ | E3 finds a single atom at ≈ source level |
| 4 | P': fixed 25 Hz gate + coasting on fast glides | 2.04–2.19 s: 75 Hz/hop vs 25 Hz; steady stubs in both outputs [M] | Explains glide events, not the steady-region ones at 530/925/1,205 Hz | E6 removes the stubs | Stubs persist with a proportional/predictive gate |
| 5 | Transient energy turned into long-lived "birdies" | Pre-echo frame does render transients as partials | Prominent events under-represented in segments (19–22% vs 33%) [M] | — | Already largely falsified for the tones |
| 6 | Oscillator phase wrap/reset bug | — | No clicks; only −69 dB frame-rate PM; +0.3 cent accuracy [M] | — | Already falsified as a cause of tones |

---

## Are you hearing it correctly?

- **Are the tones real?** Yes. They are present in both shifted outputs, absent in the original and in unity, and +14–19 dB over the mapped source [M].
- **Are your locations accurate?** Partly.
  - **1.20–1.31 s: correct.** It is the 530 Hz event.
  - **4.46–4.64 s: correct.** It is the recurring 3,441 Hz event, overlapping segment 4.459.
  - **2.12 and 5.48 s: correct** (the glide stub and the 3,500 Hz event).
  - **2.88–3.00, 4.08–4.29 and 6.31–6.88 s: mostly mis-described.** These are the segment dropout/pre-echo at 2.944, 4.181, 6.251/6.507/6.763, with only weak HF events (−16 to −34 dB) [M]. You are hearing a rhythmic hole and HF "spit", not a tone.
  - **5.81 and 7.39 s:** weak HF musical-noise events, src 7–8.5 kHz [M].
- **What you missed:** the two strongest tones in the file, **1.50–1.58 s (1,205 Hz src, loudest component in its frame)** and **1.60–1.67 s (925 Hz, +18–19 dB)**, plus the 2.04–2.19 s glide stubs at +16–20 dB [M].
- **Wrong labels:** "spectral birdies" is the codec term from Erne (AES 111th Convention Paper 5489, 2001), where "A high frequency component of an instrument... may therefore momentarily disappear and become audible again, when the bit-allocation of a following block has changed. This artifact which is called 'birdies'" [L]. "Musical noise" is the speech-enhancement term for random isolated surviving peaks [L: Berouti 1979; Boll 1979].\[21\]\[22\] Your prominent tones are neither. They are deterministic, direction-symmetric, smooth swells from phase-incoherent duplicate partials. The best descriptive term is *spurious (phantom) partials from loss of inter-atom phase coherence* [I]. The periodic component is *pre-echo plus a block-switching gain hole*, i.e. a non-perfect-reconstruction transition [L: Edler 1989; Erne 2001].\[5\]\[23\]
- **Is it limited to 1985?** No. It is most *exposed* there. Every metric reproduces on HappyMono, DrumLoopShort and take-me-out [M]. Caveat: segments recur every 0.256 s and prominent events occur ~5/s, so a listed region coinciding with *some* event is weak evidence on its own [M].

---

## 9. Experiments to run next

Each experiment changes one variable. Predictions use the measured baselines.

| ID | Change (one variable) | Metric | Pass (supports hypothesis) | Fail (falsifies) |
|---|---|---|---|---|
| E1 | Shift path at ratio 1.0, passthrough **off** (propagated phase) | Segment-locked dropout in [S+1536,S+2304]; HF ratio in [S+512,S+1536]; level at the 3 Defect-P sites | Dropout −5…−6 dB; pre-echo ≈+10 dB; 1.24–1.30 s/530 Hz, 1.50–1.58 s/1,205 Hz, 4.37–4.69 s/3,441 Hz reappear at +12…+17 dB | Dropout ≥ −1 dB, or partials < +6 dB → defects originate in the transposition, not the shared representation |
| E2 | Dump per-sample w_L, w_S, W across 5 segments | W(n) | W<0.9 over ≈S+1500…S+2400, min ≤0.3 near S+2048 | W≡1±0.01 everywhere |
| E3 | Dump all atoms within ±1.5% of 1,205 Hz, 1.40–1.65 s (amp, phase, freq, coast flag, tier) | Cluster coherence | ≥2 atoms <1 bin apart, |Σaᵢe^{jφᵢ}| ≈15 dB below Σaᵢ, ≥1 coasting | Single atom at ≈ source level |
| E4 | LF tier off, HappyMono | 80–200 Hz excess; 180–200 vs 200–220 Hz step | Excess (+3…+6 dB) and step vanish or halve | Unchanged |
| E5 | Output ceiling off, 440saw up/down | AM pk-pk; band level error | AM <0.1 dB; level error ≈0 dB | AM persists ≥0.5 dB → spectral-engine AM |
| E6 | Gate Δf_max = max(25 Hz, 0.03·f), 1985 2.04–2.19 s | Stub level at 4,160/1,852 Hz | Falls from +16–20 dB to < +6 dB | Unchanged |
| E7 | Residual transposition forced on for up-shift (or off for down) | Prominent-event ratio up:down | Ratio moves from 1.7–2× toward 1× | Ratio unchanged → asymmetry is in the sinusoidal path |

**Metrics to adopt permanently:**
- **Dropout depth** D = median over segments of [(L_shift − L_orig) over [S+1536,S+2304]] − [(L_shift − L_orig) over [S−2048,S]], per band. Lock to S, not to detected onsets: onset-locked averaging blurs the result to −2…−3 dB [M].
- **Pre-echo ratio** PE = L_shift − L_orig at 6–16 kHz over [S+512,S+1536].
- **Direction-symmetry fraction** of prominent events, with a circular-shift null.
- **Exposed holes/s.**
- **Per-cluster coherence ratio** κ = Σaᵢ² / |Σaᵢe^{jφᵢ}|². κ > 1.26 (+1 dB) flags a Defect-P hazard before rendering.

---

## 10. Suggested fix

These are ordered by audible payoff per line of code. Each lists the number it should move and the risk.

**F1: Complementary start/stop windows with a sample-accurate onset (fixes T: hole, late attack).**
```cpp
// n0: onset from a sub-hop detector (e.g. 64-sample energy envelope, half-rise of 2–16 kHz)
// Ls = 256 short frame, short hop Ls/2 = 128
// Start frame = last long frame whose centre c* <= n0 - Ls/2
float wStart(int n){            // synthesis window of start frame
  if (n <  c*)            return hann2048(n - (c* - 1024));       // normal rising half
  if (n <  n0 - Ls/2)     return 1.f;                              // flat
  if (n <  n0 + Ls/2)     return 0.5f + 0.5f*cosf(M_PI*(n-(n0-Ls/2))/Ls); // short slope
  return 0.f;
}
// first short frame's rising half covers [n0-Ls/2, n0+Ls/2): wStart + ΣwS == 1 there
// stop window = mirror image at segment end; drop the separate fixed-sum normaliser
// unit test: max_n |ΣwL + ΣwS - 1| < 1e-3 over every segment
```
This mirrors AAC's LONG_START geometry (ones → short slope → zeros) [L].\[3\] Predicted: dropout D from −5…−6 dB to ≥ −1 dB; HF onset lateness from +13.0/+10.5 ms to ≤ +3 ms; impulse density back toward 11.1/s [I]. Risk: start frames get shorter effective support, giving slight LF-resolution loss right before hits.

**F2: Truncated analysis for the start frame (fixes T': pre-echo).** Use an asymmetric analysis window wₐ(n) = hann4096(n)·[n < n0 − 64] and re-estimate the start frame's partials from pre-onset samples only. Predicted: 6–16 kHz excess in [S+512,S+1536] from +11.3 dB to ≤ +3 dB [I]. Risk: wider main lobe means more peak merging in that one frame. Rely on track continuity from the previous frame for frequencies.

**F3: Transient layer under transposition.** Choose one of two options. (a) A time-domain transient pass: cut [n0−Ls/2, n0+~10 ms], resample it by r, and splice with the F1 windows. This is the "remove transient, shift stationary part, reinsert" pattern [L: Nagel & Walther 2009; Driedger, Müller & Ewert, IEEE Signal Processing Letters 21(1):105–109, 2014, where the percussive part goes through short-frame OLA]. (b) Short-frame sinusoids with phase reset at n0 [L: Röbel 2003; Duxbury, Davies & Sandler 2002 phase locking at transients; Quatieri & McAulay 1992 onset-anchored phase].\[7\]\[24\]\[25\]\[26\] Option (a) is simpler and reproduces the unity passthrough benefit under shift.

**F4: Kill near-coincident atoms before rendering (fixes P).**
```cpp
// after joint LS, per frame: cluster atoms with |fi - fj| < 1 bin (11.72 Hz)  -> rho(1 bin)=0.5, cond<=3
for (auto& C : clusters) {
  std::complex<float> z = 0; float p = 0, fw = 0;
  for (auto& a : C) { z += a.amp*std::polar(1.f, a.phase); p += a.amp*a.amp; fw += a.amp*a.amp*a.freq; }
  if (p > 1.26f*std::norm(z)) {               // kappa > +1 dB
    Atom m; m.amp = std::abs(z); m.phase = std::arg(z); m.freq = fw/p;  // power-weighted freq
    replace(C, m);                             // one atom, coherent amplitude
  }
}
// coasting: cap at 3 frames in ALL bands; terminate a coasting track if any matched peak lies within ±1.5 bins of its bin
// LS: enforce min spacing 1 bin before the solve, or use lambda_ij scaled by 1/(1-rho_ij)
```
Predicted: prominent events on 1985 up from 50 to ≤15, and the direction-twin fraction falls; the 1,205/925/530 Hz events drop below +6 dB [I]. Risk: legitimate beating pairs (chorus, piano unisons) merge into one atom. Restrict merging to pairs with κ>1.26, where one atom is coasting or comes from the other tier.

**F5: Frequency-proportional or predictive gate (fixes P').** Δf_max(f) = max(25 Hz, 0.03·f), centred on the prediction f̂ₖ₊₁ = fₖ + (ψ̂/2π)·H/fs. Take ψ̂ from a phase-difference slope estimator [L: Marchand 2012] or LP prediction [L: Lagrange et al. 2007].\[8\]\[9\] Predicted: the 2.04–2.19 s stubs fall below +6 dB [I]. Risk: more false continuations in dense HF. Keep the 25 Hz floor below ~800 Hz.

**F6: Replace the per-sample ceiling.** This is an offline renderer, so apply a static file gain g = 0.95/max|y|, or a true-peak look-ahead limiter with ≥300 ms release. Predicted: 440saw AM from 0.8–2.4 dB to <0.1 dB and level error from −0.8…−1.0 to ≈0 dB [I]. Risk: none audible. Shifted outputs get quieter by the crest-factor increase (≤4.7 dB on saws).

**F7: Tier exclusivity at 200 Hz.** Below 200 Hz, when the LF tier is active, delete long-tier peaks and assign the 180–220 Hz overlap to exactly one tier with hysteresis. Predicted: HappyMono 180–200 Hz from +6.0 dB to ≤ +1 dB [I].

**F8: Residual and HF admission.** Require ≥3-frame track survival above 4 kHz before rendering as a sinusoid, and send everything else to the residual. Replace max(0, |X|−|M|) with over-subtraction plus a spectral floor (Berouti), or temporally smoothed gains (Cappé-style decision-directed smoothing) [L].\[17\]\[22\]\[27\]\[28\]\[29\] Predicted: bulk events from 127/s toward ≤40/s; 4–15 kHz excess from +1.9…+2.7 dB toward ≤ +1 dB [I]. Risk: some HF tonal detail (cymbal partials) moves into noise.

---

## Caveats

- All code-level statements are hypotheses built from [D] and measured signatures. Nothing here was read from main.cpp.
- The Defect-T and Defect-P mechanisms are **very likely, not proven**. E2 and E3 are the atom-level and weight-level proofs. The a₂/a₁≈0.78 figure is an illustrative fit to the measured excess, not a measured atom ratio.
- The asymmetry between up and down in prominent-event counts (1.7–2×) is unexplained (**unable to determine**).
- The cause of the HF +2 dB excess is **possible** only.
- The LF-tier-inactive-on-1985 explanation is an inference from window length versus hit spacing.
- Aliasing, chorus/phasing and decay-time metrics were not measured. Their absence from this report is not evidence that they are absent from the outputs.
- Literature quotes on Edler's PR condition and on AAC window shapes come from patent reproductions of the standard and of Edler's work, not from Edler's paywalled original. The Duxbury et al. description rests on its abstract and on secondary descriptions.

## Sources

1. [\[1902.01053\] Overlap-Add Windows with Maximum Energy Concentration for Speech and Audio Processing](https://ar5iv.labs.arxiv.org/html/1902.01053)
2. [Method and apparatus for encoding and decoding an audio signal using adaptively switched temporal resolution in the spectral domain](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8095359)
3. [Encoder, decoder and methods for backward compatible dynamic adaption of time/frequency resolution spatial-audio-object-coding](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/9734833)
4. [Encoder, decoder and methods for signal-dependent zoom-transform in spatial audio object coding](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10152978)
5. [Low-complexity spectral analysis/synthesis using selectable time resolution](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8392202)
6. [Improved time-scaling of musical audio using phase locking at transients](https://www.researchgate.net/publication/243766616_Improved_time-scaling_of_musical_audio_using_phase_locking_at_transients)
7. [(PDF) Phase Vocoder Done Right](https://www.researchgate.net/publication/358632034_Phase_Vocoder_Done_Right)
8. <https://sylvain-marchand.info/Publications/dafx12.pdf>
9. [Enhanced Partial Tracking Using Linear Prediction - Archive ouverte HAL](https://hal.science/hal-00308184)
10. [\[PDF\] Speech analysis/Synthesis based on a sinusoidal representation](https://www.semanticscholar.org/paper/Speech-analysis-Synthesis-based-on-a-sinusoidal-McAulay-Quatieri/a10b70775dfade81cd7aeea6a193e73764cef5c5)
11. [IMPROVED HIDDEN MARKOV MODEL PARTIAL TRACKING THROUGH TIME-FREQUENCY ANALYSIS](https://www.semanticscholar.org/paper/IMPROVED-HIDDEN-MARKOV-MODEL-PARTIAL-TRACKING-Kereliuk-Depalle/740e96fcf61d89e26c5079ece63b40eae3302134)
12. [Synthesis of low-peak-factor signals and binary sequences with low autocorrelation (Corresp.)](https://www.semanticscholar.org/paper/Synthesis-of-low-peak-factor-signals-and-binary-low-Schroeder/cb75785e6c3efdc3edc2717b6558070e1559d9f2)
13. [\[PDF\] Spectral Modeling Synthesis](https://www.semanticscholar.org/paper/Spectral-Modeling-Synthesis-Serra-Smith/11f17c4e16927cb5a41fad4ce4e751743659de64)
14. [SMS](https://doc.sccode.org/Classes/SMS.html)
15. [Enhanced Fuzzy Decomposition of Sound Into Sines, Transients and Noise](https://arxiv.org/html/2210.14041v2)
16. [Extending Spectral Modeling Synthesis with Transient Modeling Synthesis](https://www.semanticscholar.org/paper/Extending-Spectral-Modeling-Synthesis-with-Modeling-Verma-Meng/aeb40b67a80cf028413333036f3a66641e4045af)
17. [Improved Noise Minimum Statistics Estimation Algorithm for Using in a Speech-Passing Noise-Rejecting Headset](https://link.springer.com/article/10.1155/2010/395048)
18. [On the optimal linear filtering techniques for noise reduction - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0167639307000295)
19. [Noise suppression by spectral magnitude estimation —mechanism and theoretical limits— - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/0165168485900027)
20. [Noise suppressor having weighted gain smoothing](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/6317709)
21. [Experiments with Speech Enhancement Techniques](https://www.ee.columbia.edu/~dpwe/classes/e4810-2002-09/projects/mg2016/report.html)
22. [Noise reduction method](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/7133824)
23. <https://www.audiolabs-erlangen.de/content/resources/aesCodingTutorial/media/tutorial/aes.pdf>
24. [Thomas QUATIERI](https://www.researchgate.net/profile/Thomas-Quatieri/3)
25. [Publications](https://www.ll.mit.edu/r-d/publications?rdgroup=763&page=43)
26. <http://www.aes.org/e-lib/browse.cfm?elib=11324>
27. [A multi-band spectral subtraction method for enhancing speech corrupted by colored noise](https://www.researchgate.net/publication/229005750_A_multi-band_spectral_subtraction_method_for_enhancing_speech_corrupted_by_colored_noise)
28. [On the decision-directed estimation approach of Ephraim and Malah](https://www.researchgate.net/publication/4087373_On_the_decision-directed_estimation_approach_of_Ephraim_and_Malah)
29. [UBERSICHTSARBEIT A review of the Ephraim-Malah noise reduction algorithms](https://www.orl.uzh.ch/dga-ev/publikationen/zfaudiologie/archiv/ZfA_2001_40-1_004-015_Review.pdf)
