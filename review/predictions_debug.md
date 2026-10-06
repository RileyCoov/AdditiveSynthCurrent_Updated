# F1 / F2 / F3 on the Xcode Debug build: predictions and outcomes

Branch `review/short-block-fixes` at f5dd718 (17cddd8 + `fixes.patch`, applied unmodified).
Binary: `review/bin_debug/patched` (Xcode Debug, md5 f8d4f4ff3b3922a4d79ee9ca31e45933).

## FROZEN PREDICTIONS

**Frozen at 2026-10-06T00:29:17-0700.** This section is not edited after this time. Outcomes go in the separate section below.

State when frozen: no switch-ON shifted render exists from this build. Four switch-OFF shifted renders exist (step-3 identity checks: pristine and patched-off, 1985 −7 and DrumLoopShort +7). No metric has been run on any of them.

### Source

Copied from "Pitch shift engine independent review.md" §7.1-7.4 (Predicted column and labels, verbatim), with the full pass-rule text copied from `research_notes/.../fix_prototypes.md` "Predictions and pass rules" (written 2026-10-05T13:53:58-07:00), because the §7 tables only print the outcome of each rule. The "Release (review)" column is the review's own measurement on its cmake Release build, kept here only as the reference for the Debug-vs-Release check. It is not a prediction about this build.

Abbreviations: dn5 / up5 / dno = −7 / +7 / −12. base = patched build, all switches off.

### F1 `SHORT_BASS_AMP_HOLD_HZ=400` (§7.1)

| Item | Label | Predicted | Pass rule (fix_prototypes.md) | Release (review) |
|---|---|---|---|---|
| P1 1985 M1_dip, dn5 / up5 / dno | informed | −1.4 / −1.2 / −0.6 ±0.6 | −2.0 ≤ M1_dip ≤ +1.0 on all three | −1.40 / −1.19 / −0.56 (base −5.61 / −4.96 / −4.81) |
| P2 1985 OV 160-400 Hz | unlabelled | −1.0 ±1.0 | within [−2.5, +1.5] | −2.47 / −1.25 / +0.12 (base −4.25 / −4.50 / −3.76) |
| P2 1985 OV 400-2000 Hz vs base | unlabelled | ±0.7 | not lower than BASE − 1.0 | +0.35 / +0.11 / +0.09 |
| P3 DrumLoopShort kick RISE, Δ vs base, up5 / dn5 | blind | +3 ±3 ms; "expect FAIL" (~55%) | median Δ_F1 ≤ median Δ_BASE + 2.0 ms on BOTH up5 and dn5 | +1.33 / 0.00 ms (base +11.33 / +32.67) |
| P4 HappyMono dn5 M1_dip | blind | [−2.5, +1.0] | inside [−2.5, +1.0] | −0.85 (base −4.75) |
| P5 FairC2 CLICK, within ±3 dB of base | blind | ±3 dB | inside ±3 dB of BASE | dn5 17.27 → 10.14; up5 20.13 → 15.94 |
| P6 Piano dn5 RISE, LF / HF | blind | ±3 ms | not later than BASE + 2 ms | 0.00 / 0.00 ms |
| P7 1985 HFATT, dn5 / up5 / dno | unlabelled | ±150 smp | within ±150 smp of BASE | 0 / 0 / 0 |
| P7 1985 kick RISE not earlier than base − 1 ms | unlabelled | — | not earlier than BASE − 1 ms | dn5 −3.4, up5 +0.7, dno +1.3 ms |

### F1+F2 combined (§7.2)

COMBO rules: "M1_dip and OV as that fix alone ±0.7 dB; CLICK as F2 alone ±2 dB; kick RISE as that fix alone ±1 ms. PASS = all of that fix's and F2's pass rules hold together."

| Item | Predicted | Pass rule | Release (review) |
|---|---|---|---|
| M1_dip vs F1 ±0.7 | — | F1 P1 rule [−2.0, +1.0] | −0.89 / −0.51 / −1.94 |
| OV 160-400 vs F1 ±0.7; rule [−2.5, +1.5] | — | [−2.5, +1.5] | −3.43 / −1.26 / −2.66 |
| CLICK vs F2 ±2 dB | — | F2 P1 rule (drop ≥ 6 dB vs base on both, ≤ +10 dB) | Fair up5 5.61, dn5 3.52; Drum 9.45 / 7.82 |
| Drum kick RISE vs F1 ±1 ms; rule ≤ base + 2 ms | — | median ≤ BASE + 2.0 ms, up5 and dn5 | up5 25.00, dn5 27.17 ms (base +11.33 / +32.67) |

### F2 `SHORT_BLOCK_8=1` (§7.3)

| Item | Label | Predicted | Pass rule (fix_prototypes.md) | Release (review) |
|---|---|---|---|---|
| P1 FairC2 CLICK, up5 / dn5 | blind | drop ≥6 dB, ≤ +10 dB | median falls by ≥ 6 dB on BOTH up5 and dn5 (value ≤ +10 dB) | 20.13 → 5.79; 17.27 → 3.48 (original 4.60); max abs d2 1403 → 270, 959 → 76 LSB |
| P2 1985 M1_dip, dn5 / up5 / dno | informed | −5.3 / −5.9 ±0.5; −5.0 ±0.7 | no harm: M1_dip ≥ BASE − 1.2 dB on all three | −5.22 / −5.97 / −5.55 |
| P3 DrumLoopShort kick RISE, Δ vs base | unlabelled | ±2 ms | ≤ BASE + 2 ms | up5 −8.8, dn5 −5.5 ms |
| P4 1985 HFATT, dn5 / up5 / dno | unlabelled | ±150 smp | within ±150 smp of BASE | −20 / +16 / +32 |
| P4 Piano RISE, LF / HF | unlabelled | ±2 ms | within ±2 ms of BASE | 0.00 / +3.00 ms (n=1) |

### F3 `ATOM_MERGE_BINS=1.0` (§7.4)

Pass rule (fix_prototypes.md): "PASS if [BIRDIES] reduced by ≥ 10 % on BOTH, with the 40-400 Hz blob count not up by > 20 % and M1_dip within ±0.5 dB of BASE."

| Item | Label | Predicted | Release (review) |
|---|---|---|---|
| BIRDIES 400-3500 Hz, dn5 / up5 | blind | −15 ±15% | 21.09 → 9.35 /s (−56%); 10.00 → 3.91 /s (−61%) |
| 40-400 Hz blob count not up by >20% | blind | — | dn5 4 → 5 (+1 blob); up5 6 → 2 |
| M1_dip within ±0.5 dB of base | blind | — | −5.61 / −5.04 (base −5.61 / −4.96) |

### Unity (all)

With each switch ON, the 0-semitone render is byte-identical to switch-OFF unity from the same binary. (Already checked in step 3, before this file was written: pass on 1985, DrumLoopShort, FairC2, HappyMono, PianoSampleMono for F1, F2, F3, F1+F2, F1+F2+F3.)

### Debug-vs-Release flag rule (fixed now, before measuring)

For every value the review reports, flag it when |Debug − Release| exceeds that item's own stated ± tolerance from the tables above. Baseline values use the tolerance of the item they belong to. Where an item states no ± interval, I use these fallbacks, chosen now:

| Metric | Fallback tolerance | Taken from |
|---|---|---|
| M1_dip | ±0.5 dB | F2 P2 / F3 (tightest stated M1 interval) |
| OV (any band) | ±0.7 dB | F1 P2 400-2000 / COMBO |
| kick RISE | ±2 ms | F2 P3 |
| CLICK dB | ±2 dB | COMBO |
| CLICK max abs d2 | ±25 % of the Release value | none stated; my choice |
| HFATT | ±150 smp | F1 P7 / F2 P4 |
| BIRDIES /s | ±15 % of the Release value | F3 interval half-width |
| 40-400 Hz blobs | ±20 % of the Release value, minimum ±1 blob | F3 rule |

## OUTCOMES

Written 2026-10-06 after all renders and metrics. Binary `review/bin_debug/patched` (Xcode Debug). Metrics come from `review/scripts_debug/fixmetrics.py`, which is the review's script with path changes only (`scripts_debug/fixmetrics_path_changes.diff`). Values are in `logs_debug/metrics.jsonl`, and every render's env and md5 is in `logs_debug/renders.tsv`. Scoring is in `scripts_debug/score.py`. All values below are [M].

"Interval" says whether the Debug value fell inside the predicted interval. "Rule" is the pre-registered pass rule. None of these results is an audible verdict.

### Predicted vs measured (Debug)

| Fix | Item | Predicted | Debug measured | Interval | Rule |
|---|---|---|---|---|---|
| F1 | P1 1985 M1_dip dn5/up5/dno | −1.4/−1.2/−0.6 ±0.6 | -1.39/-1.19/-0.56 (base -5.61/-4.96/-4.81) | hit | PASS |
| F1 | P2 1985 OV 160-400 dn5/up5/dno | −1.0 ±1.0 | -2.47/-1.25/+0.12 (base -4.25/-4.50/-3.76) | miss dn5, dno | PASS |
| F1 | P2 1985 OV 400-2000 Δ vs base | ±0.7 | +0.34/+0.11/+0.09 | hit | PASS |
| F1 | P3 Drum kick RISE Δ up5/dn5 | +3 ±3 ms; "expect FAIL" | +1.33/+0.00 ms (base +11.33/+32.67) | hit | PASS |
| F1 | P4 HappyMono dn5 M1_dip | [−2.5, +1.0] | -1.20 (base -5.04) | hit | PASS |
| F1 | P5 FairC2 CLICK within ±3 dB of base | ±3 dB | dn5 17.27 → 10.14; up5 20.13 → 15.94 | miss (lower) | **FAIL** |
| F1 | P6 Piano dn5 RISE Δ LF/HF | ±3 ms | +0.00/+0.00 ms | hit | PASS |
| F1 | P7 1985 HFATT Δ dn5/up5/dno | ±150 smp | +0/+0/+0 | hit | PASS |
| F1 | P7 1985 kick RISE Δ not earlier than −1 ms | — | -3.4/+0.7/+1.3 ms | — | **FAIL (dn5)** |
| F1+F2 | M1_dip vs F1 ±0.7; rule [−2, +1] | — | -0.89/-0.51/-1.94 | miss dno | PASS |
| F1+F2 | OV 160-400 vs F1 ±0.7; rule [−2.5, +1.5] | — | -3.43/-1.26/-2.66 | miss dn5, dno | **FAIL (dn5, dno)** |
| F1+F2 | CLICK vs F2 ±2 dB | — | Fair up5 5.61, dn5 3.52; Drum 9.45/7.83 | hit | PASS |
| F1+F2 | Drum kick RISE vs F1 ±1 ms; rule ≤ base + 2 ms | — | up5 25.00, dn5 27.17 ms (base +11.33/+32.67) | miss | **FAIL (up5)** |
| F2 | P1 FairC2 CLICK up5/dn5 | drop ≥6 dB, ≤ +10 dB | 20.13 → 5.79; 17.27 → 3.48 (orig 4.60); max|d2| 1403 → 270, 959 → 76 LSB | hit | PASS |
| F2 | P2 1985 M1_dip dn5/up5/dno | −5.3/−5.9 ±0.5; −5.0 ±0.7 | -5.22/-5.97/-5.55 | hit | PASS (no harm) |
| F2 | P3 Drum kick RISE Δ up5/dn5 | ±2 ms | -8.8/-5.5 ms | miss (earlier) | PASS |
| F2 | P4 1985 HFATT Δ dn5/up5/dno | ±150 smp | -20/+16/+32 | hit | PASS |
| F2 | P4 Piano RISE Δ LF/HF (n=1) | ±2 ms | +0.00/+3.00 ms | miss HF | **FAIL (HF)** |
| F3 | BIRDIES 400-3500 Hz dn5/up5 | −15 ±15% | 21.09 → 9.35 /s (-56%); 10.00 → 3.91 /s (-61%) | miss (bigger) | PASS |
| F3 | 40-400 Hz blob count not up >20% | — | dn5 4 → 5; up5 6 → 2 | — | **FAIL (dn5)** |
| F3 | M1_dip within ±0.5 dB of base | — | -5.61/-5.04 (base -5.61/-4.96) | hit | PASS |

Margins worth knowing: F1 P2 passes dn5 by 0.03 dB (−2.47 against the −2.5 limit), the same margin as Release. The F1 P3 PASS sits below RISE's own resolution (report §7.1: median shifts under ~10 ms are not resolvable at n=7).

### Debug vs Release

I compared all 555 shared values in the review's Release `metrics.jsonl` (40 tags) with the Debug values, using the tolerances frozen above. 318 float values agree to within 0.005.

Beyond tolerance, HappyMono dn5 only:

| Value | Release | Debug | Δ | Tolerance |
|---|---|---|---|---|
| base lag (smp) | −121 | −179 | −58 | ±48 |
| base OV 30-160 | −3.85 | −4.56 | −0.71 | ±0.7 |
| base kick RISE median (ms, no rule) | −0.67 | +10.83 | +11.5 | ±2 |
| base kick RISE per hit | of 15 hits, 5-13 differ (e.g. hit 5: Release −10.7, Debug +38.8) | | | |
| F1 OV 30-160 ("+2.75 dB overshoot", §7.1 Costs) | +2.75 | +1.39 | −1.36 | ±0.7 |
| F1 kick RISE per hit | hit 5 Release +18.3 / Debug +38.8, hit 15 Release −34.7 / Debug −26.3; the others are within 1.2 ms (median within tolerance) | | | |

This matches report Appendix B.1: HappyMono is the file where Release and Debug diverge (from 5.078 s). Every pre-registered item gets the same verdict on both builds.

Render md5s: 13 of 55 shared renders are byte-identical between Debug and Release (all FairC2 except base/F2 dn5, plus base/F1/F2 Piano dn5). The other 42 differ at sample level, but their metrics agree.

### Lag check (pre-registration: |lag| ≤ 48 smp required, else reported)

base_1985_dno +59, base_happy_dn5 -179, cmb_happy_dn5 -534, f1_1985_dno +110, f1_happy_dn5 -317, f1b_drum_dn5 -456, f2_drum_dn5 +62, f2_happy_dn5 -264, f3_1985_dno +88, f3_drum_dn5 +60, f3_drum_up5 +83, f3_happy_dn5 -353

The Release run shows the same out-of-range set (1985 dno base +59 / F1 +110, HappyMono −121..−534, Drum dn5 F1b −456 / F2 +62). The report (B.3) says this estimator is unreliable at −12 and on HappyMono. F3 adds 1985 dno +88 and Drum ±7 +60 / +83; those have no Release run to compare against.

### Renders with no prediction and no Release value (F3 beyond 1985 ±7)

These were requested here but not run by the review. They are reported without a verdict:

| Render | M1_dip (base → F3) | OV 30-160 | OV 160-400 | Other |
|---|---|---|---|---|
| 1985 dno | −4.81 → −4.55 | −5.03 → −3.39 | −3.76 → −2.49 | OV 400-2000 −4.08 → −3.37 |
| Drum dn5 | −3.89 → −5.15 | −4.97 → −4.96 | −0.93 → −2.13 | CLICK 13.98 → 11.11 dB |
| Drum up5 | −4.17 → −4.30 | −4.92 → −5.89 | −3.34 → −4.04 | CLICK 10.57 → 9.19 dB |
| Fair dn5 | −5.37 → −5.19 | ≈ | ≈ | CLICK unchanged 17.27; d2 959 → 881 |
| Fair up5 | −5.33 → −5.34 | ≈ | ≈ | CLICK unchanged 20.13 |
| HappyMono dn5 | −5.04 → −4.94 | −4.56 → −4.27 | −4.18 → −4.75 | |
| Piano dn5 | −0.73 → −0.04 | | | RISE LF/HF unchanged |

Piano M1/OV values are not meaningful: they come from the decay-tail block (report §7.1).
