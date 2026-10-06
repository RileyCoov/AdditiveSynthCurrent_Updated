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

(filled after the renders; not part of the frozen section)
