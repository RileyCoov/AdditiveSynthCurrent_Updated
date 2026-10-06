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

### Debug vs Release beyond tolerance

| Value | Release (review) | Debug | Δ | Tolerance |
|---|---|---|---|---|
| base_happy_dn5 OV_30_160 (other effect) | -3.85 | -4.55839 | -0.71 | ±0.7 |
| f1_happy_dn5 OV_30_160 (other effect) | 2.75 | 1.38815 | -1.36 | ±0.7 |

### Lag check (|lag| ≤ 48 smp required by the pre-registration, else reported)

base_1985_dno +59, base_happy_dn5 -179, cmb_happy_dn5 -534, f1_1985_dno +110, f1_happy_dn5 -317, f1b_drum_dn5 -456, f2_drum_dn5 +62, f2_happy_dn5 -264, f3_1985_dno +88, f3_drum_dn5 +60, f3_drum_up5 +83, f3_happy_dn5 -353
