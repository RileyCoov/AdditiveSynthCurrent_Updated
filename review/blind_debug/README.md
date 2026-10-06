# Blind A/B pairs from the Xcode Debug build

All 34 files come from one binary: `review/bin_debug/patched`, which is the Xcode Debug build of branch `review/short-block-fixes` (17cddd8 + fixes.patch). Switches are environment variables. They are off by default and forced off at 0 semitones.

The comparisons are the same 17 as `work2/fixes/README.md`. The A/B sides were randomised again with a fresh seed, so they are unrelated to the review's Release pairs.

| Switch | Fix |
|---|---|
| `SHORT_BASS_AMP_HOLD_HZ=400` | F1 |
| `TSA_SPLIT_HZ=400` | F1b (existing switch, no new code) |
| `SHORT_BLOCK_8=1` | F2 |
| `ATOM_MERGE_BINS=1.0` | F3 |

"F1+F2" means `SHORT_BASS_AMP_HOLD_HZ=400 SHORT_BLOCK_8=1`.

## How to listen

- `pairNN_A.wav` and `pairNN_B.wav` are full renders of the same file at the same shift. One is the reference and the other is the fix, with sides assigned at random.
- The key is `review/.blind_key_debug.tsv`. Please don't open it until you have recorded your verdicts.
- Don't md5 the pair files or compare them with `review/renders_debug/`, which holds the same renders under named configs.
- The pairs are **not level-matched**. Level changes inside the drum-hit blocks are part of what the fixes do.
- For each pair, record one of: A better / B better / tie. Add a word on what you heard.
- All times below are output times. S is the start of a transient block. The short-frame block runs from S+31 ms to S+52 ms, and the S+2496 junction is at S+52 ms.

## Pairs

Tier 1 decides the main question: does removing the bass dip remove the "continuous musical noise", and does it cost the kick?

| Pair | File / shift | What is compared | What to listen for |
|---|---|---|---|
| 01 | 1985 down_fifth | baseline vs F1+F2 | The disturbance at every drum hit (every 0.256 s): a bass "chop" or dropout, and a pitch blip. Clean hits on steady bass are at 5.74, 6.76 and 7.81 s. Also listen for any added thump or boom at the hits |
| 02 | 1985 down_fifth | baseline vs F1 | Same as 01 (F1 alone, without the splice fix) |
| 03 | DrumLoopShort up_fifth | baseline vs F1 | **Kick timing and weight**, especially the first hit (0.28-0.33 s). Does the low thump land with the stick, or late and soft? |
| 04 | DrumLoopShort up_fifth | baseline vs F1+F2 | Same as 03. Also listen for a click or tick at the end of each hit |
| 05 | 1985 up_fifth | baseline vs F1+F2 | Same as 01, in the up direction |

Tier 2 covers the splice and click fix, the birdie fix, and side effects.

| Pair | File / shift | What is compared | What to listen for |
|---|---|---|---|
| 06 | FairlightAdditiveWaveTA1_C2 up_fifth | baseline vs F2 | A click at each timbre switch: 0.095, 0.201, 0.756, 1.993, 2.548, 3.785, 4.148 s |
| 07 | FairlightAdditiveWaveTA1_C2 down_fifth | baseline vs F2 | Same as 06 |
| 08 | 1985 down_fifth | **F1** vs F1+F2 | Once the dip is gone, is the splice still audible? Listen for a small tick or thump at the end of each hit block |
| 09 | 1985 down_fifth | baseline vs F3 | Short tonal chirps ("birdies", under 40 ms) in the 400-3500 Hz range, all through the song. Also listen for hollowness or loss of brightness |
| 10 | 1985 up_fifth | baseline vs F3 | Same as 09. The mid band is rougher on up |
| 11 | DrumLoopShort down_fifth | baseline vs F1 | Kick timing and weight (as in 03) |
| 12 | HappyMono down_fifth | baseline vs F1+F2 | Bass dips at the drum hits, and whether the kick gets boomy |
| 13 | FairlightAdditiveWaveTA1_C2 up_fifth | baseline vs F1 | At each timbre switch (times as in 06): the click, and whether the low/mid of the old timbre hangs on for ~20 ms |
| 14 | 1985 down_octave | baseline vs F1+F2 | As 01 (octave) |

Tier 3 is the alternative to F1, plus an onset control.

| Pair | File / shift | What is compared | What to listen for |
|---|---|---|---|
| 15 | 1985 down_fifth | baseline vs F1b | As 01. Also listen for the bass swelling slightly at the hits |
| 16 | DrumLoopShort down_fifth | baseline vs F1b | Kick: a low "pre-thump" before the stick, or a kick that is too boomy |
| 17 | PianoSampleMono down_fifth | baseline vs F1+F2 | Hammer attack at 0.04 s (should be unchanged), and the very end of the decay |
