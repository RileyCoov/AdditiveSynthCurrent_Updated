# Build determinism and the threaded joint solve (2026-11-08)

The reference renders (the `*ADD.wav` files judged by ear) were made with the
Xcode **Debug** configuration at optimization level 0. Every other build must
reproduce them **byte for byte**; a null test (render minus reference) is the
check, not listening. `cmp` on the two WAVs is enough.

## What changes the bits, and what does not

Measured on SaintSaensMonoChord, HappyMono and DrumLoopShort with the standard
`<in> <out> 16` invocation, Apple clang 17:

| Build | Result vs Debug |
|---|---|
| `-O0` (Xcode Debug as it was) | reference |
| `-O1 -fno-builtin` | identical on all three, after the RealFFT.cpp fix below |
| `-O1` alone | differs: the optimizer merges `sin(x)`+`cos(x)` into Apple's `__sincos_stret`, whose last bit differs from separate calls |
| `-O2`, `-O3`, `-Os` (with or without `-fno-builtin`) | differ, 1 LSB on most samples but up to ~77 LSB where a rounding flip changes a peak-detection decision (null peak about -53 dBFS) |
| `-ffp-contract=off` | differs: the Debug build already uses fused multiply-add, so it is part of the reference |

Sanitizers (address, undefined behaviour) and zero/pattern stack initialisation
all leave the Debug output unchanged, so none of the differences are bugs; they
are rounding-order changes.

**RealFFT.cpp**: `InvRealFFT` computes its twiddles in `float`. `cos(a)` on a
float resolves to libc++'s inline wrapper around a compiler intrinsic, which
`-fno-builtin` does not cover, so the sincos merge still happened there. The
twiddles now call `cosf`/`sinf` directly. That is the same library routine the
Debug build was already calling, so the Debug output is unchanged, and the
calls stay separate in optimized builds.

## Settings

- Xcode: both configurations set `GCC_OPTIMIZATION_LEVEL = 1` and
  `OTHER_CPLUSPLUSFLAGS = -fno-builtin`.
- CMake: `CMAKE_CXX_FLAGS_RELEASE` is `-O1`, plus `-fno-builtin` on the target.
  The battery baselines recorded before this date were rendered at `-O3` and
  may show tiny metric shifts against a new run; re-baseline once.

Speed on an M-series Mac (8 P + 4 E cores), output bit-identical in every row:

| Build | SaintSaens (3.5 s audio) | HappyMono (14 s audio) |
|---|---|---|
| Debug `-O0`, before | 11.7 s (14.3 s from Xcode) | 73.9 s |
| `-O1 -fno-builtin`, single thread | 4.5 s | 27.9 s |
| `-O1 -fno-builtin`, threaded (shipped) | 3.0 s | 18.6 s |

## Threading in joint_amp_phase

Profiling the optimized build put 1.76 s of the 4.65 s SaintSaens run in the
right-hand-side loop of `joint_amp_phase` (328 calls, 183k partials, each an
oscillator recurrence over 4096 samples). The pair and conjugate-gradient parts
are under 0.2 s together and are left sequential.

The loop now hands partials to `std::thread` workers in chunks of 4 from an
atomic counter. Each partial's two sums are accumulated in the same order on
whichever thread runs them and each writes only its own `rhs` slots, so the
result is bit-identical. `JOINT_THREADS=1` disables it; any other value sets
the worker count (default: hardware thread count). Frames with fewer than 64
partials run inline.

### Porting

`std::thread` and `std::atomic` are standard C++11 and build unchanged on any
ARM target with an operating system (Linux, Android, macOS/iOS). The speedup
scales with the number of cores actually available; on a single-core ARM
(an ADSP-SC5xx Cortex-A5/A55, for instance) set `JOINT_THREADS=1` or let the
threshold keep it inline, and on a bare-metal toolchain without a threading
runtime replace the `std::thread` block with the sequential loop under an
`#ifdef`. Bit-identity across *compilers* is a separate matter: GCC defaults to
`-ffp-contract=fast` and its libm differs from Apple's, so an ARM Linux build
will not match the Mac renders bit for bit regardless of threading. Judge the
port with the null test against the Mac reference and expect a residual of the
same order as the -O3 case above.

## Shifted renders changed on 2026-11-24 (by design)

The cluster phase lock's member rule now carries the measured member-minus-
dominant offset over to the shift (main.cpp, CLUSTER PHASE LOCK block, `cl_link`):
a member keeps its own shifted frequency and only its relation to the dominant is
held. Before, the offset advanced at the unshifted difference, which put a
member at the dominant's shifted frequency plus that difference (SaintSaens up a
fifth: the pedal's 60 and 79 Hz partials at 80 and 98 Hz instead of 90 and 118,
each with a grain-rate sideband). Chosen by ear against Re-Synth's oscillator,
which does the same. Unity renders are untouched (the block is shift-only), so
the Debug references above still hold at unity; the shifted references made
before this date differ from the current output in the member partials only.
