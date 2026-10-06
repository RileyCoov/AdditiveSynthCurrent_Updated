#!/bin/bash
# render.sh <tag> <bin: patched|pristine> <name> <semitones> [KEY=VAL ...]
# Xcode Debug binaries snapshotted to review/bin_debug; TRANS_DEBUG=1 always (log feeds fixmetrics segments).
# Every engine knob is scrubbed from the environment; only the KEY=VAL args given here reach the binary.
R=/Users/Riley/Downloads/AdditiveSynthCurrent-LongWindowFreqUsedForShortWindows/review; A=$HOME/Desktop/comparingOut/EXPERT_REVIEW_2026-10-04/audio/original
tag=$1; b=$2; name=$3; semi=$4; shift 4
BIN=$R/bin_debug/$b
out=$R/renders_debug/$tag.wav
t0=$(date +%s)
env -i HOME=$HOME PATH=/usr/bin:/bin TRANS_DEBUG=1 "$@" $BIN $A/${name}_original.wav $out 1024 $semi > $R/logs_debug/$tag.log 2>&1; rc=$?
md5=$(md5 -q $out)
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%ss\n" "$(date +%H:%M:%S)" "$tag" "$b" "$name" "$semi" "$*" "$rc" "$md5" "$(( $(date +%s)-t0 ))" >> $R/logs_debug/renders.tsv
echo "$tag rc=$rc $md5"
