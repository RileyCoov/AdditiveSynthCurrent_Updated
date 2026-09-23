//
//  SaveAddititve.h
//  AdditiveSynthChange
//
//  Created by Riley on 2/20/25.
//
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>   // uint8_t / int32_t: Apple clang pulls these in transitively; GCC does not
#include <algorithm>
#include <complex>
#include <fstream> // This is for the binary process
using namespace std;

// The addition of pragma ensures there isn't padding happening which could
//cause more memory than wanted to be used
struct BinaryFrameHeader {
    int32_t frame_index;
    int32_t start;
    int32_t stop;
    int32_t frame_size;
    uint8_t trans;
    int32_t hop_size;
    size_t window_size;
};

class PeakTrack {
public:
    int id;
    double freq_hz;
    double max_db;
    double current_db;
    int peak_bin;
    double phase;
    bool alive;
    bool edit;
    double freq_parabolic = 0.0;
    int analysis_fft_size = 4096; // FFT size that produced current_db (for correct scaling)

    // ===== FIX C: birth-confirmation gate =====
    // matched_count counts how many times this peak has been seen as a match.
    // Born peaks start at 1; each subsequent match increments. Once it reaches
    // peak_birth_confirm_frames, `confirmed` flips true and the peak is allowed
    // into frame_info for synthesis. Phantom peaks that flash for one frame
    // never reach confirmation and are silently discarded.
    int matched_count = 1;
    bool confirmed = false;
    // ===== END FIX C =====
    // Consecutive frames carried by the unmatched-coast path (reset on every
    // real match). Sub-LF-cutoff tracks are killed after a few coasts: the LF
    // spectrum is dense enough (bass + kick) that a dead note's bin keeps real
    // energy indefinitely, so uncapped coasting lets stale bass tracks pile up
    // until the output limiter crushes the whole file.
    int coast_count = 0;
    // Last two RAW frequency measurements (t-1, t-2), for the median-of-3
    // trajectory filter in the long matched path. Median over raw values so
    // the filter never feeds back on itself; 0 = not yet populated.
    // freq_slope_hz_s: this track's frequency rate of change at this frame, Hz
    // per second, filled in after tracking by a centred difference over the
    // track's own trajectory. Drives the linear-FM (chirp) term in both the
    // joint solve and synthesis -- see chirp_mode in main.cpp. 0 = stationary.
    double freq_slope_hz_s = 0.0;
    double freq_hist1 = 0.0;
    double freq_hist2 = 0.0;
    // Consecutive frames the measured magnitude has fallen (reset on any
    // rise). Distinguishes a real decay (monotone fall for many frames)
    // from vibrato/analysis ripple (alternating) for the release choice.
    int fall_streak = 0;

    PeakTrack(int _id, double _freq, double _mag, int _peak_bin, double _phase, int _analysis_fft_size = 4096) : id(_id), freq_hz(_freq), max_db(_mag), current_db(_mag), peak_bin(_peak_bin), phase(_phase), alive(true), edit(true), analysis_fft_size(_analysis_fft_size) {}
};


void save_binary(const std::string& filename, int num_frames,
                 std::vector<SynthInformation> containsSynthPlacement,
                 std::vector<std::vector<PeakTrack>> frames_peaks);
void read_binary(const std::string& filename, int num_frames,
                 std::vector<SynthInformation>& containsSynthPlacement,
                 std::vector<std::vector<PeakTrack>>& frames_peaks);

