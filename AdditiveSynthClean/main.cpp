//
//  main.cpp
//  AdditiveSynthFreqMask
//
//  Created by Riley on 12/16/24.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <string>
#include <sstream>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "RealFFT.h"
#include "AnalysisInfo.h"
#include "Phaser.h"
#include "SampleNormalizer.h"
#include "Windows.h"
#include "UtilityFuncs.h"

#include "SaveAddititve.h"
using namespace UtilityFuncs;
using namespace std;

int LONG_SIZE = 2048;
int SHORT_SIZE = 256;
const double PI = 3.14159265358979323846;



/**
 * Definition of all the windows being used in the file right here
 */
TransitionWindows windowConverter(LONG_SIZE, SHORT_SIZE);
NormalWindows normalWindowObject;
vector<float> long_hanning = normalWindowObject.HanningWindow(LONG_SIZE);
vector<float> short_hanning = normalWindowObject.HanningWindow(SHORT_SIZE);
vector<float> rectangular_window = normalWindowObject.RectangularWindow(LONG_SIZE);
vector<float> rect_fade_to_hann = windowConverter.RectToHann(LONG_SIZE);
vector<float> rect_fade_to_hann_short = windowConverter.RectToHann(SHORT_SIZE);
vector<float> long_to_short = windowConverter.createLongToShortWindow(long_hanning, short_hanning);
vector<float> short_to_long = windowConverter.createShortToLongWindow(long_to_short);

// 4096-sample Hanning window for the analysis FFT.
// The synthesis still uses 2048-sample overlap-add windows for time resolution,
// but the analysis uses 4096 actual samples for true frequency resolution
// (11.72 Hz/bin vs 23.44 Hz/bin), which resolves low harmonics like 42/84 Hz.
static const int ANALYSIS_SIZE = 4096;
vector<float> analysis_hanning = normalWindowObject.HanningWindow(ANALYSIS_SIZE);

// 8192-sample Hanning for the low-frequency analysis FFT (see lf_cutoff_hz in
// the user settings). 2.93 Hz/bin cleanly resolves bass partials ~20 Hz apart that the
// 4096 window (whose Hann main lobe spans ~47 Hz) smears into one beating lobe.
static const int LF_ANALYSIS_SIZE = 16384; // 8192 still let 39.6 Hz mask 59.7 Hz
vector<float> lf_analysis_hanning = normalWindowObject.HanningWindow(LF_ANALYSIS_SIZE);


/*
 The following code is the code that allows for reading and writing of wavs
 */
struct AudioBuffer
{
    // Float sample stream normalized between 1.0 and -1.0.
    // ie - Stereo is stored as mSamples[Left][mNumFrames] and mSamples[Right][mNumFrames]
    float** mSamples;
    int    mChannels; // 1 = mono, 2 = stereo, etc
    long mNumSamples; // total length of the audio stream in samples per channel
    long mSampleRate;
};

struct AppSettings
{
    float sampleRate;
    int blockSize; // smoothing frequency in Hertz
    string inputWavFilePath;
    string outputWavFilePath;
};

static bool readInWaveFile(const string& waveFile, AudioBuffer* buff);
static void writePCM16WaveFile(const string& waveFilePath, float** samples, size_t numSamples, short numChannels, int sampleRate);
static bool parseArgs(int argc, const char* argv[], AppSettings& settings);
static void printUsage();

double wrap_phase(double x) {
    while (x > M_PI) x -= 2.0 * M_PI;
    while (x < -M_PI) x += 2.0 * M_PI;
    return x;
}

std::vector<std::vector<float>> audioBufferToVector(const AudioBuffer& buff)
{
    std::vector<std::vector<float>> audioData;
    audioData.resize(buff.mChannels);

    for (int ch = 0; ch < buff.mChannels; ++ch) {
        audioData[ch].resize(buff.mNumSamples);
        for (int i = 0; i < buff.mNumSamples; ++i) {
            audioData[ch][i] = buff.mSamples[ch][i];
        }
    }
    return audioData;
}

/**
 *These are the contents saved frame by frame for the audio signal. During
 *implemenation thrown into "active_peaks" which holds onto the contents of
 *importance frame by frame, and all collected in a vector<vector<active_peaks>>
 */
/*
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
     
     PeakTrack(int _id, double _freq, double _mag, int _peak_bin, double _phase) : id(_id), freq_hz(_freq), max_db(_mag), current_db(_mag), peak_bin(_peak_bin), phase(_phase), alive(true), edit(true) {}
 };
 */


/**
*tracks the transients, we compare previous frame with the current frame and
 subtract differences *
 */
class TransientDetector {
public:
    vector<float> mCurrentFrame;
    vector<float> mPrevFrame;
};

/**
 *Pretty simple appraoch to detecting pitches. this is the first means of
 *appoaching the signal which gives us our base line detected peaks aboce the
 *threshold. After this there is further refinement*
 */
// Phase 1: extra detection sensitivity (in dB) as a function of frequency.
// Returns 0 dB below f_lo, ramping linearly to max_extra_db by f_hi. Used to
// keep quiet upper harmonics alive (they otherwise fall under the per-frame
// floor / threshold and darken harmonically-rich tones).
static inline double hf_sensitivity_db(double freq_hz, double max_extra_db,
                                       double f_lo = 1500.0, double f_hi = 8000.0) {
    if (max_extra_db <= 0.0 || freq_hz <= f_lo) return 0.0;
    if (freq_hz >= f_hi) return max_extra_db;
    return max_extra_db * (freq_hz - f_lo) / (f_hi - f_lo);
}

vector<int> detect_peaks(vector<double>& magnitude, double threshold,
                         int sr = 48000, int fft_size = 0, double hf_extra_db = 0.0) {
    vector<int> peaks;
    for (size_t i = 1; i < magnitude.size()-1; i++) {
        double thr = threshold;
        if (hf_extra_db > 0.0 && fft_size > 0) {
            double f = (double)i * sr / (double)fft_size;
            thr = threshold * pow(10.0, -hf_sensitivity_db(f, hf_extra_db) / 20.0);
        }
        if (magnitude[i] > magnitude[i-1] && magnitude[i] > magnitude[i+1] && magnitude[i] > thr) {
            peaks.push_back((int) i);
        }
    }
    return peaks;
}

static void filter_peaks_by_quality(const vector<double>& mag, vector<int>& peaks,
                                    double sidelobe_attenuation_dB,
                                    double floor_below_max_dB,
                                    int sr = 48000, int fft_size = 0,
                                    double hf_extra_floor_db = 0.0,
                                    int main_lobe_radius = 4) {
    if (peaks.empty()) return;

    // B1: per-frame floor.
    double max_peak_mag = 0.0;
    for (int p : peaks) {
        if (mag[p] > max_peak_mag) max_peak_mag = mag[p];
    }
    double floor_mag = max_peak_mag * pow(10.0, -floor_below_max_dB / 20.0);

    vector<int> after_floor;
    after_floor.reserve(peaks.size());
    for (int b : peaks) {
        double fm = floor_mag;
        if (hf_extra_floor_db > 0.0 && fft_size > 0) {
            double f = (double)b * sr / (double)fft_size;
            fm = max_peak_mag * pow(10.0, -(floor_below_max_dB + hf_sensitivity_db(f, hf_extra_floor_db)) / 20.0);
        }
        if (mag[b] >= fm) after_floor.push_back(b);
    }
    if (after_floor.empty()) { peaks.clear(); return; }

    //amplitude-aware sidelobe-exclusion. Sort by magnitude descending,
    //accept greedily, but block surrounding bins only at the sidelobe level.
    vector<int> by_mag = after_floor;
    sort(by_mag.begin(), by_mag.end(),
         [&mag](int a, int c) { return mag[a] > mag[c]; });

    int last_idx = (int)mag.size() - 1;
    // block_threshold[k] = the highest sidelobe-attenuated level claimed at
    // bin k by any already-accepted peak.  A candidate at bin k is rejected
    // only if mag[k] < block_threshold[k].
    vector<double> block_threshold(mag.size(), 0.0);
    double sidelobe_ratio = pow(10.0, -sidelobe_attenuation_dB / 20.0);

    vector<int> kept;
    kept.reserve(by_mag.size());
    for (int b : by_mag) {
        if (mag[b] < block_threshold[b]) continue; // looks like a sidelobe of an accepted peak

        kept.push_back(b);
        double claim_level = mag[b] * sidelobe_ratio;
        int lo = std::max(0, b - main_lobe_radius);
        int hi = std::min(last_idx, b + main_lobe_radius);
        for (int k = lo; k <= hi; k++) {
            if (block_threshold[k] < claim_level) block_threshold[k] = claim_level;
        }
    }

    sort(kept.begin(), kept.end()); // restore ascending bin order
    peaks.swap(kept);
}

/**
 *more precsise information rather than what is straight up given to us from
 *the detected peaks*
 */
void parabolic_interpolation(const vector<double>& mag_spec, const vector<int>& peak_bins, vector<double>& true_freqs, vector<double>& true_mags) {
    true_freqs.clear();
    true_mags.clear();
    double alpha = 0.0;
    double beta = 0.0;
    double gamma = 0.0;
    double denom = 0.0;
    double p = 0.0;
    double true_bin = 0.0;
    double true_mag = 0.0;
    for (int bin : peak_bins) {
        alpha = mag_spec[bin-1];
        beta = mag_spec[bin];
        gamma = mag_spec[bin+1];
        denom = alpha - 2 * beta + gamma;
        p = 0.0;
        if (denom != 0.0) {
            p = 0.5 * (alpha - gamma) / denom;
        }
        true_bin = bin + p;
        true_mag = beta - 0.25 * (alpha-gamma) * p;
        true_freqs.push_back(true_bin);
        true_mags.push_back(true_mag);
    }
}

/**
 This is really important for our updating process. This checks to see what is
 inside of active peaks, and from this, if soemthing is within the threshold of
 frequency movment in this case 2 then bascially call it the same peak and
 change the contents at the existing peak. This is to prvent things liek having
 91hz, 92hz, 91.8hz, etc. Now we have 91hz, updated with 92 info, updated with 91.8 info
 */
int find_best_match_peak_HZ(double freq_hz, const vector<PeakTrack>& active_peaks,
                         double freq_tolerance_hz = 25.0)
{
    int best_idx = -1;
    double best_diff = numeric_limits<double>::infinity();
    for (size_t i = 0; i < active_peaks.size(); i++) {
        if (active_peaks[i].alive) {
            double diff = abs(active_peaks[i].freq_hz - freq_hz);
            if (diff < freq_tolerance_hz && diff < best_diff) {
                best_diff = diff;
                best_idx = (int)i;
            }
        }
    }
    return best_idx;
}

/**
 Mutual-exclusive nearest peak<->track assignment. The per-peak greedy match
 above lets two detected peaks within the tolerance both claim the same track:
 the second overwrites the first's update and the loser is never born as its
 own track. On close-spaced bass partials (SaintSaens pedal: 59.7 Hz sits
 20 Hz from the much louder 39.6 Hz) the weaker REAL partial was effectively
 deleted (-6 dB at unity, -26 dB under shift). Assign candidate pairs globally
 by ascending |df|; each peak and each track is used at most once, so a peak
 always keeps its nearest free track and the rest are born as new tracks.
 */
vector<int> assign_peaks_to_tracks(const vector<double>& freqs_hz,
                                   const vector<PeakTrack>& active_peaks,
                                   double freq_tolerance_hz = 25.0)
{
    vector<int> peak_to_track(freqs_hz.size(), -1);
    struct Cand { double df; int pi; int ti; };
    vector<Cand> cands;
    for (size_t i = 0; i < freqs_hz.size(); i++)
        for (size_t t = 0; t < active_peaks.size(); t++)
            if (active_peaks[t].alive) {
                double df = fabs(active_peaks[t].freq_hz - freqs_hz[i]);
                if (df < freq_tolerance_hz) cands.push_back({df, (int)i, (int)t});
            }
    sort(cands.begin(), cands.end(),
         [](const Cand& a, const Cand& b) { return a.df < b.df; });
    vector<char> track_used(active_peaks.size(), 0);
    for (auto& c : cands) {
        if (peak_to_track[c.pi] != -1 || track_used[c.ti]) continue;
        peak_to_track[c.pi] = c.ti;
        track_used[c.ti] = 1;
    }
    return peak_to_track;
}

/**
 So the two functions below are specific to window switching. Single parabolic
 interpolation changes what we have in the frequency and magnitude based on the
 fact that we are having to scale a bin now to be more accurate to what the smaller
 window is going to be holding
 */
void single_parabolic_interpolation(const vector<double>& mag_spec, double bin, double& true_freq, double& true_mag) {
    int lower_bin = floor(bin);
    if (lower_bin == bin) {
        true_freq = -1.0;
        true_mag = mag_spec[bin];
        return;
    }
    if (lower_bin < 1 || lower_bin >= mag_spec.size()-1) {
        true_freq = bin;
        true_mag = mag_spec[round(bin)];
        return;
    }
    
    double alpha = mag_spec[lower_bin-1];
    double beta = mag_spec[lower_bin];
    double gamma = mag_spec[lower_bin+1];
    double denom = alpha - 2 * beta + gamma;
    double p = 0.0;
    if (denom != 0.0) {
        p = 0.5 * (alpha-gamma) / denom;
    }
    true_freq = lower_bin + p;
    true_mag = beta - 0.25 * (alpha-gamma) * p;
}


/**
 Same thing this is to make sure the phase is updated accordingly with the freqyncy
 to make sure that the scaling of window switches is done correctly.
 */
double interpolate_phase(const vector<double>& phase_spec, double bin) {
    int lower_bin = floor(bin);
    int upper_bin = ceil(bin);
    if (bin == lower_bin && bin == upper_bin) {
        return phase_spec[bin];
    }
    
    if (lower_bin < 0 || upper_bin >= phase_spec.size()) {
        return phase_spec[round(bin)];
    }
    
    double lower_phase = phase_spec[lower_bin];
    double upper_phase = phase_spec[upper_bin];
    
    double phase_diff = upper_phase - lower_phase;
    phase_diff = fmod(phase_diff + M_PI, 2*M_PI) - M_PI;
    
    double fraction = bin - lower_bin;
    double interpolated_phase = lower_phase + fraction * phase_diff;
    return fmod(interpolated_phase + M_PI, 2*M_PI) - M_PI;
}

/**
 This is how we get the transients. Just walking through the file and comparing with contents
 in seperate frames to determine energy difference.
 */
#define NoTransients 0

vector<float> transientNegotiationTactics(int num_frames, float transientThresholdDB, int hop_size, int frame_size, vector<float>&singleChannelData) {
    TransientDetector mTD;
    
    mTD.mCurrentFrame.resize(frame_size);
    mTD.mPrevFrame.resize(frame_size);
    vector<float> transientList(num_frames, 0.0f);
    int halfFFTSize = frame_size / 2;
    vector<float> frame_data(frame_size, 0.0f);
    
    for (int f = 0; f < num_frames; f++) {
        int start = f  * hop_size;
        for (int i = 0; i < frame_size; i++) {
            frame_data[i] = singleChannelData[start + i] * long_hanning[i];
        }
        
        RealFFT(frame_data.data(), frame_size);
        memcpy(mTD.mPrevFrame.data(), mTD.mCurrentFrame.data(), frame_size*sizeof(float));
        memcpy(mTD.mCurrentFrame.data(), frame_data.data(), frame_size*sizeof(float));
        // Find the difference between the bins, add it up, see if the added sum is above a threshold,
        // and mark the transientList as either 0 or 1 right over the top of the sum
        MagnitudeFFTVec(mTD.mCurrentFrame);
        for (int j=1; j<halfFFTSize; j++)
        {
            const float eps = 1.0e-12f;
            float cur = max(mTD.mCurrentFrame[j], eps);
            float prev = max(mTD.mPrevFrame[j], eps);
            float diff = 20.0f * (log10(cur) - log10(prev));
            if (diff >= 0.0f)
                transientList[f] += diff;
        }
        transientList[f] /= halfFFTSize;
        if (transientList[f] > transientThresholdDB && f > 0){
            // Require a minimum gap of 2 frames between transients.
            // The old rule (suppress if previous == 1) missed every other hit in
            // rapid drum patterns.  Now we only suppress if either of the two
            // preceding frames was already marked as a transient.
            bool too_close = false;
            for (int back = 1; back <= 2 && (f - back) >= 0; back++) {
                if (transientList[f - back] == 1.0f) { too_close = true; break; }
            }
            transientList[f] = too_close ? 0.0f : 1.0f;
        }
        else{
            transientList[f] = 0.0f;
        }
#if NoTransients > 0
        transientList[f] = 0.0f;
#endif
    }
    return transientList;
}





int main(int argc, const char * argv[]) {
    int sr = 48000;
    int frame_size = LONG_SIZE;
    int hop_size = LONG_SIZE/2;
    
    /**
     *User settings, these are the things I think users should be able to control.
     *Right now how much noise, and how sensitive they want tranisents
     *and then semitone shifts
     */
    double thresholdMultiplier = 0.00025;
    float transientThresholdDB = 7.0f;
    int pitch_shift_semi = 0;


    double peak_sidelobe_attenuation_dB = 12.0;
    double peak_floor_below_max_dB = 60.0;
    // ===== Phase 1: amplitude / HF fidelity tuning =====
    // Asymmetric amplitude smoothing (replaces the symmetric 0.7/0.3 EMA): react
    // fast to rising partials (preserves attacks + upper harmonics), smooth slower
    // on decay. amp_smooth_* = weight on the newly measured value.
    // Set both to 0.7 to recover the old symmetric behavior.
    double amp_smooth_attack  = 0.90;
    double amp_smooth_release = 0.50;
    // Faster release when pitch-shifting. The slow release (plus the 85 ms
    // analysis-window average) stretches every decay; at unity that is just
    // slightly longer reverb at the SAME pitch (benign), but under shift the
    // stretched tail rings at the SHIFTED pitch after the note has ended —
    // the remaining "ghost" energy (+2 dB in the Female gaps). Analysis-side
    // but keyed to pitch_shift_semi, so unity output is byte-identical.
    double amp_smooth_release_shift = 0.85;
    // Slope-gated fast release (all modes): when a track's new measurement
    // drops below amp_release_fast_drop x its level (-6 dB in one frame),
    // that is a real event decay, not analysis ripple — track it fast. The
    // slow release blurring genuine drops is why the DrumLoop's second
    // cymbal crash read as "the first one never died" (the inter-crash dip
    // sat +3..4 dB above the input). Gentle decays (piano, ~-0.5 dB/frame)
    // never trigger the gate and keep the smooth release.
    double amp_release_fast = 0.85;
    double amp_release_fast_drop = 0.5;
    // Extra HF detection sensitivity (dB), applied to both the detection threshold
    // and the per-frame floor, ramped in over ~1.5–8 kHz. Recovers brightness lost
    // on rich tones (piano, Fairlight, choir). Set to 0 for old behavior.
    double hf_extra_sensitivity_db = 12.0;
    // ===== Phase 2: stochastic residual (noise / air / attack-sizzle fill) =====
    // Fills the per-bin magnitude deficit max(0,|X_in|-|X_model|) that the
    // sinusoidal model can't represent, with random-phase noise, above
    // residual_hp_hz. Magnitude-domain, so it never doubles well-modeled partials.
    // At unity the model is the output itself. Under pitch shift the deficit is
    // computed against a parallel UNSHIFTED tonal render and the noise is added
    // UNSHIFTED to the shifted output: reverb/breath/room physically do not
    // change pitch with the source, and shifting them was the Female "ghost
    // echo" (the gap reverb appeared as an exact xRatio copy of the input's).
    // Set residual_noise_gain = 0 to disable.
    double residual_noise_gain = 1.3;    // makeup gain on the noise fill
    double residual_hp_hz      = 2500.0; // unity: only fill above this frequency
    // Shift mode fills much lower: the ghost reverb components live at
    // 0.5-2 kHz (Female gaps: 505-1281 Hz shifted copies). Below this the
    // deficit is dominated by tonal-model error, not real noise.
    double residual_hp_hz_shift = 300.0;
    double residual_oversub    = 1.0;    // subtract this * model magnitude
    // peak_birth_confirm_frames: a new peak must be matched in this many
    //   consecutive frames before it's allowed to synthesize. 1 = old behavior
    //   (immediate). 2 = one frame of confirmation (kills single-frame
    //   phantoms — the most audible musical-noise source). 3+ is more
    //   aggressive but adds onset latency proportional to the long hop.
    int peak_birth_confirm_frames = 2;
    // ===== Shift-mode phase hygiene (warble fix) =====
    // Only active when pitch-shifting; unity output is byte-identical.
    // Under shift each track is an independent phase-propagated oscillator, so
    // duplicate detections of one partial beat against each other, and
    // frame-to-frame frequency-estimate jitter integrates into audible FM/AM
    // locked to the analysis frame rate (48000/4096 = 11.7 Hz).
    // shift_dedup_bins: per synthesis frame, drop a track whose freq is within
    //   this many 4096-FFT bins (11.72 Hz each) of a STRONGER track. 1.5 bins
    //   (~17.6 Hz) is inside one Hann main lobe -> can only be the same partial.
    //   0 disables.
    double shift_dedup_bins = 1.5;
    // shift_dedup_max_hz: only dedupe below this frequency. The audible beat
    //   pairs are low partials (Fairlight fundamentals 65-196 Hz, choir
    //   0.5-1.2 kHz); in dense mixes (take-me-out) partials from DIFFERENT
    //   instruments legitimately fall within one main lobe above ~2 kHz, and
    //   deduping them costs ~1 dB of highs (centroid 1787 -> 1615).
    double shift_dedup_max_hz = 1500.0;
    // shift_freq_smooth_alpha: weight on the new per-frame freq measurement in
    //   the synthesis-side per-track EMA (1.0 = no smoothing). Real partials are
    //   far more stable than per-frame estimates; smoothing removes estimate
    //   jitter before it is integrated into propagated phase.
    double shift_freq_smooth_alpha = 0.4;
    // shift_rebirth_max_gap_frames: a newborn track within max(8 Hz, 3%) of a
    //   track that died fewer than this many frames ago continues that track's
    //   propagated phase instead of re-seeding from analysis phase (which is
    //   meaningless at the shifted frequency). 0 disables.
    int shift_rebirth_max_gap_frames = 4;
    // ===== LF high-resolution analysis (bass rumble fix) =====
    // The 4096 Hann (11.72 Hz/bin, ~47 Hz main lobe) cannot resolve bass
    // partials spaced ~20 Hz (SaintSaens organ pedal 39.6/59.7/78.7 Hz): the
    // composite beating lobe spawns spurious 50-110 Hz tracks +18..+33 dB
    // above the input (the low rumble) and starves the real 59.7/118.3 Hz
    // partials. Below lf_cutoff_hz, peaks are taken from a LF_ANALYSIS_SIZE
    // (16384-sample, 2.93 Hz/bin) FFT centered at the same frame position
    // instead. (8192 was tried first: its ±11.7 Hz main lobe still let the
    // 39.6 Hz partial mask the 59.7 Hz one in 87% of frames.) 0 disables.
    double lf_cutoff_hz = 200.0;
    // Per-frame floor for LF peaks, dB below the loudest LF peak. Much tighter
    // than the global 60 dB floor: real bass partials sit within ~25 dB of the
    // strongest one, while the long LF window resolves piles of -40..-50 dB
    // junk whose per-frame flicker splatters broadband skirts (+9..+15 dB
    // between partials) and overshoots the limiter.
    double lf_floor_below_max_db = 35.0;
    // Kill a sub-cutoff track after this many consecutive unmatched (coasted)
    // frames — see PeakTrack::coast_count.
    int lf_coast_max_frames = 3;
    // LF stationarity gate: the LF splice only runs when the normalized LF
    // spectral flux between consecutive LF frames is below this. Sustained
    // bass (organ pedal ~0.05) passes; moving basslines / kick transients
    // fail and keep the plain 4096 analysis for that frame.
    double lf_max_flux = 0.10;
    // ===== Trajectory smoothness (Female "shaky voice" fix) =====
    // NOTE: these change UNITY output too (first knobs that do — approved).
    // traj_median_max_hz: apply the per-track median-of-3 freq filter below
    //   this frequency (vocal/instrument fundamentals + low partials).
    double traj_median_max_hz = 2000.0;
    // Rebirth credit: a newborn peak within max(8 Hz, 2%) of a track that
    //   died fewer than this many frames ago skips birth confirmation
    //   (confirmed immediately). A held vowel was carried by 4 track ids in
    //   34 frames; every handoff cost a 2-frame confirmation dropout.
    int rebirth_credit_max_gap_frames = 4;
    //End of user settings

    // effective EMA release for this run (see amp_smooth_release_shift)
    double amp_release_eff = (pitch_shift_semi != 0) ? amp_smooth_release_shift
                                                     : amp_smooth_release;
    
    
    AppSettings settings;
    
    
    /* parse input arguments */
    if (!parseArgs(argc, argv, settings))
    {
        printUsage();
        for(int j = 1; j < argc; j++)
            printf("%s\n", argv[j]);
        return -1;
    }
    
    AudioBuffer inputWav;
    inputWav.mSamples = nullptr;
    
    if (!readInWaveFile(settings.inputWavFilePath, &inputWav)) {
        cerr << "Error: Not read" << endl;
        return -1;
    }
    
    vector<vector<float>> audioData = audioBufferToVector(inputWav);
    vector<float> singleChannelData = audioData[0];
    int lengthYouNeed = inputWav.mNumSamples;
    cout << "Read " << inputWav.mNumSamples << " samples, " << inputWav.mChannels << " channels at " << inputWav.mSampleRate << " Hz.\n";
    
    int num_frames = (int)ceil((double)singleChannelData.size() / hop_size);
    singleChannelData.resize(singleChannelData.size() + frame_size, 0.0f);
    
    
    
    /**
     * Here is where we have the transients being thrown into a list. 1 for transient, 0 for nothing. Needed for window switiching
     */
    vector<float> transientList = transientNegotiationTactics(num_frames, transientThresholdDB, hop_size, LONG_SIZE, singleChannelData);

    /**
     *Calculating the FFT over frames keeping in mind that there needs to be window switching
     */
    STFTAdjustment stftAndSynthPlacement(LONG_SIZE, SHORT_SIZE, hop_size);
    stftAndSynthPlacement.processFrames(num_frames, singleChannelData, transientList, hop_size, long_hanning, short_hanning, long_to_short, short_to_long);
    vector<vector<complex<double>>> spec = stftAndSynthPlacement.getAdjustedSTFT();
    vector<SynthInformation> containsSynthPlacement = stftAndSynthPlacement.getSynthPlacement();
    const vector<vector<complex<double>>>& transientLongSpecs = stftAndSynthPlacement.getTransientLongSpecs();
    const unordered_map<int,int>& transientLongSpecForFrame = stftAndSynthPlacement.getTransientLongSpecForFrame();
    
    //Obtain the max magnitude this is important for threshold
    bool appliedShort = false;
    num_frames = spec.size();
    double max_value = 0.0;
    for (int f = 0; f < num_frames; f++) {
        frame_size = spec[f].size();
        for (int k = 0; k < frame_size/2; k++) {
            double mag = abs(spec[f][k]);
            if (mag > max_value) {
                max_value = mag;
            }
        }
    }
    
    double threshold = thresholdMultiplier * max_value;
    double threshold_long_analysis = threshold * ((double)ANALYSIS_SIZE / (double)LONG_SIZE);
    vector<PeakTrack> active_peaks;
    int peak_id_counter = 0;
    //We don't care about the signal once its below 70 db the highest tracked magnitude it obtained
    double threshold_factor = pow(10.0, (-70.0/20));
    
    vector<vector<PeakTrack>> frames_peaks;
    int printCoutner = 1;

    // Store previous frame's phase for instantaneous frequency calculation
    vector<double> prev_phase_spec;
    int prev_frame_size = 0;
    // Previous LF-FFT phase spectrum + its sample position, for the LF
    // phase-vocoder instantaneous frequency (parabolic alone on the LF FFT
    // drifts ~±1 Hz on bass partials the input holds to ±0.04 Hz).
    vector<double> prev_lf_phase;
    int prev_lf_start = -1;
    vector<double> prev_lf_mag_flux;   // last LF magnitude spectrum, for the stationarity gate
    // Recently-died tracks (freq, death frame) for rebirth credit at birth.
    vector<pair<double,int>> recent_deaths;

    for (int frame_idx = 0; frame_idx < num_frames; frame_idx++) {
        vector<complex<double>> frameGuy = spec[frame_idx];
        frame_size = frameGuy.size()*2;
        bool is_long_window = (frame_size == LONG_SIZE);
        vector<double> mag_spec(frame_size/2, 0.0);
        vector<double> phase_spec(frame_size/2, 0.0);


        vector<double> short_unmatched_long_mag;
        vector<double> short_unmatched_long_phase;
        int short_unmatched_long_frame_size = 0;

        for (int k = 0; k < frame_size/2; k++) {
            mag_spec[k] = abs(spec[frame_idx][k]);
            //arg calcualtes the phase angle
            phase_spec[k] = arg(spec[frame_idx][k]);
        }
        // Use a true 4096-sample analysis FFT for peak detection.
        // The synthesis still uses 2048-sample overlap-add windows (for transient
        // time resolution), but the analysis uses 4096 actual samples of audio
        // centered at the same position for true 11.72 Hz/bin frequency resolution.
        // This resolves low harmonics (e.g., 42 Hz vs 84 Hz are 3.6 bins apart
        // instead of 1.8 bins at 2048) without affecting the transient window switching.
        vector<double> analysis_mag;          // filled for long frames, used in unmatched update
        vector<double> analysis_phase_store;  // filled for long frames, stored as prev_phase_spec
        vector<double> lf_mag, lf_phase;      // LF FFT of this frame (long frames, lf_cutoff_hz > 0)
                                              // — the unmatched-coast source for sub-cutoff tracks
        int lf_start_frame = -1;              // sample position of this frame's LF window

        if (is_long_window) {
            // Center the 4096 analysis window at the same point as the 2048 frame center
            int frame_start = containsSynthPlacement[frame_idx].start;
            int frame_center = frame_start + frame_size / 2;
            int analysis_start = frame_center - ANALYSIS_SIZE / 2;

            vector<float> analysis_frame(ANALYSIS_SIZE, 0.0f);
            for (int i = 0; i < ANALYSIS_SIZE; i++) {
                int idx = analysis_start + i;
                if (idx >= 0 && idx < (int)singleChannelData.size())
                    analysis_frame[i] = singleChannelData[idx] * analysis_hanning[i];
            }
            RealFFT(analysis_frame.data(), ANALYSIS_SIZE);

            // Extract magnitude and phase from 4096-point FFT
            int a_half = ANALYSIS_SIZE / 2;
            analysis_mag.resize(a_half, 0.0);
            analysis_phase_store.resize(a_half, 0.0);
            analysis_mag[0] = fabs((double)analysis_frame[0]);
            analysis_phase_store[0] = (analysis_frame[0] >= 0) ? 0.0 : M_PI;
            for (int k = 1; k < a_half; k++) {
                double re = analysis_frame[k];
                double im = analysis_frame[ANALYSIS_SIZE - k];
                analysis_mag[k] = sqrt(re * re + im * im);
                analysis_phase_store[k] = atan2(im, re);
            }


            vector<int> peaks = detect_peaks(analysis_mag, threshold_long_analysis,
                                             sr, ANALYSIS_SIZE, hf_extra_sensitivity_db);
            filter_peaks_by_quality(analysis_mag, peaks,
                                    peak_sidelobe_attenuation_dB, peak_floor_below_max_dB,
                                    sr, ANALYSIS_SIZE, hf_extra_sensitivity_db);
            vector<double> freqs, mags;
            if (peaks.size() > 0) {
                parabolic_interpolation(analysis_mag, peaks, freqs, mags);

                vector<double> phases;
                for (int i : peaks)
                    phases.push_back(analysis_phase_store[i]);

                // Convert 4096-domain bins to Hz, apply instantaneous frequency
                vector<double> freqs_hz(freqs.size(), 0.0);
                for (size_t i = 0; i < freqs.size(); i++) {
                    double parabolic_freq_hz = freqs[i] * ((double)sr / (double)ANALYSIS_SIZE);

                    if (frame_idx > 0 && prev_frame_size == frame_size &&
                        (int)prev_phase_spec.size() == a_half) {
                        int peak_bin_a = peaks[i];
                        double current_phase = analysis_phase_store[peak_bin_a];
                        double previous_phase = prev_phase_spec[peak_bin_a];

                        double phase_diff = current_phase - previous_phase;
                        while (phase_diff > M_PI) phase_diff -= 2.0 * M_PI;
                        while (phase_diff < -M_PI) phase_diff += 2.0 * M_PI;

                        double expected = 2.0 * M_PI * freqs[i] * hop_size / (double)ANALYSIS_SIZE;

                        double deviation = phase_diff - expected;
                        while (deviation > M_PI) deviation -= 2.0 * M_PI;
                        while (deviation < -M_PI) deviation += 2.0 * M_PI;

                        double correction = deviation * sr / (2.0 * M_PI * hop_size);
                        double inst_freq = parabolic_freq_hz + correction;

                        if (parabolic_freq_hz < 150.0) {
                            freqs_hz[i] = 0.95 * inst_freq + 0.05 * parabolic_freq_hz;
                        } else if (parabolic_freq_hz < 500.0) {
                            freqs_hz[i] = 0.9 * inst_freq + 0.1 * parabolic_freq_hz;
                        } else if (parabolic_freq_hz < 2000.0) {
                            freqs_hz[i] = 0.5 * inst_freq + 0.5 * parabolic_freq_hz;
                        } else {
                            freqs_hz[i] = 0.2 * inst_freq + 0.8 * parabolic_freq_hz;
                        }
                    } else {
                        freqs_hz[i] = parabolic_freq_hz;
                    }
                }

                // ===== LF high-resolution replacement (see lf_cutoff_hz) =====
                // Below the cutoff, discard the 4096-derived peaks and re-detect
                // from an 8192 FFT centered at the same position. Appended
                // entries are converted to the 4096 magnitude scale and phase
                // reference, so tracking, the amplitude EMA and the synthesis
                // phase offset all work unchanged.
                if (lf_cutoff_hz > 0.0) {
                    int lf_start = frame_center - LF_ANALYSIS_SIZE / 2;
                    vector<float> lf_frame(LF_ANALYSIS_SIZE, 0.0f);
                    for (int i = 0; i < LF_ANALYSIS_SIZE; i++) {
                        int idx = lf_start + i;
                        if (idx >= 0 && idx < (int)singleChannelData.size())
                            lf_frame[i] = singleChannelData[idx] * lf_analysis_hanning[i];
                    }
                    RealFFT(lf_frame.data(), LF_ANALYSIS_SIZE);
                    int lf_half = LF_ANALYSIS_SIZE / 2;
                    lf_mag.assign(lf_half, 0.0);
                    lf_phase.assign(lf_half, 0.0);
                    for (int k = 1; k < lf_half; k++) {
                        double re = lf_frame[k];
                        double im = lf_frame[LF_ANALYSIS_SIZE - k];
                        lf_mag[k] = sqrt(re * re + im * im);
                        lf_phase[k] = atan2(im, re);
                    }

                    // LF stationarity gate. The 341 ms window only tells the
                    // truth when the bass is sustained (organ pedal); on a
                    // moving bassline it sees two consecutive notes at once,
                    // detects BOTH, and the doubled LF energy piles up until
                    // the limiter crushes the file (HappyMono -16 dB). Same
                    // philosophy as the long/short window switch, one tier up:
                    // only replace the 4096 peaks when the LF spectrum is
                    // stationary, measured as normalized spectral flux over
                    // the sub-cutoff bins between consecutive LF frames.
                    int flux_hi = (int)(lf_cutoff_hz * (double)LF_ANALYSIS_SIZE / sr) + 2;
                    if (flux_hi > lf_half) flux_hi = lf_half;
                    double lf_flux = 1e9; // no history -> treat as moving
                    if ((int)prev_lf_mag_flux.size() == lf_half) {
                        double num = 0.0, den = 1e-12;
                        for (int k = 1; k < flux_hi; k++) {
                            num += fabs(lf_mag[k] - prev_lf_mag_flux[k]);
                            den += prev_lf_mag_flux[k];
                        }
                        lf_flux = num / den;
                    }
                    prev_lf_mag_flux = lf_mag;
                    if (getenv("LF_DEBUG"))
                        fprintf(stderr, "F%d flux %.3f\n", frame_idx, lf_flux);
                    if (lf_flux > lf_max_flux) {
                        // moving bass: keep the 4096 peaks, disable LF coast
                        lf_mag.clear();
                        lf_phase.clear();
                    } else {

                    size_t w2 = 0;
                    for (size_t i = 0; i < freqs_hz.size(); i++) {
                        if (freqs_hz[i] >= lf_cutoff_hz) {
                            peaks[w2] = peaks[i]; freqs[w2] = freqs[i];
                            mags[w2] = mags[i]; phases[w2] = phases[i];
                            freqs_hz[w2] = freqs_hz[i]; w2++;
                        }
                    }
                    peaks.resize(w2); freqs.resize(w2); mags.resize(w2);
                    phases.resize(w2); freqs_hz.resize(w2);

                    double threshold_lf = threshold *
                        ((double)LF_ANALYSIS_SIZE / (double)LONG_SIZE);
                    vector<int> lf_peaks = detect_peaks(lf_mag, threshold_lf,
                                                        sr, LF_ANALYSIS_SIZE, 0.0);
                    {   // only below-cutoff candidates compete in the quality filter
                        vector<int> tmp;
                        for (int b : lf_peaks)
                            if (b * (double)sr / (double)LF_ANALYSIS_SIZE < lf_cutoff_hz)
                                tmp.push_back(b);
                        lf_peaks.swap(tmp);
                    }
                    filter_peaks_by_quality(lf_mag, lf_peaks,
                                            peak_sidelobe_attenuation_dB,
                                            lf_floor_below_max_db,
                                            sr, LF_ANALYSIS_SIZE, 0.0);
                    if (getenv("LF_DEBUG")) {
                        fprintf(stderr, "F%d:", frame_idx);
                        for (int b : lf_peaks)
                            fprintf(stderr, " %.1f/%.1f",
                                    b * (double)sr / LF_ANALYSIS_SIZE,
                                    20.0 * log10(lf_mag[b] + 1e-12));
                        fprintf(stderr, "\n");
                    }
                    if (!lf_peaks.empty()) {
                        vector<double> lf_freqs, lf_mags;
                        parabolic_interpolation(lf_mag, lf_peaks, lf_freqs, lf_mags);
                        double bin_to_4096 = (double)ANALYSIS_SIZE / (double)LF_ANALYSIS_SIZE;
                        for (size_t i = 0; i < lf_peaks.size(); i++) {
                            double f_hz = lf_freqs[i] * ((double)sr / (double)LF_ANALYSIS_SIZE);
                            if (f_hz <= 0.0 || f_hz >= lf_cutoff_hz) continue;
                            // Phase-vocoder instantaneous frequency, same as the
                            // 4096 path: parabolic alone drifts ~±1 Hz on bass
                            // partials the input holds to ±0.04 Hz.
                            int kb = lf_peaks[i];
                            if (prev_lf_start >= 0 &&
                                (int)prev_lf_phase.size() == lf_half) {
                                double dt = (double)(lf_start - prev_lf_start);
                                if (dt > 0) {
                                    double dphi = lf_phase[kb] - prev_lf_phase[kb];
                                    double expected = 2.0 * M_PI * kb * dt /
                                                      (double)LF_ANALYSIS_SIZE;
                                    double dev = dphi - expected;
                                    while (dev > M_PI) dev -= 2.0 * M_PI;
                                    while (dev < -M_PI) dev += 2.0 * M_PI;
                                    double inst = kb * ((double)sr / (double)LF_ANALYSIS_SIZE)
                                                + dev * sr / (2.0 * M_PI * dt);
                                    f_hz = 0.95 * inst + 0.05 * f_hz;
                                }
                            }
                            double omega = 2.0 * M_PI * f_hz / (double)sr;
                            // phase reference: the LF window starts (LF-4096)/2
                            // samples before the 4096 one (same frame center)
                            double ph = wrap_phase(lf_phase[kb] +
                                omega * (double)((LF_ANALYSIS_SIZE - ANALYSIS_SIZE) / 2));
                            peaks.push_back((int)round(f_hz * (double)ANALYSIS_SIZE / (double)sr));
                            freqs.push_back(f_hz * (double)ANALYSIS_SIZE / (double)sr);
                            mags.push_back(lf_mags[i] * bin_to_4096); // Hann mag scales with N
                            phases.push_back(ph);
                            freqs_hz.push_back(f_hz);
                        }
                    }
                    lf_start_frame = lf_start;
                    } // end stationarity-gated LF splice
                }

                vector<int> peak_to_track = assign_peaks_to_tracks(freqs_hz, active_peaks);
                for (size_t i = 0; i < peaks.size(); i++) {
                    double f_hz = freqs_hz[i];
                    double m_db = mags[i];
                    // Convert 4096-domain bin to 2048-domain for peak tracking
                    int p_bin = (int)round((double)peaks[i] / 2.0);
                    double ph = phases[i];

                    int match_idx = peak_to_track[i];
                    if (match_idx != -1) {
                        // Always mark as matched so the unmatched-update path
                        // (which uses the 2048 synthesis spectrum) never overwrites
                        // parameters that were measured from the 4096 analysis spectrum.
                        active_peaks[match_idx].edit = true;

                        // Trajectory de-jitter: median-of-3 over the RAW per-frame
                        // measurements, below traj_median_max_hz. On the Female
                        // held vowel the raw trajectory jumps 9.4 Hz/frame (real
                        // vibrato slope 5.8, max spike 66 Hz) — the rendered
                        // "shaky voice". A median kills single-frame spikes but
                        // passes monotone vibrato ramps EXACTLY (no depth loss,
                        // unlike an EMA).
                        double f_raw = f_hz;
                        {
                            PeakTrack &tk = active_peaks[match_idx];
                            if (f_hz < traj_median_max_hz &&
                                tk.freq_hist1 > 0.0 && tk.freq_hist2 > 0.0 &&
                                fabs(f_hz - tk.freq_hist1) < 0.10 * tk.freq_hist1) {
                                double a1 = f_hz, a2 = tk.freq_hist1, a3 = tk.freq_hist2;
                                double lo_ = fmin(a1, fmin(a2, a3));
                                double hi_ = fmax(a1, fmax(a2, a3));
                                f_hz = a1 + a2 + a3 - lo_ - hi_; // median
                            }
                            tk.freq_hist2 = tk.freq_hist1;
                            tk.freq_hist1 = f_raw;
                        }
                        // Conditional frequency smoothing for low-freq tracked peaks
                        if (f_hz < 200.0 && active_peaks[match_idx].freq_hz > 0.0) {
                            double freq_delta_pct = fabs(f_hz - active_peaks[match_idx].freq_hz) / active_peaks[match_idx].freq_hz;
                            if (freq_delta_pct < 0.05)
                                f_hz = 0.3 * f_hz + 0.7 * active_peaks[match_idx].freq_hz;
                        }


                        if (getenv("FTRAJ_DEBUG") && f_hz > 380.0 && f_hz < 480.0)
                            fprintf(stderr, "T %d %d %.3f %.4f %.5f\n",
                                    frame_idx, active_peaks[match_idx].id,
                                    f_hz, ph, m_db);
                        active_peaks[match_idx].freq_hz = f_hz;
                        active_peaks[match_idx].peak_bin = p_bin;
                        active_peaks[match_idx].phase = ph;
                        //light EMA (0.7 new / 0.3 old)
                        // tracks dynamics within a few frames but kills the
                        //as musical noise
                        { double _a = (m_db > active_peaks[match_idx].current_db) ? amp_smooth_attack : ((m_db < amp_release_fast_drop * active_peaks[match_idx].current_db) ? amp_release_fast : amp_release_eff); active_peaks[match_idx].current_db = _a * m_db + (1.0 - _a) * active_peaks[match_idx].current_db; }
                        active_peaks[match_idx].analysis_fft_size = ANALYSIS_SIZE;

                        if (active_peaks[match_idx].current_db > active_peaks[match_idx].max_db)
                            active_peaks[match_idx].max_db = active_peaks[match_idx].current_db;

                        double thresholdDB = threshold_factor * active_peaks[match_idx].max_db;
                        if (active_peaks[match_idx].current_db < thresholdDB)
                            active_peaks[match_idx].alive = false;

                        active_peaks[match_idx].coast_count = 0;
                        active_peaks[match_idx].matched_count++;
                        if (active_peaks[match_idx].matched_count >= peak_birth_confirm_frames) {
                            active_peaks[match_idx].confirmed = true;
                        }
                    } else {
                        PeakTrack newPeak(peak_id_counter++, f_hz, m_db, p_bin, ph);
                        if (peak_birth_confirm_frames <= 1) newPeak.confirmed = true;
                        // Rebirth credit: continuation of a just-died track at
                        // ~the same freq is not a phantom — confirm immediately
                        // so the handoff doesn't punch a 2-frame dropout into a
                        // sustained note.
                        if (!newPeak.confirmed)
                            for (auto &d : recent_deaths)
                                if (fabs(d.first - f_hz) <= fmax(8.0, 0.02 * d.first)) {
                                    newPeak.confirmed = true;
                                    break;
                                }
                        active_peaks.push_back(newPeak);
                    }
                }
            }
        } else {
            // For short frames at a transient: use the side long-window FFT that was
            // computed at the transient position.  This gives full frequency resolution
            // for peak detection instead of the coarse short-window spectrum.
            auto it = transientLongSpecForFrame.find(frame_idx);
            if (it != transientLongSpecForFrame.end()) {
                const vector<complex<double>>& long_spec = transientLongSpecs[it->second];
                short_unmatched_long_frame_size = (int)long_spec.size() * 2; // LONG_SIZE
                short_unmatched_long_mag.assign(long_spec.size(), 0.0);
                short_unmatched_long_phase.assign(long_spec.size(), 0.0);
                for (int k = 0; k < (int)long_spec.size(); k++) {
                    short_unmatched_long_mag[k]   = abs(long_spec[k]);
                    short_unmatched_long_phase[k] = arg(long_spec[k]);
                }
                // Local aliases for the rest of this block (matched-peak path)
                int long_frame_size = short_unmatched_long_frame_size;
                vector<double>& long_mag_spec   = short_unmatched_long_mag;
                vector<double>& long_phase_spec = short_unmatched_long_phase;


                // Transient long-specs are 2048-point, same scale as `spec`,
                vector<int> peaks = detect_peaks(long_mag_spec, threshold,
                                                 sr, long_frame_size, hf_extra_sensitivity_db);
                filter_peaks_by_quality(long_mag_spec, peaks,
                                        peak_sidelobe_attenuation_dB, peak_floor_below_max_dB,
                                        sr, long_frame_size, hf_extra_sensitivity_db);
                vector<double> freqs, mags;
                if (peaks.size() > 0) {
                    parabolic_interpolation(long_mag_spec, peaks, freqs, mags);
                    vector<double> s_freqs_hz(freqs.size(), 0.0);
                    for (size_t i = 0; i < freqs.size(); i++)
                        s_freqs_hz[i] = freqs[i] * ((double)sr / (double)long_frame_size);
                    vector<int> peak_to_track = assign_peaks_to_tracks(s_freqs_hz, active_peaks);
                    for (size_t i = 0; i < peaks.size(); i++) {
                        double f_hz  = s_freqs_hz[i];
                        double m_db  = mags[i];
                        int    p_bin = peaks[i];
                        double ph    = long_phase_spec[p_bin];

                        int match_idx = peak_to_track[i];
                        if (match_idx != -1) {
                            active_peaks[match_idx].freq_hz    = f_hz;
                            //temporal smoothing on current_db, same as long path.
                            { double _a = (m_db > active_peaks[match_idx].current_db) ? amp_smooth_attack : ((m_db < amp_release_fast_drop * active_peaks[match_idx].current_db) ? amp_release_fast : amp_release_eff); active_peaks[match_idx].current_db = _a * m_db + (1.0 - _a) * active_peaks[match_idx].current_db; }
                            active_peaks[match_idx].peak_bin   = p_bin;
                            active_peaks[match_idx].phase      = ph;
                            active_peaks[match_idx].edit       = true;
                            active_peaks[match_idx].analysis_fft_size = long_frame_size;

                            if (active_peaks[match_idx].current_db > active_peaks[match_idx].max_db)
                                active_peaks[match_idx].max_db = active_peaks[match_idx].current_db;
                            double thresholdDB = threshold_factor * active_peaks[match_idx].max_db;
                            if (active_peaks[match_idx].current_db < thresholdDB)
                                active_peaks[match_idx].alive = false;


                            active_peaks[match_idx].coast_count = 0;
                        active_peaks[match_idx].matched_count++;
                            if (active_peaks[match_idx].matched_count >= peak_birth_confirm_frames) {
                                active_peaks[match_idx].confirmed = true;
                            }

                        } else {
                            PeakTrack newPeak(peak_id_counter++, f_hz, m_db, p_bin, ph, long_frame_size);

                            if (peak_birth_confirm_frames <= 1) newPeak.confirmed = true;
                            if (!newPeak.confirmed)   // rebirth credit (see long path)
                                for (auto &d : recent_deaths)
                                    if (fabs(d.first - f_hz) <= fmax(8.0, 0.02 * d.first)) {
                                        newPeak.confirmed = true;
                                        break;
                                    }
                            active_peaks.push_back(newPeak);
                        }
                    }
                }
            }
        }

        for (auto &ap : active_peaks) {
            if (ap.alive && !ap.edit) {
                double scaled_bin;
                if (is_long_window && lf_cutoff_hz > 0.0 && !lf_mag.empty() &&
                    ap.freq_hz > 0.0 && ap.freq_hz < lf_cutoff_hz) {
                    // Sub-cutoff tracks must coast on the LF spectrum. On the
                    // 4096 spectrum a weak bass track (e.g. 94.5 Hz) sits on the
                    // huge unresolved composite lobe of its strong neighbours
                    // (78.7/118.3 Hz), so a single unmatched frame inflated it
                    // by tens of dB via the fast-attack EMA — inaudible at unity
                    // (analysis-locked phases re-sum the error) but a loud
                    // independent partial under pitch shift. Same lesson as the
                    // short-frame fix: every update path needs a magnitude
                    // source with resolution comparable to the matched path.
                    double lf_bin = ap.freq_hz * (double)LF_ANALYSIS_SIZE / (double)sr;
                    ap.coast_count++;
                    if (ap.coast_count > lf_coast_max_frames) {
                        ap.alive = false;
                    } else if (lf_bin >= 1.0 && lf_bin < (double)lf_mag.size() - 1.0) {
                        double true_freq, true_mag;
                        single_parabolic_interpolation(lf_mag, lf_bin, true_freq, true_mag);
                        // convert to the 4096 magnitude scale used track-wide
                        true_mag *= (double)ANALYSIS_SIZE / (double)LF_ANALYSIS_SIZE;
                        { double _a = (true_mag > ap.current_db) ? amp_smooth_attack : ((true_mag < amp_release_fast_drop * ap.current_db) ? amp_release_fast : amp_release_eff); ap.current_db = _a * true_mag + (1.0 - _a) * ap.current_db; }
                        double f_new = (true_freq >= 0.0)
                            ? true_freq * ((double)sr / (double)LF_ANALYSIS_SIZE) : ap.freq_hz;
                        double omega = 2.0 * M_PI * f_new / (double)sr;
                        ap.phase = wrap_phase(interpolate_phase(lf_phase, lf_bin) +
                            omega * (double)((LF_ANALYSIS_SIZE - ANALYSIS_SIZE) / 2));
                        ap.freq_hz = f_new;
                        ap.analysis_fft_size = ANALYSIS_SIZE;
                        if (ap.current_db > ap.max_db) ap.max_db = ap.current_db;
                        double thresholdDB = threshold_factor * ap.max_db;
                        if (ap.current_db < thresholdDB) ap.alive = false;
                    } else {
                        ap.alive = false;
                    }
                } else if (is_long_window && !analysis_mag.empty()) {
                    // Scale peak_bin from 2048 domain to 4096 domain
                    // (A coast cap like the LF one was tried here — it sharpened
                    // the DrumLoop inter-crash dip slightly more but cost dense
                    // mixes ~1.5 dB of real energy. The slope-gated release gets
                    // most of the contrast win without that cost.)
                    scaled_bin = ap.peak_bin * 2.0;
                    if (scaled_bin >= 0 && scaled_bin < (int)analysis_mag.size() - 1) {
                        double true_freq, true_mag;
                        single_parabolic_interpolation(analysis_mag, scaled_bin, true_freq, true_mag);
                        { double _a = (true_mag > ap.current_db) ? amp_smooth_attack : ((true_mag < amp_release_fast_drop * ap.current_db) ? amp_release_fast : amp_release_eff); ap.current_db = _a * true_mag + (1.0 - _a) * ap.current_db; }
                        ap.phase = interpolate_phase(analysis_phase_store, scaled_bin);
                        if (true_freq >= 0.0) {
                            ap.freq_hz = true_freq * ((double)sr / (double)ANALYSIS_SIZE);
                        }
                        ap.analysis_fft_size = ANALYSIS_SIZE;
                        if (ap.current_db > ap.max_db) ap.max_db = ap.current_db;
                        double thresholdDB = threshold_factor * ap.max_db;
                        if (ap.current_db < thresholdDB) ap.alive = false;
                    } else {
                        ap.alive = false;
                    }
                } else if (!is_long_window && !short_unmatched_long_mag.empty()) {
                    scaled_bin = ap.peak_bin; // already in 2048 domain
                    if (scaled_bin >= 0 && scaled_bin < (int)short_unmatched_long_mag.size() - 1) {
                        double true_freq, true_mag;
                        single_parabolic_interpolation(short_unmatched_long_mag, scaled_bin, true_freq, true_mag);

                        { double _a = (true_mag > ap.current_db) ? amp_smooth_attack : ((true_mag < amp_release_fast_drop * ap.current_db) ? amp_release_fast : amp_release_eff); ap.current_db = _a * true_mag + (1.0 - _a) * ap.current_db; }
                        ap.phase = interpolate_phase(short_unmatched_long_phase, scaled_bin);
                        if (true_freq >= 0.0) {
                            ap.freq_hz = true_freq * ((double)sr / (double)short_unmatched_long_frame_size);
                        }
                        ap.analysis_fft_size = short_unmatched_long_frame_size; // 2048
                        if (ap.current_db > ap.max_db) ap.max_db = ap.current_db;
                        double thresholdDB = threshold_factor * ap.max_db;
                        if (ap.current_db < thresholdDB) ap.alive = false;
                    } else {
                        ap.alive = false;
                    }
                } else {
                    // Legacy fallback: long-frame without analysis_mag (shouldn't
                    // happen) or short-frame outside any transient region (also
                    // shouldn't happen — short frames only exist at transients).
                    // Coast on the last known magnitude rather than reading from
                    // the coarse spectrum.
                    double thresholdDB = threshold_factor * ap.max_db;
                    if (ap.current_db < thresholdDB) ap.alive = false;
                }
            }
            ap.edit = false;
        }
        
        //Remove the unalive
        {
            vector<PeakTrack> temp;
            for (auto &p : active_peaks) {
                if (p.alive) temp.push_back(p);
                // Only CONFIRMED tracks earn rebirth credit — otherwise every
                // dying phantom lets the next noise peak skip confirmation and
                // the birdie defense is void.
                else if (p.confirmed) recent_deaths.push_back({p.freq_hz, frame_idx});
            }
            active_peaks.swap(temp);
            // prune expired death records
            vector<pair<double,int>> keep;
            for (auto &d : recent_deaths)
                if (frame_idx - d.second <= rebirth_credit_max_gap_frames)
                    keep.push_back(d);
            recent_deaths.swap(keep);
        }
        vector<PeakTrack> frame_info;
        for (auto &ap : active_peaks) {
            if (!ap.confirmed) continue;
            frame_info.push_back(ap);
        }
        printCoutner++;
        frames_peaks.push_back(frame_info);

        // Store phase spectrum for next frame's instantaneous frequency calculation.
        // For long frames, use the zero-padded phase spectrum so the instantaneous
        // frequency calculation in the next long frame compares at the finer resolution.
        if (is_long_window && !analysis_phase_store.empty()) {
            prev_phase_spec = std::move(analysis_phase_store);
        } else {
            prev_phase_spec = phase_spec;
        }
        prev_frame_size = frame_size;
        if (is_long_window && !lf_phase.empty() && lf_start_frame != -1) {
            prev_lf_phase = std::move(lf_phase);
            prev_lf_start = lf_start_frame;
        }

    }
    
    
    //==========================================================================
    // SYNTHESIS — Center-Relative Oscillator Bank with 4096 Analysis
    //==========================================================================
    frame_size = LONG_SIZE;
    float total_length = (float)lengthYouNeed;
    vector<float> synthesized_signal(total_length, 0.0f);
    vector<float> window_sum(total_length, 0.0f);
    unordered_map<int, double> synth_phase_by_track;
    unordered_set<int> synth_phase_initialized;
    
    vector<float> frame_signal(LONG_SIZE, 0.0f);
    vector<float> frame_signal_short(SHORT_SIZE, 0.0f);

    // Shift-mode residual support: a parallel UNSHIFTED tonal render (unity
    // phase rule, same windows/OLA/normalization). The stochastic-residual
    // block measures the noise the model can't represent against THIS, then
    // adds it unshifted to the shifted output — reverb/room must not shift.
    bool want_unity_model = (pitch_shift_semi != 0 && residual_noise_gain > 0.0);
    vector<float> synth_unity(want_unity_model ? (size_t)total_length : 0, 0.0f);
    vector<float> frame_unity(want_unity_model ? LONG_SIZE : 0, 0.0f);
    vector<float> frame_unity_short(want_unity_model ? SHORT_SIZE : 0, 0.0f);

    double nyquist = (double)sr / 2.0;

    vector<int> chordIntervals = {0};

    // Shift-mode: per-frame duplicate-track suppression. Two tracks closer than
    // ~1.5 analysis bins sit inside one Hann main lobe and are the same partial;
    // at unity their analysis-derived phases sum coherently (harmless), but under
    // shift they become independent oscillators that beat (the up-shift warble).
    // Done here on the synthesis input — not in analysis — so the weaker track
    // isn't sent to the unmatched-coast path where it would survive on the
    // stronger partial's leakage skirt. Compare amplitudes in the synthesis
    // domain (current_db / fft_size) since tracks can carry different fft sizes.
    if (pitch_shift_semi != 0 && shift_dedup_bins > 0.0) {
        double dedup_hz = shift_dedup_bins * (double)sr / (double)ANALYSIS_SIZE;
        // Below the LF cutoff the peaks come from the LF_ANALYSIS_SIZE FFT,
        // which genuinely resolves partials a few Hz apart (SaintSaens has
        // real ones 7-10 Hz apart at 150-167 Hz) — the dedupe radius must
        // match the resolution of the analysis that produced the track, or
        // it merges real neighbours.
        double dedup_hz_lf = shift_dedup_bins * (double)sr / (double)LF_ANALYSIS_SIZE;
        for (auto &fp : frames_peaks) {
            if (fp.size() < 2) continue;
            vector<char> drop(fp.size(), 0);
            for (size_t i = 0; i < fp.size(); i++) {
                if (drop[i]) continue;
                if (fp[i].freq_hz > shift_dedup_max_hz) continue;
                double amp_i = fp[i].current_db / (double)fp[i].analysis_fft_size;
                for (size_t j = i + 1; j < fp.size(); j++) {
                    if (drop[j]) continue;
                    if (fp[j].freq_hz > shift_dedup_max_hz) continue;
                    double radius = (lf_cutoff_hz > 0.0 &&
                                     fp[i].freq_hz < lf_cutoff_hz &&
                                     fp[j].freq_hz < lf_cutoff_hz) ? dedup_hz_lf : dedup_hz;
                    if (fabs(fp[i].freq_hz - fp[j].freq_hz) >= radius) continue;
                    double amp_j = fp[j].current_db / (double)fp[j].analysis_fft_size;
                    // Only drop a clearly weaker twin (>6 dB down). Near-equal
                    // neighbours are real distinct partials; a hard keep-stronger
                    // rule re-decides every frame as measured freqs wander across
                    // the radius, gating the loser on/off at frame rate (audible
                    // splatter). True duplicates (leakage skirts) sit 15-30 dB
                    // below their partial and are still removed.
                    double hi = (amp_j > amp_i) ? amp_j : amp_i;
                    double lo = (amp_j > amp_i) ? amp_i : amp_j;
                    if (lo > hi * 0.501 /* -6 dB */) continue;
                    if (amp_j > amp_i) { drop[i] = 1; break; }
                    drop[j] = 1;
                }
            }
            size_t w = 0;
            for (size_t i = 0; i < fp.size(); i++)
                if (!drop[i]) fp[w++] = fp[i];
            fp.erase(fp.begin() + w, fp.end());
        }
    }

    // Shift-mode per-track synthesis state (see settings comment):
    // smoothed freq trajectory + last propagated phase/position for rebirth.
    unordered_map<int, double> shift_freq_smoothed;
    struct ShiftTrackMemo { double freq_hz; double phase; long long pos; int frame_idx; };
    unordered_map<int, ShiftTrackMemo> shift_track_memo;


    for (int frame_idx = 0; frame_idx < num_frames; frame_idx++) {
        SynthInformation current_information = containsSynthPlacement[frame_idx];
        int frame_size = current_information.size;
        bool is_long = (frame_size == LONG_SIZE);
        int center = frame_size / 2;
        
        double t_center_abs = (double)(current_information.start + center) / (double)sr;
        
        fill(frame_signal.begin(), frame_signal.end(), 0.0f);
        fill(frame_signal_short.begin(), frame_signal_short.end(), 0.0f);
        if (want_unity_model) {
            fill(frame_unity.begin(), frame_unity.end(), 0.0f);
            fill(frame_unity_short.begin(), frame_unity_short.end(), 0.0f);
        }

        bool shorter = false;

        // Shift-mode: ids alive in THIS frame (rebirth candidates must be dead).
        unordered_set<int> cur_frame_ids;
        if (pitch_shift_semi != 0)
            for (auto &p : frames_peaks[frame_idx]) cur_frame_ids.insert(p.id);

        for (auto &peak : frames_peaks[frame_idx]) {
            double freq = peak.freq_hz;
            double mag = 4.0 * peak.current_db / (double)peak.analysis_fft_size;
            double phase = peak.phase;

            // Parallel unshifted tonal model for the shift-mode residual:
            // exactly the unity branch (raw freq, analysis-derived phase).
            if (want_unity_model) {
                double inc_u = 2.0 * M_PI * freq / (double)sr;
                int off_u = peak.analysis_fft_size / 2 - frame_size / 2;
                double ph0_u = wrap_phase(phase + inc_u * off_u);
                if (is_long) {
                    for (int n = 0; n < frame_size; n++)
                        frame_unity[n] += (float)(mag * cos(ph0_u + inc_u * n));
                } else {
                    for (int n = 0; n < frame_size; n++)
                        frame_unity_short[n] += (float)(mag * cos(ph0_u + inc_u * n));
                }
            }


            // Shift-mode: smooth the per-track freq trajectory before it is
            // integrated into propagated phase. The 10% guard keeps genuine
            // jumps (new partial segment on a reused id) unsmoothed.
            if (pitch_shift_semi != 0 && shift_freq_smooth_alpha < 1.0) {
                auto itf = shift_freq_smoothed.find(peak.id);
                if (itf != shift_freq_smoothed.end() &&
                    fabs(freq - itf->second) < 0.10 * itf->second) {
                    // Below the LF cutoff smooth much harder: those tracks come
                    // from the long LF FFT, so real freq movement is already
                    // slow, and per-frame estimate steps at 20-200 Hz turn into
                    // audible sideband splatter a few Hz from strong partials.
                    // Between the LF cutoff and traj_median_max_hz the analysis
                    // median has already de-jittered the trajectory — smoothing
                    // it AGAIN here stacks lag on lag and made Female down5
                    // jerk 6.5 -> 15.2; skip (a = 1 keeps the measurement).
                    double a = (lf_cutoff_hz > 0.0 && freq < lf_cutoff_hz)
                                   ? 0.15
                                   : (freq < traj_median_max_hz
                                          ? 1.0 : shift_freq_smooth_alpha);
                    freq = a * freq + (1.0 - a) * itf->second;
                }
                shift_freq_smoothed[peak.id] = freq;
            }

            for (int interval : chordIntervals) {
                double shiftFactor = pow(2.0, (interval + pitch_shift_semi) / 12.0);
                double shiftedFreq = freq * shiftFactor;

                if (shiftedFreq >= nyquist || shiftedFreq <= 0.0) continue;

                // R5: when pitch-shifting, fade partials out with a raised cosine
                // over 0.9*Nyquist -> Nyquist instead of the hard cull, so up-shifted
                // content near the top doesn't turn harsh/brittle as partials pop in
                // and out. Gated to shift mode -> unity-pitch output is unchanged.
                double mag_syn = mag;
                if (pitch_shift_semi != 0 || interval != 0) {
                    double aa_lo = 0.9 * nyquist;
                    if (shiftedFreq > aa_lo)
                        mag_syn = mag * 0.5 * (1.0 + cos(M_PI * (shiftedFreq - aa_lo) / (nyquist - aa_lo)));
                }

                double phase_inc = 2.0 * M_PI * shiftedFreq / (double)sr;

                // The analysis FFT is centered at the frame center, so its
                // phase reference (sample 0 of the FFT) is at
                //   frame_center - analysis_fft_size/2
                // Synthesis starts at frame_start = frame_center - frame_size/2.
                // Correct for this offset so the oscillator aligns with the
                // true signal phase at the synthesis start position.
                int phase_offset_samples = peak.analysis_fft_size / 2 - frame_size / 2;

                double phase0;
                if (pitch_shift_semi == 0 && interval == 0) {
                    // No pitch shift: derive phase directly from analysis each
                    // frame.  This locks the oscillator to the measured signal
                    // phase and prevents cumulative drift.  OLA windowing
                    // handles smooth blending between adjacent frames.
                    phase0 = wrap_phase(phase + phase_inc * phase_offset_samples);
                } else {
                    // Pitch shifting: propagate oscillator phase across frames
                    // so the shifted partial stays continuous at the new freq.
                    if (!synth_phase_initialized.count(peak.id)) {
                        // Rebirth continuity: if a track at ~this freq died a few
                        // frames ago, continue its propagated phase. Re-seeding
                        // from analysis phase is a random jump at the shifted
                        // freq that beats against the OLA overlap.
                        bool inherited = false;
                        if (shift_rebirth_max_gap_frames > 0) {
                            int best_id = -1;
                            double best_df = 1e18;
                            for (auto &kv : shift_track_memo) {
                                if (cur_frame_ids.count(kv.first)) continue;
                                const ShiftTrackMemo &m = kv.second;
                                if (frame_idx - m.frame_idx > shift_rebirth_max_gap_frames) continue;
                                double df = fabs(m.freq_hz - freq);
                                if (df <= max(8.0, 0.03 * m.freq_hz) && df < best_df) {
                                    best_df = df;
                                    best_id = kv.first;
                                }
                            }
                            if (best_id != -1) {
                                const ShiftTrackMemo &m = shift_track_memo[best_id];
                                double dead_inc = 2.0 * M_PI * (m.freq_hz * shiftFactor) / (double)sr;
                                long long gap = (long long)current_information.start - m.pos;
                                synth_phase_by_track[peak.id] =
                                    wrap_phase(m.phase + dead_inc * (double)gap);
                                shift_freq_smoothed[peak.id] = m.freq_hz;
                                shift_track_memo.erase(best_id);
                                inherited = true;
                            }
                        }
                        if (!inherited) {
                            synth_phase_by_track[peak.id] =
                                wrap_phase(phase + phase_inc * phase_offset_samples);
                        }
                        synth_phase_initialized.insert(peak.id);
                    }
                    phase0 = synth_phase_by_track[peak.id];
                }

                for (int n = 0; n < frame_size; n++) {
                    double sample_phase = phase0 + phase_inc * n;
                    if (is_long) {
                        frame_signal[n] += static_cast<float>(mag_syn * cos(sample_phase));
                    } else {
                        frame_signal_short[n] += static_cast<float>(mag_syn * cos(sample_phase));
                    }
                }

                // Advance propagated phase for pitch-shift mode
                if (pitch_shift_semi != 0 || interval != 0) {
                    int next_delta_samples = 0;
                    if (frame_idx + 1 < num_frames) {
                        next_delta_samples =
                            containsSynthPlacement[frame_idx + 1].start - current_information.start;
                    }
                    synth_phase_by_track[peak.id] =
                        wrap_phase(phase0 + phase_inc * next_delta_samples);
                    // Remember state for rebirth continuity (phase refers to the
                    // next frame's start position). Unshifted freq is stored.
                    shift_track_memo[peak.id] = {
                        freq, synth_phase_by_track[peak.id],
                        (long long)current_information.start + next_delta_samples,
                        frame_idx };
                }

                shorter = !is_long;
            }
        }
        
        // Synthesis window
        const vector<float>* ola_window_ptr = nullptr;
        if (frame_idx == 0 && !current_information.trans) {
            ola_window_ptr = is_long ? &rect_fade_to_hann : &rect_fade_to_hann_short;
        } else {
            ola_window_ptr = &current_information.windowApplied;
        }
        const vector<float>& ola_window = *ola_window_ptr;
        
        if (is_long) {
            for (int i = 0; i < frame_size; i++)
                frame_signal[i] *= ola_window[i];
        } else {
            for (int i = 0; i < frame_size; i++)
                frame_signal_short[i] *= ola_window[i];
        }
        if (want_unity_model) {
            if (is_long) {
                for (int i = 0; i < frame_size; i++)
                    frame_unity[i] *= ola_window[i];
            } else {
                for (int i = 0; i < frame_size; i++)
                    frame_unity_short[i] *= ola_window[i];
            }
        }

        // Overlap-add
        int start = current_information.start;
        int end = current_information.stop;
        if (end > (int)synthesized_signal.size()) {
            end = (int)synthesized_signal.size();
        }

        if (is_long) {
            for (int i = start; i < end; i++) {
                int local = i - start;
                synthesized_signal[i] += frame_signal[local];
                window_sum[i] += ola_window[local];
            }
        } else {
            for (int i = start; i < end; i++) {
                int local = i - start;
                synthesized_signal[i] += frame_signal_short[local];
                window_sum[i] += ola_window[local];
            }
        }
        if (want_unity_model) {
            const vector<float>& fu = is_long ? frame_unity : frame_unity_short;
            for (int i = start; i < end; i++)
                synth_unity[i] += fu[i - start];
        }
    }
    
    for (size_t i = 0; i < synthesized_signal.size(); i++) {
        if (window_sum[i] > 1.0e-8f) {
            synthesized_signal[i] /= window_sum[i];
            if (want_unity_model) synth_unity[i] /= window_sum[i];
        }
    }
    // Bridge overlap-add coverage notches at the short->long window switch. The
    // short_to_long transition window has a 448-sample zero prefix and the last
    // short frame's Hann tapers to zero at the junction, leaving ~2 samples with
    // window_sum ~ 0 (a COLA hole) -> a 2-sample dropout to zero = an impulsive
    // click on every transient. Linearly interpolate any run of near-zero-coverage
    // samples from the nearest well-covered neighbours (normal audio is untouched).
    {
        const float cov_thresh = 1.0e-3f;
        int N = (int)synthesized_signal.size();
        int i = 0;
        while (i < N) {
            if (window_sum[i] < cov_thresh) {
                int a = i - 1; int b = i;
                while (b < N && window_sum[b] < cov_thresh) b++;
                if (a >= 0 && b < N) {
                    float va = synthesized_signal[a], vb = synthesized_signal[b];
                    for (int k = i; k < b; k++) {
                        float t = (float)(k - a) / (float)(b - a);
                        synthesized_signal[k] = va + (vb - va) * t;
                    }
                    if (want_unity_model) {
                        float ua = synth_unity[a], ub = synth_unity[b];
                        for (int k = i; k < b; k++) {
                            float t = (float)(k - a) / (float)(b - a);
                            synth_unity[k] = ua + (ub - ua) * t;
                        }
                    }
                }
                i = b;
            } else {
                i++;
            }
        }
    }

    //==========================================================================
    // RESIDUAL MIX — blend original signal back in during transient regions
    //==========================================================================
    {
        vector<float> residual_mask(total_length, 0.0f);

        // Mark transient regions at frame resolution.
        // Each transient frame spans its synthesis placement window.
        for (int f = 0; f < (int)containsSynthPlacement.size(); f++) {
            // Walk the original transientList; every frame that was marked 1.0
            // should get residual.
            int frame_start = containsSynthPlacement[f].start;
            int orig_frame_idx = frame_start / hop_size;
            if (orig_frame_idx < 0) orig_frame_idx = 0;
            if (orig_frame_idx >= (int)transientList.size()) orig_frame_idx = (int)transientList.size() - 1;

            bool is_transient_region = (transientList[orig_frame_idx] == 1.0f);
            // Also check neighbors — the short frames around a transient need coverage
            if (!is_transient_region && orig_frame_idx > 0)
                is_transient_region = (transientList[orig_frame_idx - 1] == 1.0f);
            if (!is_transient_region && orig_frame_idx + 1 < (int)transientList.size())
                is_transient_region = (transientList[orig_frame_idx + 1] == 1.0f);

            if (is_transient_region) {
                int start = containsSynthPlacement[f].start;
                int stop  = containsSynthPlacement[f].stop;
                if (stop > (int)total_length) stop = (int)total_length;
                for (int i = start; i < stop; i++) {
                    residual_mask[i] = 1.0f;
                }
            }
        }

        int crossfade_len = (int)(0.002 * sr);
        if (crossfade_len < 4) crossfade_len = 4;
        for (int i = 1; i < (int)total_length; i++) {
            if (residual_mask[i] == 1.0f && residual_mask[i-1] == 0.0f) {
                int fade_start = max(0, i - crossfade_len);
                for (int j = fade_start; j < i; j++) {
                    float t = (float)(j - fade_start) / (float)crossfade_len;
                    float gain = 0.5f * (1.0f - cosf((float)M_PI * t));
                    if (gain > residual_mask[j]) residual_mask[j] = gain;
                }
            }
            if (residual_mask[i] == 0.0f && residual_mask[i-1] == 1.0f) {
                int fade_end = min((int)total_length, i + crossfade_len);
                for (int j = i; j < fade_end; j++) {
                    float t = (float)(j - i) / (float)crossfade_len;
                    float gain = 0.5f * (1.0f + cosf((float)M_PI * t));
                    if (gain > residual_mask[j]) residual_mask[j] = gain;
                }
            }
        }

        //output = (1 - mask) * additive + mask * original
        // R6(a): only paste the UNSHIFTED original at unity pitch. Under a pitch
        // shift this stamps original-pitch transients/onsets over the shifted tone
        // (the original note bleeding through, especially a fifth away), so skip it.
        if (pitch_shift_semi == 0)
        for (int i = 0; i < (int)total_length; i++) {
            float m = residual_mask[i];
            if (m > 0.0f) {
                float orig = (i < (int)singleChannelData.size()) ? singleChannelData[i] : 0.0f;
                synthesized_signal[i] = (1.0f - m) * synthesized_signal[i] + m * orig;
            }
        }

        // Count transient samples for diagnostics
        int transient_samples = 0;
        for (int i = 0; i < (int)total_length; i++) {
            if (residual_mask[i] > 0.01f) transient_samples++;
        }
    }

    //==========================================================================
    // STOCHASTIC RESIDUAL — fill the per-bin magnitude deficit the sinusoidal
    // model cannot represent (broadband noise, air, cymbal shimmer, attack
    // sizzle, reverb tails) with random-phase noise. Computed in the magnitude
    // domain, so at well-modelled partials the deficit is ~0 and no tonal
    // doubling occurs.
    //   unity: deficit vs the current synthesized_signal (which already
    //     includes the transient original-blend) — behavior unchanged.
    //   shift: deficit vs the parallel UNSHIFTED tonal render, and the noise
    //     is added UNSHIFTED to the shifted output. Reverb/breath/room do not
    //     change pitch with the source; pitch-shifting them was the "ghost
    //     echo" (gap reverb appeared as an exact xRatio copy of the input's).
    //==========================================================================
    if (residual_noise_gain > 0.0 && (pitch_shift_semi == 0 || want_unity_model)) {
        const int RN = 1024;           // residual FFT size
        const int RH = RN / 4;         // 75% overlap for smooth noise OLA
        const int RL = (int)total_length;
        const vector<float>& res_model =
            (pitch_shift_semi == 0) ? synthesized_signal : synth_unity;
        double hp_hz = (pitch_shift_semi == 0) ? residual_hp_hz : residual_hp_hz_shift;
        vector<float> rwin(RN);
        for (int i = 0; i < RN; i++)
            rwin[i] = 0.5f * (1.0f - cos(2.0 * M_PI * i / (RN - 1)));
        int hp_bin = (int)ceil(hp_hz * RN / (double)sr);

        vector<float> res_signal(RL, 0.0f);
        vector<float> res_wsum(RL, 0.0f);
        vector<float> fin(RN), fsy(RN), fout(RN);
        vector<double> mag_in(RN/2 + 1), mag_sy(RN/2 + 1), rmag(RN/2 + 1), sm(RN/2 + 1);

        // small deterministic PRNG for reproducible white phase
        unsigned int rng = 2463534242u;
        auto frand = [&rng]() {
            rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
            return (rng >> 8) * (1.0 / 16777216.0); // [0,1)
        };

        int n_res_frames = (RL >= RN) ? ((RL - RN) / RH + 1) : 0;
        for (int fr = 0; fr < n_res_frames; fr++) {
            int start = fr * RH;
            for (int i = 0; i < RN; i++) {
                int idx = start + i;
                float xi = (idx < (int)singleChannelData.size()) ? singleChannelData[idx] : 0.0f;
                float xs = (idx < RL) ? res_model[idx] : 0.0f;
                fin[i] = xi * rwin[i];
                fsy[i] = xs * rwin[i];
            }
            RealFFT(fin.data(), RN);
            RealFFT(fsy.data(), RN);
            mag_in[0]    = fabs(fin[0]);      mag_sy[0]    = fabs(fsy[0]);
            mag_in[RN/2] = fabs(fin[RN/2]);   mag_sy[RN/2] = fabs(fsy[RN/2]);
            for (int k = 1; k < RN/2; k++) {
                mag_in[k] = hypot((double)fin[k], (double)fin[RN-k]);
                mag_sy[k] = hypot((double)fsy[k], (double)fsy[RN-k]);
            }
            for (int k = 0; k <= RN/2; k++) {
                double d = mag_in[k] - residual_oversub * mag_sy[k];
                rmag[k] = (d > 0.0 && k >= hp_bin) ? d : 0.0;
            }
            // frequency-smooth (5-bin) so isolated partial deficits become a
            // broadband noise floor rather than tonal spikes
            for (int k = 0; k <= RN/2; k++) {
                double s = 0.0; int c = 0;
                for (int j = -2; j <= 2; j++) {
                    int kk = k + j;
                    if (kk >= 0 && kk <= RN/2) { s += rmag[kk]; c++; }
                }
                sm[k] = s / c;
            }
            for (int i = 0; i < RN; i++) fout[i] = 0.0f;
            for (int k = 1; k < RN/2; k++) {
                double th = 2.0 * M_PI * frand();
                fout[k]    = (float)(sm[k] * cos(th));
                fout[RN-k] = (float)(sm[k] * sin(th));
            }
            fout[0] = 0.0f; fout[RN/2] = 0.0f;   // no DC / Nyquist noise
            InvRealFFT(fout.data(), RN);
            for (int i = 0; i < RN; i++) {
                int idx = start + i;
                if (idx >= 0 && idx < RL) {
                    res_signal[idx] += fout[i] * rwin[i];
                    res_wsum[idx]   += rwin[i] * rwin[i];
                }
            }
        }
        // Normalize with a window-sum FLOOR. At frame edges res_wsum -> 0 (Hann is
        // ~0 there); dividing by it would create enormous spikes that slam the
        // output limiter and crush the whole file. Flooring at 0.3*peak bounds the
        // amplification and simply fades the residual out at the very edges.
        float wmax = 0.0f;
        for (int i = 0; i < RL; i++) if (res_wsum[i] > wmax) wmax = res_wsum[i];
        float wfloor = 0.3f * wmax;
        if (wfloor > 0.0f) {
            for (int i = 0; i < RL; i++) {
                float denom = res_wsum[i] > wfloor ? res_wsum[i] : wfloor;
                synthesized_signal[i] += residual_noise_gain * res_signal[i] / denom;
            }
        }
    }

    // Output limiting — windowed limiter. The old whole-file 0.95/max scale
    // meant ONE hot sample anywhere cost the entire file its level: up-shifted
    // renders (whose propagated phases occasionally align constructively) lost
    // a uniform 2-4 dB. Instead, dip only the few ms around each overshoot:
    // per-sample required gain, a 5 ms moving minimum (erosion doubles as
    // lookahead so the dip is fully in place at the peak), then a moving
    // average no wider than the erosion plateau, which smooths the gain curve
    // without lifting it at the peaks — so the ceiling still holds exactly.
    float max_val = 0.0f;
    for (auto val : synthesized_signal) {
        float abs_val = fabs(val);
        if (abs_val > max_val) max_val = abs_val;
    }
    if (max_val > 1.0f) {
        const float ceiling = 0.95f;
        int N = (int)synthesized_signal.size();
        int A = sr * 5 / 1000; // 5 ms erosion radius
        vector<float> g(N, 1.0f);
        for (int i = 0; i < N; i++) {
            float a = fabs(synthesized_signal[i]);
            if (a > ceiling) g[i] = ceiling / a;
        }
        // sliding-window minimum over [i-A, i+A] (monotonic deque)
        vector<float> gmin(N, 1.0f);
        {
            vector<int> dq(N + 2 * A + 2);
            int head = 0, tail = 0;
            for (int i = 0; i < N + A; i++) {
                if (i < N) {
                    while (tail > head && g[dq[tail - 1]] >= g[i]) tail--;
                    dq[tail++] = i;
                }
                int lo = i - 2 * A;
                while (tail > head && dq[head] < lo) head++;
                int out = i - A;
                if (out >= 0 && out < N) gmin[out] = g[dq[head]];
            }
        }
        // centered moving average, width 2A+1 (== the erosion plateau width)
        {
            double acc = 0.0;
            int W = 2 * A + 1;
            vector<float> gs(N, 1.0f);
            for (int i = 0; i < N + A; i++) {
                if (i < N) acc += gmin[i];
                if (i - W >= 0) acc -= gmin[i - W];
                int c = i - A;
                if (c >= 0 && c < N) {
                    int n_in = min(i, N - 1) - max(0, i - W + 1) + 1;
                    gs[c] = (float)(acc / n_in);
                }
            }
            for (int i = 0; i < N; i++)
                synthesized_signal[i] *= gs[i];
        }
        // safety: erosion+average guarantees the ceiling analytically; the
        // clamp only catches float rounding
        for (auto &val : synthesized_signal) {
            if (val > 0.999f) val = 0.999f;
            if (val < -0.999f) val = -0.999f;
        }
    }
     for (int i = 0; i < synthesized_signal.size(); i++) {
         if (synthesized_signal[i] >= 1 || synthesized_signal[i] <= -1) {
             cout << synthesized_signal[i] << " and " << i << endl;
         }
     }
    
    
    float** outputBuffer = new float*[1];
    outputBuffer[0] = new float[synthesized_signal.size()];
    for (size_t i = 0; i < synthesized_signal.size(); i++) {
        outputBuffer[0][i] = synthesized_signal[i];
    }
    
    writePCM16WaveFile(settings.outputWavFilePath, outputBuffer, synthesized_signal.size(), 1, sr);
    delete [] outputBuffer[0];
    delete [] outputBuffer;
    cout << "Wrote out the syntehsied outout to: " << settings.outputWavFilePath << endl;
    
    return 0;
}

void printUsage()
{
    cout << "Usage: AdditiveSynthFreqMask <input wav> <output wav> <block size in samples>" << endl;
    cout << "        <input wav>  Path to the input wav file." << endl;
    cout << "        <output_wav>  Path to the output wav file." << endl;
    cout << "        <block size> block size in samples." << endl;
}

bool parseArgs(int argc, const char* argv[], AppSettings& settings)
{
    if (argc != 4)
    {
        return false;
    }
    
    settings.inputWavFilePath = argv[1];
    //settings.inputCSVPath = argv[2];
    settings.outputWavFilePath = argv[2];
    //settings.sampleRate = atof(argv[4]);
    settings.blockSize = atof(argv[3]);
    
    return true;
}

bool readInWaveFile(const string& waveFile, AudioBuffer *buff)
{
    bool retVal = false;
    // Buffers etc..
    char ChunkID[4];
    char Format[4];
    char Subchunk1ID[4];
    char Subchunk2ID[4];
    
    int ChunkSize;
    int Subchunk1Size;
    int SampleRate;
    int ByteRate;
    int Subchunk2Size;
    
    long NumSamples;
    short AudioFormat;
    short NumChannels;
    short BlockAlign;
    short BitsPerSample;
    short *Data = NULL;
    
    FILE *fhandle= NULL;
    
    // Read the wave file
    fhandle = fopen(waveFile.c_str(), "rb");
    if (fhandle == NULL)
    {
        cerr << "\tWARNING: Failed to open wave file - " << waveFile << endl;
        return retVal;
    }
    
    fread(ChunkID,1,4,fhandle);
    fread(&ChunkSize,4,1,fhandle);
    fread(Format,1,4,fhandle);
    fread(Subchunk1ID,1,4,fhandle);
    fread(&Subchunk1Size,4,1,fhandle);
    fread(&AudioFormat,2,1,fhandle);
    fread(&NumChannels,2,1,fhandle);
    fread(&SampleRate,4,1,fhandle);
    buff->mChannels = NumChannels;
    buff->mSampleRate = SampleRate;
    buff->mSamples=new float*[NumChannels];
    
    fread(&ByteRate,4,1,fhandle);
    fread(&BlockAlign,2,1,fhandle);
    fread(&BitsPerSample,2,1,fhandle);
    fread(&Subchunk2ID,1,4,fhandle);
    fread(&Subchunk2Size,4,1,fhandle);
    size_t numShorts = Subchunk2Size / (BitsPerSample / 8);
    Data = new short [numShorts]; // Create an element for every sample
    memset(Data, 0, numShorts * sizeof(short));
    
    NumSamples = (Subchunk2Size / (BitsPerSample / 8)) / NumChannels;
    buff->mNumSamples = NumSamples;
    size_t samplesRead = fread(Data, sizeof(short), Subchunk2Size / (BitsPerSample / 8), fhandle); // Reading raw audio data
    if (samplesRead != numShorts)
    {
        cerr << "\tWARNING: Failed to read all samples in wav file - " << waveFile << endl;
        return retVal;
    }
    
    for (int ch = 0; ch < NumChannels; ch++)
    {
        buff->mSamples[ch] = new float[NumSamples];
    }
    fclose(fhandle);
    
    CSampleNormalizer Norm = CSampleNormalizer(BitsPerSample, NumChannels);
    Norm.Normalize(buff->mSamples, Data, NumSamples);
    
    delete[] Data;
    
    retVal = true;
    return retVal;
}

void writePCM16WaveFile(const string& waveFilePath, float** samples, size_t numSamples, short numChannels, int sampleRate)
{
    const int Subchunk2Size = (int)numSamples * numChannels * 2;
    const int Subchunk1Size = 16;
    const int ChunkSize = Subchunk1Size + Subchunk2Size;
    const short AudioFormat = 1;
    const int ByteRate = sampleRate * numChannels * 2;
    const short BlockAlign = 4;
    const short BitsPerSample = 16;
    short *Data = NULL;
    FILE *fhandle= NULL;
    
    Data = new short[numChannels * numSamples];
    CSampleNormalizer Norm = CSampleNormalizer(BitsPerSample, numChannels);
    Norm.ConvertToShort(Data, samples, numSamples);
    
    // Write the processed file
    fhandle=fopen(waveFilePath.c_str(), "wb");
    if (fhandle == NULL)
    {
        cerr << "\tWARNING: Failed to open output file - " << waveFilePath << endl;
        delete[] Data;
        return;
    }
    
    fwrite("RIFF", 1, 4, fhandle);
    fwrite(&ChunkSize, 4, 1, fhandle);
    fwrite("WAVE", 1, 4, fhandle);
    fwrite("fmt ", 1, 4, fhandle);
    fwrite(&Subchunk1Size, 4, 1, fhandle);
    fwrite(&AudioFormat, 2, 1, fhandle);
    fwrite(&numChannels, 2, 1, fhandle);
    fwrite(&sampleRate, 4, 1, fhandle);
    fwrite(&ByteRate, 4, 1, fhandle);
    fwrite(&BlockAlign, 2, 1, fhandle);
    fwrite(&BitsPerSample, 2, 1, fhandle);
    fwrite("data", 1, 4, fhandle);
    fwrite(&Subchunk2Size, 4, 1, fhandle);
    fwrite(Data, 2, Subchunk2Size / 2, fhandle);
    fclose(fhandle);
    
    delete [] Data;
}




