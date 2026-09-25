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

// Spectral-reassignment companion windows (used only when phase_mode==1). Built
// once from the analysis Hann window:
//   analysis_th[n] = (n - centre) * h[n]   (time-weighted, in samples)
//   analysis_dh[n] = d/dn h[n]             (central difference)
// These let the analysis measure the phase/frequency of a partial that CHIRPS
// within the 85 ms window under vibrato — the stationary phase read (even after
// the 5ae358e fractional-bin fix) is biased for a moving partial, scrambling the
// harmonics' relative phases (the "phasey/double voice"). See vibrato_rig.py.
vector<float> analysis_th = [](){
    vector<float> v(ANALYSIS_SIZE, 0.0f);
    double c = (ANALYSIS_SIZE - 1) / 2.0;
    for (int n = 0; n < ANALYSIS_SIZE; n++)
        v[n] = (float)((n - c) * analysis_hanning[n]);
    return v;
}();
vector<float> analysis_dh = [](){
    vector<float> v(ANALYSIS_SIZE, 0.0f);
    for (int n = 0; n < ANALYSIS_SIZE; n++) {
        float prev = (n > 0) ? analysis_hanning[n - 1] : 0.0f;
        float next = (n < ANALYSIS_SIZE - 1) ? analysis_hanning[n + 1] : 0.0f;
        v[n] = 0.5f * (next - prev);
    }
    return v;
}();

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
    int pitchShiftSemi = 0; // optional CLI override; 0 = unity (default)
    int synthMode = 0;      // optional CLI override; 0 = OLA (default), 1 = MQ
    double residualScale = 1.0; // diagnostic: multiplies residual_noise_gain (0 = off)
    int unityDedup = 0;         // diagnostic: run duplicate-partial dedup at unity too
    int residualMode = 0;       // 0 = stochastic residual, 1 = true-phase (unity)
    double residualHpHz = -1.0; // diagnostic: override residual_hp_hz (<0 = keep setting)
    int ampMode = 0;            // 0 = parabolic peak amplitude, 1 = energy-integrated (vibrato)
    int phaseMode = 0;          // 0 = stationary phase (default), 1 = reassigned (chirp-aware)
    int pitchSync = 1;          // 1 = auto (default): self-engages on monophonic + modulated
                                // (singing) material, exact passthrough elsewhere. 0 = force off.
    int jointMode = 2;          // 2 = joint least-squares amp/phase at the FINAL track
                                // frequencies (default; ear-approved Sep 2026).
                                // 1 = joint solve inside analysis (a near no-op: the
                                //     tracker then moves the frequencies -- see
                                //     joint_amp_phase). 0 = legacy per-peak estimate.
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

// ===================== Oscillator-bank (MQ) synthesis =====================
// One continuous oscillator per track across its whole life, instead of the
// legacy path's per-frame constant-frequency render overlap-added at hop 1024.
// The legacy overlap sums two constant-frequency copies of a partial whose
// frequency moved between frames -> a comb (the saw shape ripple, the Female
// "double voice" under vibrato, the steady-tone wobble). A single continuous
// oscillator has no overlap to comb.
//
// A node is one (time, freq, amp, phase) measurement at a frame centre.
//   match_phase == true  (unity): cubic phase interpolation matching the
//     MEASURED phase at both node ends (McAulay-Quatieri, with the added-cycles
//     M* term) -> the waveform shape is pinned to the analysis, fixing the saw.
//   match_phase == false (shift): phase is PROPAGATED by integrating a linear
//     frequency ramp between nodes (measured phase is only the birth seed) ->
//     every partial stays a smooth continuous trajectory, killing the wobble.
struct MQNode { double t; double freq; double amp; double phase; };

static void mq_synthesize(const std::unordered_map<int, std::vector<MQNode>>& tracks,
                          std::vector<float>& out, bool match_phase, int sr,
                          double nyquist, int fade) {
    const double TWO_PI = 2.0 * M_PI;
    const int N = (int)out.size();
    for (const auto& kv : tracks) {
        std::vector<MQNode> nodes = kv.second;
        std::sort(nodes.begin(), nodes.end(),
                  [](const MQNode& a, const MQNode& b) { return a.t < b.t; });
        if (nodes.empty()) continue;

        // Birth: fade amp 0 -> first node over `fade` samples, meeting its phase.
        {
            const MQNode& n0 = nodes.front();
            if (n0.freq > 0.0 && n0.freq < nyquist) {
                double w0 = TWO_PI * n0.freq / sr;
                int t0 = (int)llround(n0.t);
                int s0 = std::max(0, t0 - fade);
                for (int t = s0; t < t0 && t < N; t++) {
                    double a = n0.amp * (double)(t - s0) / (double)std::max(1, t0 - s0);
                    out[t] += (float)(a * cos(n0.phase + w0 * (t - n0.t)));
                }
            }
        }

        double phi = nodes.front().phase; // running phase for the propagate path
        for (size_t k = 0; k + 1 < nodes.size(); k++) {
            const MQNode& A = nodes[k];
            const MQNode& B = nodes[k + 1];
            int tA = (int)llround(A.t), tB = (int)llround(B.t);
            double T = (double)(tB - tA);
            if (T <= 0.0 || A.freq <= 0.0 || B.freq <= 0.0) { phi = B.phase; continue; }
            double w0 = TWO_PI * A.freq / sr, w1 = TWO_PI * B.freq / sr;
            if (match_phase) {
                double M = round((A.phase + 0.5 * (w0 + w1) * T - B.phase) / TWO_PI);
                double x = B.phase + TWO_PI * M - A.phase - w0 * T;
                double y = w1 - w0;
                double a2 = 3.0 * x / (T * T) - y / T;
                double a3 = -2.0 * x / (T * T * T) + y / (T * T);
                for (int t = tA; t < tB; t++) {
                    if (t < 0 || t >= N) continue;
                    double dt = (double)(t - tA);
                    double theta = A.phase + w0 * dt + a2 * dt * dt + a3 * dt * dt * dt;
                    double a = A.amp + (B.amp - A.amp) * (dt / T);
                    out[t] += (float)(a * cos(theta));
                }
            } else {
                for (int t = tA; t < tB; t++) {
                    if (t < 0 || t >= N) continue;
                    double dt = (double)(t - tA);
                    double theta = phi + w0 * dt + (w1 - w0) * dt * dt / (2.0 * T);
                    double a = A.amp + (B.amp - A.amp) * (dt / T);
                    out[t] += (float)(a * cos(theta));
                }
                phi = wrap_phase(phi + 0.5 * (w0 + w1) * T);
            }
        }

        // Death: fade last node amp -> 0 over `fade` samples.
        {
            const MQNode& nL = nodes.back();
            if (nL.freq > 0.0 && nL.freq < nyquist) {
                double wL = TWO_PI * nL.freq / sr;
                double phiL = match_phase ? nL.phase : phi;
                int tL = (int)llround(nL.t);
                int e = std::min(N, tL + fade);
                for (int t = tL; t < e; t++) {
                    if (t < 0) continue;
                    double a = nL.amp * (1.0 - (double)(t - tL) / (double)std::max(1, e - tL));
                    out[t] += (float)(a * cos(phiL + wL * (t - tL)));
                }
            }
        }
    }
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

// shapeThresholdDB > 0 additionally flags SPECTRAL-CHANGE events: a frame whose
// spectrum changes shape sharply while its LEVEL stays flat. The positive-flux sum
// below is an onset detector -- it answers "did energy arrive" -- and a wavetable
// switch answers no: measured on Fairlight C2 the level is flat to +-1 dB across the
// event while the spectrum moves 6-10 dB per bin, so the detector fires once in three
// seconds and every switch is analysed with an 85 ms window and smeared. A listener
// identified that blind as "too smooth, like we sanded out the fine details"
// (docs 5b.3). The shape measure normalises both frames to unit energy first, so it
// is blind to level by construction and catches exactly what the flux sum cannot.
vector<float> transientNegotiationTactics(int num_frames, float transientThresholdDB, int hop_size, int frame_size, vector<float>&singleChannelData, float shapeThresholdDB = 0.0f) {
    TransientDetector mTD;
    
    mTD.mCurrentFrame.resize(frame_size);
    mTD.mPrevFrame.resize(frame_size);
    vector<float> transientList(num_frames, 0.0f);
    vector<float> shape_series(num_frames, 0.0f);
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

        // Level-invariant spectral shape change. Normalise each spectrum to unit
        // total magnitude, then take the mean absolute dB difference over the bins
        // that carry the signal (a floor keeps near-empty bins from dominating).
        float shape_change = 0.0f;
        if (shapeThresholdDB > 0.0f && f > 0) {
            double sc = 0.0, sp = 0.0;
            for (int j = 1; j < halfFFTSize; j++) { sc += mTD.mCurrentFrame[j]; sp += mTD.mPrevFrame[j]; }
            if (sc > 1e-12 && sp > 1e-12) {
                double peak = 0.0;
                for (int j = 1; j < halfFFTSize; j++)
                    peak = max(peak, max((double)mTD.mCurrentFrame[j] / sc, (double)mTD.mPrevFrame[j] / sp));
                const double floor_rel = peak * 1.0e-3;     // -60 dB below the peak bin
                double acc = 0.0; int cnt = 0;
                for (int j = 1; j < halfFFTSize; j++) {
                    double a = (double)mTD.mCurrentFrame[j] / sc;
                    double b = (double)mTD.mPrevFrame[j] / sp;
                    if (a < floor_rel && b < floor_rel) continue;
                    a = max(a, floor_rel); b = max(b, floor_rel);
                    acc += fabs(20.0 * log10(a / b)); cnt++;
                }
                if (cnt > 0) shape_change = (float)(acc / cnt);
            }
        }

        shape_series[f] = shape_change;
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

    // Adaptive pass for the shape events. An ABSOLUTE threshold cannot work here:
    // on a dense mix the spectrum changes shape constantly, so a fixed threshold
    // fires on nearly every frame (measured: take-me-out went from 0.3% to 83.6%
    // of samples served by the transient original-blend, i.e. the engine stopped
    // being an additive engine). What distinguishes a wavetable switch is that it
    // is a SPIKE against that file's own background of spectral movement. So flag
    // a frame only when its shape change stands out from the local median by
    // shapeThresholdDB, and cap how many frames may qualify.
    if (shapeThresholdDB > 0.0f && num_frames > 4) {
        const int HALF = 12;                 // ~0.5 s of context each side
        const int MAX_FRAC_PCT = 8;          // never mark more than 8% this way
        vector<pair<float,int>> cand;
        vector<float> win;
        for (int f = 1; f < num_frames; f++) {
            if (transientList[f] == 1.0f) continue;
            int a = max(1, f - HALF), b = min(num_frames, f + HALF + 1);
            win.assign(shape_series.begin() + a, shape_series.begin() + b);
            if (win.size() < 5) continue;
            nth_element(win.begin(), win.begin() + win.size() / 2, win.end());
            float med = win[win.size() / 2];
            if (shape_series[f] > med + shapeThresholdDB)
                cand.push_back({shape_series[f] - med, f});
        }
        sort(cand.begin(), cand.end(), [](const pair<float,int>& a, const pair<float,int>& b){
            return a.first > b.first;                      // strongest first
        });
        int budget = max(1, num_frames * MAX_FRAC_PCT / 100);
        for (auto &c : cand) {
            if (budget <= 0) break;
            int f = c.second;
            bool too_close = false;
            for (int back = 1; back <= 2; back++) {
                if (f - back >= 0 && transientList[f - back] == 1.0f) too_close = true;
                if (f + back < num_frames && transientList[f + back] == 1.0f) too_close = true;
            }
            if (too_close) continue;
            transientList[f] = 1.0f; budget--;
        }
    }
    return transientList;
}


// ===================== Joint amplitude/phase estimation =========================
// The analysis reads each partial's amplitude and phase INDEPENDENTLY from its own
// spectral peak. At a 4096-point Hann the main lobe is ~47 Hz wide, so in any dense
// spectrum neighbouring partials contaminate each other's peak height and phase. That
// information is entangled, and no better per-peak reader can recover it — which is why
// energy-integrated amplitude (amp_mode=1) and spectral reassignment (phaseMode) both
// failed: both are still per-peak estimators.
//
// Measured headroom (docs/engine-direction-2026-09.md §4b.3): holding the engine's OWN
// frequencies and re-solving only amplitude and phase, INDEPENDENT least squares gains
// ~nothing (-6.5..+9.5 dB) while JOINT least squares gains +11..+22 dB on every class.
//
// So: model the frame as a sum of sinusoids at the already-detected frequencies,
//     x[n] ~= sum_i a_i*cos(w_i n) + b_i*sin(w_i n),
// and solve for (a_i, b_i) together in the Hann-weighted least-squares sense. The normal
// matrix entries are closed-form in the window's DTFT, the system is banded (partials
// only interact within a few bins), and it is solved by conjugate gradient starting from
// the existing independent estimate — so it strictly refines today's answer.
//
// Regularisation is centred on the independent estimate x0 rather than on zero:
//     minimise ||Ax - y||^2_w + lambda*||x - x0||^2
// so lambda -> infinity reproduces the current engine exactly and lambda -> 0 is the full
// joint solve. That bounds the worst case and gives a single safety knob.

// sum_{n=0}^{N-1} e^{j w n}
static inline void dirichlet_sum(double w, int N, double &re, double &im) {
    double s = sin(w * 0.5);
    double c = cos(w * (N - 1) * 0.5), sn = sin(w * (N - 1) * 0.5);
    double mag = (fabs(s) < 1e-12) ? (double)N : (sin(w * N * 0.5) / s);
    re = mag * c; im = mag * sn;
}

// Wc(w) = sum_n h[n] cos(w n), Ws(w) = sum_n h[n] sin(w n), for the SYMMETRIC Hann
// h[n] = 0.5 - 0.5*cos(2*pi*n/(N-1)) built by NormalWindows::HanningWindow.
// Verified against brute-force sums to ~1e-12 relative.
static inline void hann_dtft(double w, int N, double &Wc, double &Ws) {
    double wm = 2.0 * M_PI / (double)(N - 1);
    double r0, i0, rp, ip, rm, im_;
    dirichlet_sum(w, N, r0, i0);
    dirichlet_sum(w + wm, N, rp, ip);
    dirichlet_sum(w - wm, N, rm, im_);
    Wc = 0.5 * r0 - 0.25 * rp - 0.25 * rm;
    Ws = 0.5 * i0 - 0.25 * ip - 0.25 * im_;
}

// Overwrites mags[] (FFT-magnitude scale, i.e. A*N/4 to match the synthesis inverse
// 4*current_db/N) and phases[] (radians at n=0 of the analysis window, the existing
// convention) with the joint solution. xw = the Hann-WINDOWED analysis frame.
// slopes_hz_s (optional, may be null or all-zero) turns each atom into a linear-FM
// chirp: phase(n) = w*n + pi*(slope/sr^2)*(n - N/2)^2, i.e. the instantaneous frequency
// passes through freqs_hz[i] at the window CENTRE. The closed-form Hann-DTFT Gram is
// only valid for stationary atoms, so any pair involving a chirped atom is integrated
// numerically against `win` instead; stationary-stationary pairs keep the closed form
// and material without vibrato pays nothing. See chirp_mode.
static void joint_amp_phase(const vector<double>& xw, int N, int sr,
                            const vector<double>& freqs_hz,
                            vector<double>& mags, vector<double>& phases,
                            double band_bins, double reg, int max_iters,
                            int max_neighbors = 48,
                            const vector<double>* slopes_hz_s = nullptr,
                            const vector<float>* win = nullptr) {
    const int K = (int)freqs_hz.size();
    if (K == 0) return;

    vector<double> w(K);
    for (int i = 0; i < K; i++) w[i] = 2.0 * M_PI * freqs_hz[i] / (double)sr;

    // Quadratic phase coefficient per atom, 0 when stationary.
    vector<double> q(K, 0.0);
    bool any_chirp = false;
    if (slopes_hz_s && win && (int)slopes_hz_s->size() == K && (int)win->size() >= N) {
        for (int i = 0; i < K; i++) {
            double sl = (*slopes_hz_s)[i];
            if (sl != 0.0) { q[i] = M_PI * sl / ((double)sr * (double)sr); any_chirp = true; }
        }
    }
    const double half = 0.5 * (double)N;
    // Windowed atom samples, built once and reused by every numeric inner product.
    // cos/sin of the full phase INCLUDING the chirp, times the analysis window.
    vector<vector<double>> atom_c, atom_s;
    if (any_chirp) {
        atom_c.assign(K, {}); atom_s.assign(K, {});
        for (int i = 0; i < K; i++) {
            if (q[i] == 0.0) continue;                       // stationary: closed form
            atom_c[i].resize(N); atom_s[i].resize(N);
            for (int n = 0; n < N; n++) {
                double d = (double)n - half;
                double ph = w[i] * n + q[i] * d * d;
                double ww = (double)(*win)[n];
                atom_c[i][n] = cos(ph) * ww;
                atom_s[i][n] = sin(ph) * ww;
            }
        }
        // Stationary atoms still need samples when paired with a chirped one.
        for (int i = 0; i < K; i++) {
            if (q[i] != 0.0) continue;
            bool needed = false;
            for (int j = 0; j < K && !needed; j++)
                if (q[j] != 0.0 && fabs(freqs_hz[j] - freqs_hz[i]) <=
                        band_bins * (double)sr / (double)N) needed = true;
            if (!needed) continue;
            atom_c[i].resize(N); atom_s[i].resize(N);
            for (int n = 0; n < N; n++) {
                double ph = w[i] * n, ww = (double)(*win)[n];
                atom_c[i][n] = cos(ph) * ww;
                atom_s[i][n] = sin(ph) * ww;
            }
        }
    }
    // xw is already windowed ONCE; the atoms above carry a second window factor so
    // that <x*w, a*w> matches the Hann-DTFT convention the closed forms use.

    // ---- right-hand side: r_c = sum xw[n]cos(w n), r_s = sum xw[n]sin(w n).
    // Oscillator recurrence, re-seeded periodically so phase error cannot accumulate
    // over 4096 samples.
    vector<double> rhs(2 * K, 0.0);
    const int RESEED = 512;
    for (int i = 0; i < K; i++) {
        double ac = 0.0, as = 0.0;
        if (q[i] != 0.0) {                       // chirped: direct, no recurrence
            for (int n = 0; n < N; n++) {
                double d = (double)n - half;
                double ph = w[i] * n + q[i] * d * d;
                ac += xw[n] * cos(ph); as += xw[n] * sin(ph);
            }
        } else {
            double cw = cos(w[i]), sw = sin(w[i]);
            double c = 1.0, s = 0.0;
            for (int n = 0; n < N; n++) {
                if ((n & (RESEED - 1)) == 0) { c = cos(w[i] * n); s = sin(w[i] * n); }
                ac += xw[n] * c; as += xw[n] * s;
                double nc = c * cw - s * sw;
                s = s * cw + c * sw; c = nc;
            }
        }
        rhs[2 * i] = ac; rhs[2 * i + 1] = as;
    }

    // ---- diagonal blocks
    vector<double> dcc(K), dcs(K), dsc(K), dss(K);
    for (int i = 0; i < K; i++) {
        if (q[i] != 0.0) {
            double cc = 0.0, cs = 0.0, ss = 0.0;
            const vector<double> &ai = atom_c[i], &bi = atom_s[i];
            for (int n = 0; n < N; n++) { cc += ai[n]*ai[n]; cs += ai[n]*bi[n]; ss += bi[n]*bi[n]; }
            dcc[i] = cc; dcs[i] = cs; dsc[i] = cs; dss[i] = ss;
            continue;
        }
        double WcD, WsD, WcS, WsS;
        hann_dtft(0.0, N, WcD, WsD);
        hann_dtft(2.0 * w[i], N, WcS, WsS);
        dcc[i] = 0.5 * (WcD + WcS);
        dcs[i] = 0.5 * (WsS - WsD);
        dsc[i] = 0.5 * (WsS + WsD);
        dss[i] = 0.5 * (WcD - WcS);
    }

    // ---- banded off-diagonal blocks. freqs_hz is not sorted, so index by a sorted
    // order and only pair partials within band_bins of each other.
    vector<int> ord(K);
    for (int i = 0; i < K; i++) ord[i] = i;
    sort(ord.begin(), ord.end(), [&](int a, int b) { return freqs_hz[a] < freqs_hz[b]; });
    double band_hz = band_bins * (double)sr / (double)N;

    struct Pair { int i, j; double cc, cs, sc, ss; };
    vector<Pair> pairs;
    pairs.reserve((size_t)K * 8);
    for (int oi = 0; oi < K; oi++) {
        int i = ord[oi];
        int used = 0;
        for (int oj = oi + 1; oj < K && used < max_neighbors; oj++) {
            int j = ord[oj];
            if (freqs_hz[j] - freqs_hz[i] > band_hz) break;
            Pair p;
            p.i = i; p.j = j;
            if (q[i] != 0.0 || q[j] != 0.0) {
                const vector<double> &ci = atom_c[i], &si_ = atom_s[i];
                const vector<double> &cj = atom_c[j], &sj = atom_s[j];
                double cc = 0.0, cs = 0.0, sc = 0.0, ss = 0.0;
                for (int n = 0; n < N; n++) {
                    cc += ci[n]*cj[n]; cs += ci[n]*sj[n];
                    sc += si_[n]*cj[n]; ss += si_[n]*sj[n];
                }
                p.cc = cc; p.cs = cs; p.sc = sc; p.ss = ss;
            } else {
                double WcD, WsD, WcS, WsS;
                hann_dtft(w[i] - w[j], N, WcD, WsD);
                hann_dtft(w[i] + w[j], N, WcS, WsS);
                p.cc = 0.5 * (WcD + WcS);
                p.cs = 0.5 * (WsS - WsD);
                p.sc = 0.5 * (WsS + WsD);
                p.ss = 0.5 * (WcD - WcS);
            }
            pairs.push_back(p);
            used++;
        }
    }

    // ---- starting point / regularisation centre = the existing independent estimate
    vector<double> x0(2 * K);
    for (int i = 0; i < K; i++) {
        double A = 4.0 * mags[i] / (double)N;
        x0[2 * i]     =  A * cos(phases[i]);
        x0[2 * i + 1] = -A * sin(phases[i]);
    }

    double dmean = 0.0;
    for (int i = 0; i < K; i++) dmean += 0.5 * (dcc[i] + dss[i]);
    dmean /= (double)K;
    const double lam = reg * dmean;

    // y = (G + lam I) v
    auto matvec = [&](const vector<double>& v, vector<double>& y) {
        for (int i = 0; i < K; i++) {
            double a = v[2 * i], b = v[2 * i + 1];
            y[2 * i]     = dcc[i] * a + dcs[i] * b + lam * a;
            y[2 * i + 1] = dsc[i] * a + dss[i] * b + lam * b;
        }
        for (const Pair &p : pairs) {
            double ai = v[2 * p.i], bi = v[2 * p.i + 1];
            double aj = v[2 * p.j], bj = v[2 * p.j + 1];
            y[2 * p.i]     += p.cc * aj + p.cs * bj;
            y[2 * p.i + 1] += p.sc * aj + p.ss * bj;
            // transposed block: Gcc/Gss symmetric, Gcs <-> Gsc swap
            y[2 * p.j]     += p.cc * ai + p.sc * bi;
            y[2 * p.j + 1] += p.cs * ai + p.ss * bi;
        }
    };

    // rhs of the regularised system: r + lam*x0
    vector<double> b(2 * K);
    for (int k = 0; k < 2 * K; k++) b[k] = rhs[k] + lam * x0[k];

    // ---- conjugate gradient from x0
    vector<double> x = x0, r(2 * K), p(2 * K), Ap(2 * K);
    matvec(x, Ap);
    double rr = 0.0;
    for (int k = 0; k < 2 * K; k++) { r[k] = b[k] - Ap[k]; p[k] = r[k]; rr += r[k] * r[k]; }
    double rr0 = rr;
    for (int it = 0; it < max_iters && rr > 1e-14 * rr0; it++) {
        matvec(p, Ap);
        double pAp = 0.0;
        for (int k = 0; k < 2 * K; k++) pAp += p[k] * Ap[k];
        if (!(pAp > 0.0)) break;                     // lost positive-definiteness
        double alpha = rr / pAp;
        double rr_new = 0.0;
        for (int k = 0; k < 2 * K; k++) {
            x[k] += alpha * p[k];
            r[k] -= alpha * Ap[k];
            rr_new += r[k] * r[k];
        }
        double beta = rr_new / rr;
        for (int k = 0; k < 2 * K; k++) p[k] = r[k] + beta * p[k];
        rr = rr_new;
    }

    // ---- write back in the engine's conventions
    for (int i = 0; i < K; i++) {
        double a = x[2 * i], bb = x[2 * i + 1];
        double A = sqrt(a * a + bb * bb);
        if (!std::isfinite(A)) continue;             // never emit NaN into a track
        mags[i] = A * (double)N / 4.0;
        phases[i] = wrap_phase(atan2(-bb, a));
    }
}

// ===================== Pitch-synchronous (vibrato-demodulated) analysis =========
// The fixed-window STFT mis-measures partials that MOVE under vibrato (the Female
// "double voice": relative phases scrambled -> vibrato_rig.py shape_corr 0.47). Fix:
// warp time so the fundamental is constant (partials become stationary, which the
// analysis reconstructs at ~1.0), run the WHOLE engine on the warped signal, then
// un-warp the output. This is a wrapper around the unchanged engine; validated in
// Python (rig 0.469->0.998, real Female unity 0.706->0.995). Enabled by pitchSync
// arg 12; default 0 leaves the signal untouched.

// Compact iterative radix-2 complex FFT. sign=-1 forward, +1 inverse (unnormalized).
static void fft_radix2(vector<complex<double>>& a, int sign) {
    int n = (int)a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = sign * 2.0 * M_PI / len;
        complex<double> wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            complex<double> w(1.0, 0.0);
            for (int k = 0; k < len / 2; k++) {
                complex<double> u = a[i + k];
                complex<double> v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

// Auto-engage gates. Pitch-sync should run permanently but only ACT on monophonic
// pitched material that is actually modulated (a singer / solo instrument with
// vibrato or melody). Two cheap measures from per-frame FFT autocorrelation:
//   clarity = median periodicity (peak/zero-lag) -> monophonicity. Low for
//             polyphony (chords, mixes) and percussion.
//   span    = p95-p5 of f0 in cents over periodic frames -> pitch movement. ~0
//             for steady tones (saw/sine/held notes), large for singing.
// Engage only when clarity AND span are both high; otherwise the caller returns
// false and the signal passes through UNWARPED (an exact passthrough).
static void pitchsync_gates(const vector<float>& x, int sr,
                            double& clarity_out, double& span_out) {
    const int win = 2048, hop = 512;
    int lo = sr / 500, hi = std::min(win - 1, sr / 70);   // 70..500 Hz period search
    vector<complex<double>> buf(win);
    vector<double> cls, f0s;
    for (int s = 0; s + win <= (int)x.size(); s += hop) {
        double e = 0.0;
        for (int i = 0; i < win; i++) e += (double)x[s + i] * x[s + i];
        if (sqrt(e / win) < 3e-3) continue;                // skip near-silence
        for (int i = 0; i < win; i++) buf[i] = complex<double>((double)x[s + i], 0.0);
        fft_radix2(buf, -1);
        for (int i = 0; i < win; i++) buf[i] = complex<double>(norm(buf[i]), 0.0);
        fft_radix2(buf, +1);                                // autocorrelation (real part)
        double ac0 = buf[0].real();
        if (ac0 <= 0.0) continue;
        double peak = -1.0; int plag = lo;
        for (int L = lo; L < hi; L++) if (buf[L].real() > peak) { peak = buf[L].real(); plag = L; }
        double clarity = peak / ac0;
        cls.push_back(clarity);
        if (clarity > 0.6 && plag > 0) f0s.push_back((double)sr / plag);
    }
    clarity_out = 0.0; span_out = 0.0;
    if (cls.empty()) return;
    { vector<double> t(cls); size_t m = t.size() / 2;
      std::nth_element(t.begin(), t.begin() + m, t.end()); clarity_out = t[m]; }
    if (f0s.size() < 6) return;
    vector<double> t(f0s); size_t m = t.size() / 2;
    std::nth_element(t.begin(), t.begin() + m, t.end()); double med = t[m];
    if (med <= 0.0) return;
    vector<double> cents(f0s.size());
    for (size_t i = 0; i < f0s.size(); i++) cents[i] = 1200.0 * log2(f0s[i] / med);
    vector<double> c5(cents), c95(cents);
    size_t i5 = (size_t)(0.05 * cents.size()), i95 = (size_t)(0.95 * cents.size());
    std::nth_element(c5.begin(), c5.begin() + i5, c5.end());
    std::nth_element(c95.begin(), c95.begin() + i95, c95.end());
    span_out = c95[i95] - c5[i5];
}

// Track the fundamental's unwrapped instantaneous phase (band-limited analytic
// signal via FFT = bandpass + Hilbert in one step), and from it build tau[i] =
// the WARPED sample position of each original sample i, such that the fundamental
// advances at a constant f0_ref in warped time. Returns false (skip warp) if no
// plausible fundamental is found in 80-400 Hz.
static bool pitchsync_analyze(const vector<float>& x, int sr, vector<double>& tau,
                              bool force = false) {
    int M = (int)x.size();
    if (M < 4096) return false;

    // Auto-engage gates: only act on monophonic (clarity) + modulated (span)
    // material. Everything else returns false here -> exact passthrough. Thresholds
    // overridable via env for calibration. force=true skips the gates (manual on).
    double clarity = 0.0, span = 0.0;
    pitchsync_gates(x, sr, clarity, span);
    double clar_thr = 0.80, span_thr = 30.0;
    if (const char* e = getenv("PS_CLAR")) clar_thr = atof(e);
    if (const char* e = getenv("PS_SPAN")) span_thr = atof(e);
    bool engage = (clarity >= clar_thr && span >= span_thr);
    if (getenv("PS_DEBUG"))
        fprintf(stderr, "[PS] clarity=%.3f span=%.1f cents -> %s%s\n",
                clarity, span, (force ? "FORCED" : (engage ? "ENGAGE" : "pass")),
                (force && !engage) ? " (gates would pass)" : "");
    if (!force && !engage) return false;

    int N = 1; while (N < M) N <<= 1;
    vector<complex<double>> A(N, complex<double>(0.0, 0.0));
    for (int i = 0; i < M; i++) A[i] = complex<double>((double)x[i], 0.0);
    fft_radix2(A, -1);

    // coarse f0: strongest bin in 80-400 Hz
    int klo0 = std::max(1, (int)floor(80.0 * N / sr));
    int khi0 = std::min(N / 2 - 1, (int)ceil(400.0 * N / sr));
    int kpk = klo0; double best = -1.0;
    for (int k = klo0; k <= khi0; k++) {
        double m = norm(A[k]);
        if (m > best) { best = m; kpk = k; }
    }
    double f0c = (double)kpk * sr / N;
    if (f0c < 50.0 || f0c > 500.0) return false;

    // band-limited analytic signal: keep +freqs in [0.6,1.6]*f0c (doubled), zero else
    int klo = std::max(1, (int)floor(0.6 * f0c * N / sr));
    int khi = std::min(N / 2 - 1, (int)ceil(1.6 * f0c * N / sr));
    for (int k = 0; k < N; k++) {
        if (k >= klo && k <= khi) A[k] *= 2.0;
        else A[k] = complex<double>(0.0, 0.0);
    }
    fft_radix2(A, +1);
    double invN = 1.0 / (double)N;

    // fundamental-band amplitude + unwrapped instantaneous phase per sample
    vector<double> amp(M), phi(M);
    amp[0] = abs(A[0]) * invN;
    double prev_raw = atan2(A[0].imag(), A[0].real());
    double acc = prev_raw; phi[0] = acc;
    for (int i = 1; i < M; i++) {
        amp[i] = abs(A[i]) * invN;
        double raw = atan2(A[i].imag(), A[i].real());   // scale cancels in atan2
        double d = raw - prev_raw;
        while (d > M_PI) d -= 2.0 * M_PI;
        while (d < -M_PI) d += 2.0 * M_PI;
        acc += d; phi[i] = acc; prev_raw = raw;
    }

    // Smoothed fundamental envelope (~8 ms box) for a voicing decision. An inhale
    // / breath / consonant has no fundamental, so its analytic phase is just noise;
    // warping there resamples the breath erratically (the "crunchy inhale"). We
    // gate those spans to pass through UNWARPED.
    vector<double> ps(M + 1, 0.0);
    for (int i = 0; i < M; i++) ps[i + 1] = ps[i] + amp[i];
    int rad = std::max(1, sr / 125);
    vector<double> amp_s(M, 0.0);
    for (int i = 0; i < M; i++) {
        int lo = std::max(0, i - rad), hi = std::min(M, i + rad + 1);
        amp_s[i] = (ps[hi] - ps[lo]) / (double)(hi - lo);
    }
    // voicing threshold = fraction of a high percentile of the envelope
    double ref_amp;
    { vector<double> tmp(amp_s); int q = std::min(M - 1, (int)(0.90 * M));
      std::nth_element(tmp.begin(), tmp.begin() + q, tmp.end()); ref_amp = tmp[q]; }
    double thr = 0.12 * ref_amp;

    // f0_ref from voiced samples only (mean instantaneous frequency)
    double sumf = 0.0; long cntf = 0;
    for (int i = 1; i < M; i++) if (amp_s[i] > thr) {
        double inst = (phi[i] - phi[i - 1]) * sr / (2.0 * M_PI);
        if (inst > 0.5 * f0c && inst < 2.0 * f0c) { sumf += inst; cntf++; }
    }
    double f0ref = (cntf > 0) ? sumf / (double)cntf : f0c;
    if (!(f0ref > 1.0)) return false;
    double scale = (double)sr / (2.0 * M_PI * f0ref);

    // voiced weight w in [0,1], smoothed (~4 ms) so voiced<->unvoiced transitions
    // ramp instead of stepping
    vector<double> wv(M, 0.0);
    for (int i = 0; i < M; i++) wv[i] = (amp_s[i] > thr) ? 1.0 : 0.0;
    for (int i = 0; i < M; i++) ps[i + 1] = ps[i] + wv[i];
    int wr = std::max(1, sr / 250);
    vector<double> w(M, 0.0);
    for (int i = 0; i < M; i++) {
        int lo = std::max(0, i - wr), hi = std::min(M, i + wr + 1);
        w[i] = (ps[hi] - ps[lo]) / (double)(hi - lo);
    }

    // integrate the gated warp rate -> tau (warped-sample positions). Voiced: rate
    // ~ f0(t)/f0_ref (removes vibrato). Unvoiced: rate 1 (identity, breath passes
    // through). Fully-voiced signals (rig, sustained note) are unchanged.
    tau.assign(M, 0.0);
    for (int i = 1; i < M; i++) {
        double rate = (phi[i] - phi[i - 1]) * scale;
        if (rate < 0.5) rate = 0.5;
        if (rate > 2.0) rate = 2.0;
        double r = w[i] * rate + (1.0 - w[i]);
        tau[i] = tau[i - 1] + r;
    }
    for (int i = 1; i < M; i++) if (tau[i] < tau[i - 1]) tau[i] = tau[i - 1];
    return true;
}

// Resample x onto the uniform warped-time grid (vibrato removed): xw[j] = x at the
// original position whose warped coordinate is j (linear interp of j through tau).
// Catmull-Rom cubic resample of v at fractional position (i + t), t in [0,1).
// 4-point local, no dependencies; ~50x lower resampling error than linear on
// bright content (linear interp of the ~4% warp resample was audible as a rattle
// on bright/high sung notes — sustained/dark notes hid it). Ends clamp to edge.
static inline float catmull_rom(const vector<float>& v, int i, double t) {
    int n = (int)v.size();
    auto at = [&](int k) { return (double)v[std::max(0, std::min(n - 1, k))]; };
    double p0 = at(i - 1), p1 = at(i), p2 = at(i + 1), p3 = at(i + 2);
    double t2 = t * t, t3 = t2 * t;
    return (float)(0.5 * (2.0 * p1 + (-p0 + p2) * t
                          + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
                          + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3));
}

static vector<float> pitchsync_warp(const vector<float>& x, const vector<double>& tau) {
    int M = (int)x.size();
    int W = (int)floor(tau[M - 1]);
    vector<float> xw(std::max(0, W), 0.0f);
    int i = 0;
    for (int j = 0; j < W; j++) {
        while (i + 1 < M && tau[i + 1] < (double)j) i++;
        if (i + 1 >= M) { xw[j] = x[M - 1]; continue; }
        double denom = tau[i + 1] - tau[i];
        double frac = denom > 1e-12 ? ((double)j - tau[i]) / denom : 0.0;
        xw[j] = catmull_rom(x, i, frac);
    }
    return xw;
}

// Un-warp the (warped-domain) output back to real time: out[i] = yw at warped
// position tau[i] (linear interp on the uniform warped grid).
static void pitchsync_unwarp(const vector<float>& yw, const vector<double>& tau,
                             vector<float>& out) {
    int M = (int)tau.size();
    int W = (int)yw.size();
    out.assign(M, 0.0f);
    if (W == 0) return;
    for (int i = 0; i < M; i++) {
        double pos = tau[i];
        if (pos <= 0.0) { out[i] = yw[0]; continue; }
        int j = (int)floor(pos);
        if (j >= W - 1) { out[i] = yw[W - 1]; continue; }
        double frac = pos - (double)j;
        out[i] = catmull_rom(yw, j, frac);
    }
}
// ================================================================================


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
    // transientShapeThresholdDB: see transientNegotiationTactics. How far a frame's
    // level-blind spectral SHAPE change must stand out above the local median before
    // it counts as an event, in dB. Adaptive rather than absolute because a dense mix
    // changes shape constantly: an absolute threshold of 6 dB took take-me-out from
    // 0.3% to 83.6% of samples served by the transient original-blend while its real
    // (passthrough-excluded) SRR did not move at all -- a pure measurement artefact.
    // At 3 dB above the local median, both Fairlights gain +2.6 dB of REAL SRR and
    // every other file in the corpus is untouched. 0 = off (onset detection only).
    float transientShapeThresholdDB = 3.0f;
    int pitch_shift_semi = 0;
    // Synthesis engine: 0 = legacy per-frame OLA (byte-identical to before),
    // 1 = continuous-phase oscillator bank (McAulay-Quatieri). See mq_synthesize.
    // A/B switch; MQ supports chordIntervals == {0} only (no chord/interval mode).
    int synth_mode = 0;


    double peak_sidelobe_attenuation_dB = 12.0;
    double peak_floor_below_max_dB = 60.0;
    // ===== Phase 1: amplitude / HF fidelity tuning =====
    // Asymmetric amplitude smoothing (replaces the symmetric 0.7/0.3 EMA): react
    // fast to rising partials (preserves attacks + upper harmonics), smooth slower
    // on decay. amp_smooth_* = weight on the newly measured value.
    // Set both to 0.7 to recover the old symmetric behavior.
    double amp_smooth_attack  = 0.90;
    double amp_smooth_release = 0.50;
    // Fast-release gates. The slow release blurring genuine decays caused
    // both the shifted "ghost tail" (Female gaps +2 dB) and the DrumLoop
    // "first crash never died" contrast loss, so real decays must track
    // fast — but a BLANKET fast release under shift made the engine chase
    // vibrato ripple frame by frame (the 440 vibrato-saw down5 wobbled
    // audibly, ±1.7 dB at ~1 Hz). Two discriminators; either selects fast:
    //  - drop gate (all modes): new measurement below amp_release_fast_drop
    //    x current level (-6 dB in one frame) = a real event decay;
    //  - streak gate (shift only): measurement fallen for
    //    shift_release_streak_frames consecutive frames = monotone decay
    //    (reverb tails); vibrato/ripple alternates and never builds a streak.
    double amp_release_fast = 0.85;
    double amp_release_fast_drop = 0.5;
    int shift_release_streak_frames = 3;
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
    // Shift-mode fill floor. 300 Hz was tried (to chase reverb) but the
    // deficit is nonzero AT the voice's own partials during vibrato (tonal-
    // model imperfection), so the fill painted an unshifted voice-shaped
    // noise ghost — audibly a second voice at the original pitch (+25 dB at
    // her partials vs before). The cymbal/air fill this exists for lives
    // above ~6 kHz, so a high floor keeps the win without the ghost.
    double residual_hp_hz_shift = 2500.0;
    double residual_oversub    = 1.0;    // subtract this * model magnitude
    // Residual reconstruction: 0 = stochastic (random-phase noise shaped to the
    // spectral deficit); 1 = true-phase (add the real model-subtracted residual
    // above residual_hp_hz, preserving its temporal structure -> restores the
    // voice's natural breath/consonant detail instead of decorrelated haze).
    // Mode 1 applies at unity only (unshifted true residual over a shifted tonal
    // body would double the HF at the wrong pitch). Overridable via CLI arg 9.
    int residual_mode = 0;
    // peak_birth_confirm_frames: a new peak must be matched in this many
    //   consecutive frames before it's allowed to synthesize. 1 = old behavior
    //   (immediate). 2 = one frame of confirmation (kills single-frame
    //   phantoms — the most audible musical-noise source). 3+ is more
    //   aggressive but adds onset latency proportional to the long hop.
    int peak_birth_confirm_frames = 2;
    // file_start_confirm_frames: FILE-START CREDIT. Every track is newborn in
    //   frame 0, so with 2-frame confirmation the engine renders ~nothing for
    //   the first hop or two (~30-40 ms): measured 10 ms out/in level ratios of
    //   0.08 0.03 0.09 on DrumLoop and 0.00 0.15 0.72 on a steady sine, and on
    //   files that begin on a hit that gap was 80-90% of the whole-file
    //   residual (DrumLoop SRR 8.9 whole vs 18.5 interior). Confirmation
    //   exists to reject single-frame phantoms that flicker against the
    //   PREVIOUS frame; in frame 0 there is no previous frame, so a peak born
    //   there is not that failure mode -- but a noise-floor peak born in a
    //   quiet first frame IS (confirming frame-0 births outright put a -38 dB
    //   burst at t=0 on the choir, whose input starts at -64 dB). So the
    //   credit is RETROACTIVE: a track born in a frame < this value is held
    //   back, and if it passes normal confirmation later (i.e. it is still
    //   there in the next frame) its held-back frame-0 measurement is added to
    //   that frame's render list. Analysis completes before synthesis, so this
    //   is exact lookahead, not a heuristic: the same tracks render as before,
    //   they just render from the frame they were first seen in. 0 = off.
    //   Env override FILE_START_CONFIRM for A/B.
    int file_start_confirm_frames = 1;
    // file_start_cap_ms: the file-start credit above releases frame-0 births, and
    //   frame 0 uses a RECT-to-Hann synthesis window (so a file that begins mid-
    //   note is not faded in). Together those render a partial at full amplitude
    //   from sample 0 even when the note it was measured from starts later in the
    //   frame: on the Female line the input is -83.7 dBFS at t=0 and the engine
    //   emitted -28.3, thirty milliseconds before the real onset. A listener
    //   identified that blind as "a double hit in the first utterance" (docs 5b.2).
    //   The engine must not emit energy before the input does, so cap the output's
    //   short-time level to the input's over this opening window. Same principle
    //   as residual_env_cap, and scoped to the only place the rect window exists.
    //   0 = off.
    double file_start_cap_ms = 120.0;
    double file_start_cap_headroom_db = 3.0;
    // transient_short_amp: at a transient the engine renders 256-sample (5.3 ms)
    //   frames, but takes BOTH frequency and amplitude for them from a
    //   2048-point long-window FFT centred on the frame (AnalysisInfo.cpp
    //   m_transientLongSpecs, added so the spectrum evolves through the hit).
    //   A 2048 window is 42.7 ms, so an amplitude read from it is smeared by
    //   +-21 ms -- and 21 ms is exactly the measured pre-onset smear on
    //   shifted drums (docs 4d.6b). Frequency genuinely needs the long window
    //   (256 bins = 187 Hz); amplitude does not. 1 = rescale each transient
    //   frame's long-window amplitudes by the ratio of the frame's OWN
    //   short-window energy to the long window's, so the spectral shape and
    //   frequencies stay long-window but the LEVEL follows the 5.3 ms frame.
    //   Env TRANS_SHORT_AMP.
    int transient_short_amp = 1;
    // ===== Chirped (linear-FM) estimation and synthesis =====
    // chirp_mode: a partial that MOVES is not a sinusoid over a 42.7 ms frame.
    //   On three voices with independent vibrato (battery/polyvoice_rig.py, the
    //   choir case) a stationary basis at the TRUE frequencies caps at 24.4 dB
    //   while a chirped basis reaches 46.2 -- 22 dB of the gap is stationarity
    //   alone. Each track gets a frequency slope df/dt from a centred difference
    //   of its own trajectory (available post-tracking, exact lookahead), and
    //   that slope adds a quadratic term to the phase: 2pi(f*t + df*t^2/2).
    //   The chirp MUST be in both the solve and the render -- chirping only one
    //   measures worse than chirping neither (1.4 and 20.0 dB vs 24.4), the same
    //   law that separated joint_mode 1 from 2. 0 = off, 1 = on.
    //
    //   DEFAULT OFF, and the reason is the analysis window, not this code.
    //   Measured on the rig: across a 4096 window (85 ms) a vibrato partial's
    //   frequency sweeps a median 75.8 Hz -- 6.5 bins -- and after fitting the
    //   best straight line through it, 8.77 Hz of curvature REMAINS. A linear-FM
    //   atom is the wrong model at this window length, so the chirp cannot pay:
    //   the rig moves 13.11 -> 13.26 dB. The same residual at 2048 is 2.09 Hz and
    //   at 1024 is 0.52 Hz, so this becomes worth switching on only alongside a
    //   shorter analysis window for moving partials. Kept, measured and gated
    //   rather than deleted, because it is the validated half of that pair.
    //   Env CHIRP / CHIRP_MIN_SLOPE.
    int chirp_mode = 0;
    // residual_env_cap: the stochastic fill is built from 1024-sample (21.3 ms)
    //   STFT frames, so at a sharp onset it spreads the hit's noise energy about
    //   +-10 ms. On a pure click that IS the whole remaining defect: rig pre-echo
    //   is -23.8 dB with the fill switched off and -6.5 dB with it on, starting
    //   15 ms before the hit. The tonal model is not the problem there.
    //   Noise carries no phase structure worth protecting, so shape it in the
    //   time domain instead: the fill may never carry more local energy than the
    //   input does at that instant. Envelope measured over a short moving RMS;
    //   where the input is silent the fill is silent. 0 = off.
    //   Env RESID_ENV_CAP / RESID_ENV_MS.
    int residual_env_cap = 1;
    // residual_shift_mode: where the stochastic noise fill sits under pitch shift.
    //   0 = at the ORIGINAL pitch (the shipped behaviour, 4d85b5f: room, breath and
    //       air do not transpose when you move a note).
    //   1 = TRANSPOSED with the tonal content.
    //   A listener reports up-shift worse than down-shift while eight objective
    //   measures say the opposite (docs 5c). One candidate explanation is that this
    //   choice inverts the natural arrangement: down-shift leaves the air ABOVE the
    //   harmonics, as real sources are built, while up-shift moves the harmonics up
    //   PAST a stationary noise bed. The residual is only -36 dB on the drum loop
    //   overall but -18 dB on cymbals, and a drum loop's top end is cymbals. This
    //   knob exists to A/B that hypothesis, not as a settled improvement.
    //   Env RESID_SHIFT.
    int residual_shift_mode = 0;
    double residual_env_ms = 3.0;
    // Below this |slope| a track is treated as stationary, so material without
    // vibrato keeps the closed-form Gram and pays nothing. Set from measurement,
    // not theory: per-frame frequency jitter of ~1 Hz over a 42.7 ms centred
    // difference is already ~23 Hz/s, so a low threshold chirps STEADY partials
    // on estimator noise -- at 20 Hz/s the rig's steady case fell 30.14 -> 24.51
    // dB. At 300 Hz/s it is back to 30.13 while real vibrato (200-2700 Hz/s
    // on the rig) still qualifies.
    double chirp_min_slope_hz_s = 300.0;
    //   NOT generalised to every birth (tried Sep 21): releasing each onset's
    //   first observation renders the frame that only PARTLY contains the hit,
    //   which adds pre-echo (DrumLoop up5 pre-onset +7.7 -> +9.1 dB) and does
    //   not sharpen the attack (+4.0 -> +4.5 ms). The shifted-attack softness
    //   is the 4096-frame smear, not confirmation latency.
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
    // between partials) and overshoots the limiter. 28 (was 35): the
    // SaintSaens onset flutter was flickering 20-23 Hz detections at
    // -30..-35 dB rel max — each birth/death is a 341 ms infrasonic thump —
    // while the real pedal partials sit at -9..-25 rel.
    double lf_floor_below_max_db = 28.0;
    // Ignore LF candidates below this frequency: the first bins of the LF
    // FFT are DC drift / subsonic leakage, never a playable partial, and
    // their flicker was the other half of the onset flutter.
    double lf_min_hz = 15.0;
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
    // ===== Joint amplitude/phase estimation (jointMode, CLI arg 13) =====
    // See joint_amp_phase() above for the mechanism and the measured headroom.
    // joint_band_bins: partials within this many 4096-FFT bins of each other are
    //   solved together. The Hann main lobe is +/-2 bins; beyond ~6 bins the
    //   coupling is below -50 dB, so 8 is generous. Larger = more coupling
    //   captured, but more chance of the solve explaining unmodelled energy with
    //   the wrong partials.
    double joint_band_bins = 8.0;
    // joint_reg: Tikhonov weight as a fraction of the mean diagonal, centred on
    //   the INDEPENDENT estimate — so large values reproduce the current engine
    //   and 0 is the unconstrained joint solve. This is the safety knob: raise it
    //   if a class gets worse.
    double joint_reg = 1e-3;
    // joint_iters: conjugate-gradient iteration cap. The system is strongly
    //   diagonally dominant for well-separated partials, so this converges fast;
    //   the cap only bites on dense, closely-spaced frames.
    int joint_iters = 60;
    // joint_smooth: re-apply the tracker's asymmetric amplitude EMA to the jointly
    //   solved amplitudes (joint_mode 2 only, which runs after tracking and so
    //   would otherwise emit raw per-frame measurements). 1 = on.
    //   Mode 2 = MEDIAN-OF-3 instead. The EMA is a low-pass on the amplitude
    //   trajectory, so it cannot tell estimator noise from real modulation and
    //   attenuates both: ablating it entirely gains +2.3 dB mean residual SRR
    //   (Fairlight C2 +5.3, Piano +4.5, take-me-out +3.3) at the price of the
    //   variance it was controlling. A median-of-3 kills single-frame spikes
    //   while passing monotone ramps and periodic modulation exactly -- the same
    //   argument, and the same fix, that e5c2b41 applied to the FREQUENCY
    //   trajectory -- but MEASURED WORSE (-4.3 dB mean): amplitude modulates far
    //   faster than frequency does (Fairlight carries ~17.6 Hz AM against a
    //   46.9 Hz frame rate, under three frames per cycle), so a median-of-3
    //   destroys real modulation instead of preserving it. Kept as a knob, not
    //   a candidate.
    //   Mode 3 = LEVEL-GATED EMA, the one that works: smoothing exists to
    //   control estimator variance, and variance is a function of a partial's
    //   own SNR. Partials within joint_smooth_raw_db of the frame's loudest are
    //   measured well enough to pass through raw; quieter ones keep the EMA.
    //   0 = raw, 1 = EMA (legacy), 2 = median-of-3, 3 = level-gated EMA.
    int joint_smooth = 1;
    double joint_smooth_raw_db = 20.0;
    // ===== Harmonic phase coherence (shift mode) =====
    // Under pitch shift each track is an INDEPENDENT phase-propagated oscillator,
    // so per-partial frequency errors integrate into relative-phase drift: every
    // partial's own amplitude stays flat while the summed waveform slowly morphs,
    // and a morphing crest factor is heard as the level moving around. Measured by
    // waveform_shape_consistency: 300hzSaw is 0.999 at unity but 0.79 up5 / 0.82
    // down5 (the source is 0.999).
    //
    // Fix: group each frame's tracks into harmonic stacks and let members derive
    // frequency and phase from their stack ROOT's trajectory, so a stack's relative
    // phases are rigid by construction. The joint estimator cannot help here — it
    // improves the phase MEASUREMENT and the shift path never uses it.
    //
    // shift_harmonic_lock: 0 disables (falls back to independent propagation).
    int shift_harmonic_lock = 1;
    // harmonic_lock_tol: a member must sit within this RELATIVE tolerance of
    //   k*f_root (the test is |r - k| < tol*k, i.e. |r/k - 1| < tol). 0.4% is tight
    //   enough that a piano's stretched partials (~7.7% sharp by k=20 for typical
    //   inharmonicity) and chord intervals (SaintSaens 78.7/39.6 = 1.987, 0.6% off
    //   the octave) both FAIL and stay independent — only genuine harmonic stacks
    //   lock. Widening this past ~1% starts capturing piano and must not be done
    //   without re-checking the piano and chord guards. Measured sweep at +5
    //   semitones (shape_consistency): the 300 Hz saw reaches 0.9990 at EVERY
    //   tolerance from 0.0005 up, because its harmonics are exact — but the piano
    //   degrades progressively (lock-off 0.7442; 0.0005 -> 0.7453 untouched,
    //   0.002 -> 0.7082, 0.004 -> 0.6767). A piano's low partials are only ~0.1-0.3%
    //   sharp, so a loose tolerance captures them and forces them exactly harmonic.
    //   0.0005 buys the entire saw win with no measurable effect on the piano.
    double harmonic_lock_tol = 0.0005;
    // harmonic_root_max_hz: only look for stack roots below this. A root is a
    //   fundamental, not an upper partial.
    double harmonic_root_max_hz = 1200.0;
    // ----- root validity (added after the 440sawtooth artifact) -----
    // The grouping used to take the lowest-frequency track as root with no further
    // test. On 440sawtooth down5 that was sub-audio LF junk at 21.7 Hz (amplitude
    // rank 463, 56 dB below the loudest partial), and hundreds of real partials were
    // quantised onto its grid. These three tests reject that while leaving a genuine
    // fundamental (300hzSaw's root: 300 Hz, amplitude rank 1) untouched.
    // harmonic_root_min_hz: below this a "root" is rumble, not a fundamental. The LF
    //   tier tracks down to lf_min_hz = 15 Hz, which is what supplied the bad root.
    double harmonic_root_min_hz = 50.0;
    // harmonic_root_min_rel_db: the root must be within this many dB of the frame's
    //   loudest partial. A fundamental carries real energy; the bad root was -56 dB.
    double harmonic_root_min_rel_db = 40.0;
    // harmonic_lock_max_dev_hz: absolute cap on |f_member - k*f_root|, on top of the
    //   relative tolerance. 0.05% of 20 kHz is 10 Hz and measured drags reached 11 Hz;
    //   a genuine harmonic sits far closer than that to k*f0.
    double harmonic_lock_max_dev_hz = 2.0;
    // harmonic_min_members: a stack must have at least this many members, AND at
    //   least one of them at k = 2 or 3 (see the harmonic-support test).
    int harmonic_min_members = 3;
    // harmonic_lock_warmup_frames: how many frames a root must have been
    //   frequency-steady before its stack may lock. Originally this reused the
    //   steady-tone lock's 8 frames, which meant the lock switched on mid-note at
    //   t = 0.17 s and the shape converged to the locked one over the next ~5
    //   frames — a crest excursion the limiter turned into an audible 1-2 dB dip.
    //   With roots now validated the warm-up is no longer needed, and removing it
    //   is measurably better: engaging from the first frame means there is no
    //   mid-note mode switch at all. On 440sawtooth down5 the worst level dip goes
    //   -2.18 dB -> -0.61 dB and the crest factor 2.187 -> 1.983 against an input
    //   crest of 1.731 (an ideal sawtooth is sqrt(3) = 1.732), while the shape
    //   metric is unchanged. 0 = engage as soon as a valid stack is found.
    int harmonic_lock_warmup_frames = 0;
    // harmonic_stack_min_energy_frac: a stack (root + members) must carry at least
    //   this fraction of the frame's total partial energy before it may lock.
    //   Rationale: the lock exists to stop ONE harmonic series drifting apart. On a
    //   dense polyphonic mix there is no single series — many overlap, the greedy
    //   grouping is speculative, and locking cost 17 battery regressions at
    //   warmup 0 (vs 9 with a warm-up). Gating on how much of the frame the stack
    //   actually explains targets the cases the lock is FOR (a saw is ~1.0 of the
    //   frame; one stack in take-me-out is a small slice) without a time delay,
    //   which is what reintroduced the audible onset step.
    double harmonic_stack_min_energy_frac = 0.30;
    // harmonic_rps_mode: RELATIVE PHASE SHIFT synthesis (Saratxaga 2009).
    //   The shipped lock captures theta_k = phi_k - k*phi_root ONCE and freezes it, so a
    //   member whose frequency is even slightly off k*f_root drifts away from its root --
    //   which is why the tolerance has to be 0.05% AND under 2 Hz, and why the lock renders
    //   only 2.5% of DrumLoop's track-frames, 2.7% of the Female's and 5.4% of the Piano's
    //   (HLOCK_STATS; 69% on the synthetic saw it was validated on). 97% of real partials
    //   therefore free-run, which is the crest-drift and warble mechanism (docs 5e).
    //   RPS re-measures theta_k from the ANALYSIS phases EVERY frame instead, so a member's
    //   own frequency deviation is absorbed rather than accumulated. Two consequences:
    //     - the frequency tolerance can be much looser (harmonic_lock_tol_rps), and
    //     - the member keeps its OWN frequency. The shipped lock also quantises members to
    //       exactly k*f_root, which is wrong for stretched-partial sources like piano; RPS
    //       only re-anchors the PHASE each frame and leaves pitch alone.
    //   0 = shipped frozen-offset lock, 1 = RPS. Env RPS / RPS_TOL.
    int harmonic_rps_mode = 0;
    double harmonic_lock_tol_rps = 0.02;      // 2%, vs 0.0005 frozen
    // RPS re-anchors a member as phi_k = k*phi_root + theta_k, so any error in the ROOT's
    // propagated phase is multiplied by k. The frozen lock has the same exposure, which is
    // why it needed a steady root. Cap k to bound the amplification. Env RPS_MAXK.
    int harmonic_rps_max_k = 0;               // 0 = no cap
    //End of user settings
    
    
    AppSettings settings;
    
    
    /* parse input arguments */
    if (!parseArgs(argc, argv, settings))
    {
        printUsage();
        for(int j = 1; j < argc; j++)
            printf("%s\n", argv[j]);
        return -1;
    }
    // Optional CLI overrides. Absent => 0 = unity / OLA, so behavior is
    // identical to the legacy invocation.
    pitch_shift_semi = settings.pitchShiftSemi;
    synth_mode = settings.synthMode;
    // Diagnostics (bisect the "phasey/hollow double voice"):
    residual_noise_gain *= settings.residualScale; // 0 => stochastic residual off
    int unity_dedup = settings.unityDedup;         // run dedup at unity too
    residual_mode = settings.residualMode;         // 0 = stochastic, 1 = true-phase
    if (settings.residualHpHz >= 0.0) {            // diagnostic hp override (both bands)
        residual_hp_hz = settings.residualHpHz;
        residual_hp_hz_shift = settings.residualHpHz;
    }
    int amp_mode = settings.ampMode;               // 0 = parabolic peak, 1 = energy-integrated
    int phase_mode = settings.phaseMode;           // 0 = stationary phase, 1 = reassigned (chirp-aware)
    int joint_mode = settings.jointMode;           // 0 = per-peak amp/phase, 1 = joint least squares
    // JOINT_* env overrides for knob sweeps (diagnostics only; unset = the defaults
    // in the user-settings block).
    if (const char* e = getenv("BIRTH_CONFIRM")) peak_birth_confirm_frames = atoi(e); // diagnostic
    if (const char* e = getenv("FILE_START_CONFIRM")) file_start_confirm_frames = atoi(e);
    if (const char* e = getenv("TRANS_SHAPE_DB")) transientShapeThresholdDB = atof(e);
    if (const char* e = getenv("FILE_START_CAP_MS")) file_start_cap_ms = atof(e);
    if (const char* e = getenv("FILE_START_CAP_DB")) file_start_cap_headroom_db = atof(e);
    if (const char* e = getenv("TRANS_SHORT_AMP")) transient_short_amp = atoi(e);
    if (const char* e = getenv("CHIRP")) chirp_mode = atoi(e);
    if (const char* e = getenv("RESID_ENV_CAP")) residual_env_cap = atoi(e);
    if (const char* e = getenv("RESID_SHIFT")) residual_shift_mode = atoi(e);
    if (const char* e = getenv("RESID_ENV_MS")) residual_env_ms = atof(e);
    if (const char* e = getenv("LF_CUTOFF")) lf_cutoff_hz = atof(e);       // diagnostic
    if (const char* e = getenv("LF_MAX_FLUX")) lf_max_flux = atof(e);     // diagnostic
    if (const char* e = getenv("CHIRP_MIN_SLOPE")) chirp_min_slope_hz_s = atof(e);
    if (const char* e = getenv("JOINT_BAND"))  joint_band_bins = atof(e);
    if (const char* e = getenv("JOINT_REG"))   joint_reg       = atof(e);
    if (const char* e = getenv("JOINT_ITERS")) joint_iters     = atoi(e);
    if (const char* e = getenv("JOINT_SMOOTH")) joint_smooth   = atoi(e);
    if (const char* e = getenv("JOINT_SMOOTH_RAW_DB")) joint_smooth_raw_db = atof(e);
    if (const char* e = getenv("HLOCK"))       shift_harmonic_lock = atoi(e);
    if (const char* e = getenv("HLOCK_TOL"))   harmonic_lock_tol   = atof(e);
    if (const char* e = getenv("HLOCK_ROOT_MIN_HZ")) harmonic_root_min_hz     = atof(e);
    if (const char* e = getenv("HLOCK_ROOT_DB"))     harmonic_root_min_rel_db = atof(e);
    if (const char* e = getenv("HLOCK_MAXDEV"))      harmonic_lock_max_dev_hz = atof(e);
    if (const char* e = getenv("HLOCK_WARMUP"))      harmonic_lock_warmup_frames = atoi(e);
    if (const char* e = getenv("HLOCK_EFRAC"))       harmonic_stack_min_energy_frac = atof(e);
    if (const char* e = getenv("RPS"))               harmonic_rps_mode = atoi(e);
    if (const char* e = getenv("RPS_TOL"))           harmonic_lock_tol_rps = atof(e);
    if (const char* e = getenv("RPS_MAXK"))          harmonic_rps_max_k = atoi(e);
    if (const char* e = getenv("JOINT_NOEMA")) {                 // 1 = bypass the amplitude EMA
        if (atoi(e)) { amp_smooth_attack = 1.0; amp_smooth_release = 1.0;
                       amp_release_fast = 1.0; }
    }
    // Energy-integrated amplitude (amp_mode==1): a partial moving under vibrato
    // smears across bins, so its peak height under-reads its true amplitude
    // (steady partials reconstruct at ~1.0, vibrato partials at ~0.70). Measure
    // the integrated main-lobe energy over ±AMP_W bins instead and convert to an
    // equivalent-sinusoid amplitude, calibrated once against a unit on-bin tone
    // through the same window/FFT so it matches the parabolic scale on steady tones.
    const int AMP_W = 6;
    double amp_Eref_main = 1.0;
    {
        vector<float> cal(ANALYSIS_SIZE, 0.0f);
        int calbin = ANALYSIS_SIZE / 8;            // arbitrary on-bin frequency
        for (int i = 0; i < ANALYSIS_SIZE; i++)
            cal[i] = (float)(sin(2.0 * M_PI * calbin * i / ANALYSIS_SIZE)) * analysis_hanning[i];
        RealFFT(cal.data(), ANALYSIS_SIZE);
        double E = 0.0;
        for (int k = calbin - AMP_W; k <= calbin + AMP_W; k++) {
            double re = cal[k], im = cal[ANALYSIS_SIZE - k];
            E += re * re + im * im;
        }
        if (E > 0.0) amp_Eref_main = E;            // energy for a unit-amplitude tone
    }

    AudioBuffer inputWav;
    inputWav.mSamples = nullptr;
    
    if (!readInWaveFile(settings.inputWavFilePath, &inputWav)) {
        cerr << "Error: Not read" << endl;
        return -1;
    }
    
    vector<vector<float>> audioData = audioBufferToVector(inputWav);
    vector<float> singleChannelData = audioData[0];

    // ===== Pitch-synchronous analysis (wrapper): warp the vibrato out here, run
    // the whole engine on the warped signal, un-warp the output just before write.
    vector<double> ps_tau;              // warped-sample position per ORIGINAL sample
    if (settings.pitchSync != 0) {
        if (pitchsync_analyze(singleChannelData, sr, ps_tau, settings.pitchSync == 2)) {
            size_t before = singleChannelData.size();
            singleChannelData = pitchsync_warp(singleChannelData, ps_tau);
            cout << "[pitch-sync] warped " << before << " -> "
                 << singleChannelData.size() << " samples\n";
        } else {
            ps_tau.clear();             // no plausible fundamental -> process normally
            cout << "[pitch-sync] no fundamental found; skipping warp\n";
        }
    }
    int lengthYouNeed = inputWav.mNumSamples;
    cout << "Read " << inputWav.mNumSamples << " samples, " << inputWav.mChannels << " channels at " << inputWav.mSampleRate << " Hz.\n";
    
    int num_frames = (int)ceil((double)singleChannelData.size() / hop_size);
    singleChannelData.resize(singleChannelData.size() + frame_size, 0.0f);
    
    
    
    /**
     * Here is where we have the transients being thrown into a list. 1 for transient, 0 for nothing. Needed for window switiching
     */
    vector<float> transientList = transientNegotiationTactics(num_frames, transientThresholdDB, hop_size, LONG_SIZE, singleChannelData, transientShapeThresholdDB);
    if (getenv("TRANS_DEBUG")) {
        fprintf(stderr, "transient frames (s):");
        for (int f = 0; f < (int)transientList.size(); f++)
            if (transientList[f] == 1.0f) fprintf(stderr, " %.3f", (double)(f * hop_size) / sr);
        fprintf(stderr, "\n");
    }

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

    // File-start credit state (see file_start_confirm_frames).
    vector<pair<int, PeakTrack>> file_start_pending;
    unordered_set<int> ever_confirmed_ids;
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
            // RealFFT consumes analysis_frame in place; the joint solve needs the
            // windowed samples themselves, so keep a copy when it is enabled.
            vector<double> analysis_windowed;
            if (joint_mode != 0) {
                analysis_windowed.assign(ANALYSIS_SIZE, 0.0);
                for (int i = 0; i < ANALYSIS_SIZE; i++)
                    analysis_windowed[i] = (double)analysis_frame[i];
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

            // Reassignment tier: FFT the SAME samples through the time-weighted
            // (th) and derivative (dh) companion windows, so the phase loop below
            // can compute each partial's reassigned time/frequency (chirp-aware).
            // Only when phase_mode==1; same RealFFT packing (re=[k], im=[N-k]).
            vector<float> th_frame, dh_frame;
            if (phase_mode != 0) {
                th_frame.assign(ANALYSIS_SIZE, 0.0f);
                dh_frame.assign(ANALYSIS_SIZE, 0.0f);
                for (int i = 0; i < ANALYSIS_SIZE; i++) {
                    int idx = analysis_start + i;
                    double s = (idx >= 0 && idx < (int)singleChannelData.size())
                                   ? (double)singleChannelData[idx] : 0.0;
                    th_frame[i] = (float)(s * analysis_th[i]);
                    dh_frame[i] = (float)(s * analysis_dh[i]);
                }
                RealFFT(th_frame.data(), ANALYSIS_SIZE);
                RealFFT(dh_frame.data(), ANALYSIS_SIZE);
            }


            vector<int> peaks = detect_peaks(analysis_mag, threshold_long_analysis,
                                             sr, ANALYSIS_SIZE, hf_extra_sensitivity_db);
            filter_peaks_by_quality(analysis_mag, peaks,
                                    peak_sidelobe_attenuation_dB, peak_floor_below_max_dB,
                                    sr, ANALYSIS_SIZE, hf_extra_sensitivity_db);
            vector<double> freqs, mags;
            if (peaks.size() > 0) {
                parabolic_interpolation(analysis_mag, peaks, freqs, mags);

                // Energy-integrated amplitude (vibrato-robust): replace each
                // peak height with the equivalent-sinusoid amplitude from the
                // integrated main-lobe energy, so partials smeared by vibrato
                // keep their true amplitude. Frequency (freqs) is unchanged.
                if (amp_mode == 1) {
                    for (size_t i = 0; i < peaks.size(); i++) {
                        int b = peaks[i];
                        double E = 0.0;
                        for (int k = b - AMP_W; k <= b + AMP_W; k++)
                            if (k >= 1 && k < a_half) E += analysis_mag[k] * analysis_mag[k];
                        double A = sqrt(E / amp_Eref_main);
                        mags[i] = A * (double)ANALYSIS_SIZE / 4.0;
                    }
                }

                // Per-partial phase, corrected for the fractional bin offset:
                // the DFT phase of a symmetric window at the integer peak bin
                // carries a residue of pi*(true_bin - int_bin). Reading it raw
                // gave every harmonic a different phase offset (up to +-90°) —
                // waveform shape was scrambled even at unity (440 saw shape
                // correlation 0.80; 0.9997 with the correction in simulation).
                // Spectral reassignment (Auger-Flandrin) for chirping partials.
                // phase_mode is a bitmask so the two corrections can be A/B'd
                // independently from the CLI without recompiling:
                //   bit 0 (1) = reassigned PHASE   (chirp term on top of baseline)
                //   bit 1 (2) = reassigned FREQUENCY (w_hat replaces PV inst-freq)
                //   3 = both.  reassigned_hz[i] >= 0 means a valid w_hat exists.
                //
                // FINDING (vibrato_rig.py, Aug 2): this does NOT fix the vibrato
                // artifact. Kept opt-in (like amp_mode=1) as a validated dead end +
                // debug harness (RA_DEBUG). Baseline VIBRATO shape_corr 0.469;
                // mode1 (phase) 0.206, mode2 (freq) 0.419, mode3 0.208 — all worse.
                // Why: w_hat already matches parabolic/PV to ~1 Hz (frequency was
                // never the problem), and t_hat is large (±40..190 samples, sign
                // varying per harmonic — RA_DEBUG), so w_hat*t_hat is a multi-radian
                // per-partial term that scrambles relative phase in EITHER sign. The
                // linear-chirp model is itself wrong here: an 85 ms (4096) window
                // spans ~half a 5.5 Hz vibrato cycle, so the partial curves rather
                // than chirps linearly. Real fix is pitch-synchronous / vibrato-
                // demodulated analysis (make partials stationary before analysis).
                static const char* RA_DBG = getenv("RA_DEBUG");
                vector<double> reassigned_hz(peaks.size(), -1.0);
                vector<double> phases;
                for (size_t i = 0; i < peaks.size(); i++) {
                    int k = peaks[i];
                    // Baseline (stationary) phase: fractional-bin corrected. Exactly
                    // right for a steady partial (t_hat -> 0), keeps STEADY ~1.0.
                    double phi = wrap_phase(analysis_phase_store[k] -
                                            M_PI * (freqs[i] - peaks[i]));
                    if (phase_mode != 0 && k >= 1 && k < a_half) {
                        double pre = analysis_frame[k];
                        double pim = analysis_frame[ANALYSIS_SIZE - k];
                        double d = pre * pre + pim * pim;
                        if (d > 1e-20) {
                            double thre = th_frame[k], thim = th_frame[ANALYSIS_SIZE - k];
                            double dhre = dh_frame[k], dhim = dh_frame[ANALYSIS_SIZE - k];
                            double re_th = thre * pre + thim * pim;   // Re(Xth . conj(p))
                            double im_dh = dhim * pre - dhre * pim;   // Im(Xdh . conj(p))
                            double w_k = 2.0 * M_PI * k / (double)ANALYSIS_SIZE;
                            double w_hat = w_k - im_dh / d;           // reassigned freq (rad/sample)
                            double t_hat = re_th / d;                 // reassigned time (samples, rel. centre)
                            reassigned_hz[i] = w_hat * (double)sr / (2.0 * M_PI);
                            if (RA_DBG && frame_idx == 30 && i < 8)
                                fprintf(stderr, "RA f%d p%zu bin%d  parab=%.2f  w_hat=%.2f  t_hat=%.1f\n",
                                        frame_idx, i, k,
                                        freqs[i] * ((double)sr / (double)ANALYSIS_SIZE),
                                        reassigned_hz[i], t_hat);
                            if (phase_mode & 1)
                                phi = wrap_phase(phi + w_hat * t_hat);
                        }
                    }
                    phases.push_back(phi);
                }

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

                        if (pitch_shift_semi != 0) {
                            // Under shift, freq errors are INTEGRATED into
                            // propagated phase: parabolic bias of ±0.2 Hz on
                            // high partials made adjacent harmonics' relative
                            // phases slide at sub-Hz rates — an audible slow
                            // volume wobble on steady tones (300 saw down5,
                            // 15% envelope swing with every partial's own
                            // amplitude flat). The PV inst-freq is far more
                            // accurate for matched steady partials; trust it
                            // across the band. Unity keeps the original blend
                            // (phase is re-locked each frame there, so freq
                            // bias is harmless and the blend guards noisier
                            // inst estimates on transients).
                            freqs_hz[i] = (parabolic_freq_hz < 150.0)
                                ? 0.95 * inst_freq + 0.05 * parabolic_freq_hz
                                : 0.9 * inst_freq + 0.1 * parabolic_freq_hz;
                        } else if (parabolic_freq_hz < 150.0) {
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
                    // Reassigned frequency (bit 1): a chirping partial's true
                    // instantaneous frequency at the frame, replacing the PV/parabolic
                    // blend. Only where reassignment produced a valid estimate.
                    if ((phase_mode & 2) && reassigned_hz[i] > 0.0 &&
                        fabs(reassigned_hz[i] - parabolic_freq_hz) < 30.0)
                        freqs_hz[i] = reassigned_hz[i];
                }

                // JOINT AMP/PHASE: re-solve amplitude and phase for ALL of this
                // frame's partials together, at the frequencies just settled above
                // (so the estimation basis matches what will be synthesised).
                // Overwrites mags[]/phases[] in place; everything downstream —
                // tracking, the amplitude EMA, the shift path — is untouched.
                if (joint_mode != 0 && !analysis_windowed.empty())
                    joint_amp_phase(analysis_windowed, ANALYSIS_SIZE, sr, freqs_hz,
                                    mags, phases, joint_band_bins, joint_reg, joint_iters);

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
                    {   // only in-band candidates compete in the quality filter
                        vector<int> tmp;
                        for (int b : lf_peaks) {
                            double fb = b * (double)sr / (double)LF_ANALYSIS_SIZE;
                            if (fb >= lf_min_hz && fb < lf_cutoff_hz)
                                tmp.push_back(b);
                        }
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
                            // samples before the 4096 one (same frame center);
                            // minus the fractional-bin phase residue (see the
                            // 4096 path comment).
                            double ph = wrap_phase(lf_phase[kb] -
                                M_PI * (lf_freqs[i] - kb) +
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
                        { active_peaks[match_idx].fall_streak = (m_db < active_peaks[match_idx].current_db) ? active_peaks[match_idx].fall_streak + 1 : 0; double _a = (m_db > active_peaks[match_idx].current_db) ? amp_smooth_attack : ((m_db < amp_release_fast_drop * active_peaks[match_idx].current_db || (pitch_shift_semi != 0 && active_peaks[match_idx].fall_streak >= shift_release_streak_frames)) ? amp_release_fast : amp_smooth_release); active_peaks[match_idx].current_db = _a * m_db + (1.0 - _a) * active_peaks[match_idx].current_db; }
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
                // Level from THIS short frame, shape/frequency from the long
                // window (see transient_short_amp). Both spectra are Parseval-
                // scaled by their own window, so the ratio of their total
                // magnitude energy is the level correction.
                // Only at a REAL onset. The detector fires on spectral change with
                // no level rise -- every "transient" on Fairlight C3 and the Female
                // is a false positive by the level test -- and correcting level
                // there perturbs a sustained note for no reason. Require a >=6 dB
                // rise across this frame's own span in the source.
                if (transient_short_amp) {
                    double e_short = 0.0, e_long = 0.0;
                    for (int k = 0; k < (int)mag_spec.size(); k++) e_short += mag_spec[k] * mag_spec[k];
                    for (int k = 0; k < (int)short_unmatched_long_mag.size(); k++)
                        e_long += short_unmatched_long_mag[k] * short_unmatched_long_mag[k];
                    // per-sample energy density, so the different window lengths compare
                    double d_short = e_short / (double)mag_spec.size();
                    double d_long  = e_long  / (double)short_unmatched_long_mag.size();
                    if (d_long > 1e-20 && d_short >= 0.0) {
                        double g = sqrt(d_short / d_long);
                        // ATTENUATE ONLY. Pre-echo is the frame being handed more
                        // level than it holds (g < 1); that is the whole defect.
                        // Boosting (g > 1) would act on the frames AFTER a hit and,
                        // more to the point, on tonal files: the detector fires on
                        // spectral change with no level rise (every "transient" on
                        // Fairlight C3 and the Female is a false positive by the
                        // level test), and there g wanders either side of 1. Capping
                        // at 1 makes those frames a no-op and keeps the correction
                        // where it belongs.
                        if (getenv("TSA_DEBUG"))
                            fprintf(stderr, "tsa frame %d g=%.3f (%.1f dB)\n", frame_idx, g, 20*log10(g+1e-12));
                        // Attenuate only. Pre-echo is a frame handed more level
                        // than it holds; g > 1 would be the engine inventing level
                        // from a short window that cannot resolve a low partial.
                        // (Verified bit-identical to an uncapped +12 dB version on
                        // the whole corpus -- g > 1 never occurred in practice.)
                        if (g > 1.0) g = 1.0;
                        if (g < 0.0625) g = 0.0625;    // never cut more than -24 dB
                        for (int k = 0; k < (int)short_unmatched_long_mag.size(); k++)
                            short_unmatched_long_mag[k] *= g;
                    }
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
                        // fractional-bin phase residue corrected, as in the
                        // 4096 path
                        double ph    = wrap_phase(long_phase_spec[p_bin] -
                                                  M_PI * (freqs[i] - p_bin));

                        int match_idx = peak_to_track[i];
                        if (match_idx != -1) {
                            active_peaks[match_idx].freq_hz    = f_hz;
                            //temporal smoothing on current_db, same as long path.
                            { active_peaks[match_idx].fall_streak = (m_db < active_peaks[match_idx].current_db) ? active_peaks[match_idx].fall_streak + 1 : 0; double _a = (m_db > active_peaks[match_idx].current_db) ? amp_smooth_attack : ((m_db < amp_release_fast_drop * active_peaks[match_idx].current_db || (pitch_shift_semi != 0 && active_peaks[match_idx].fall_streak >= shift_release_streak_frames)) ? amp_release_fast : amp_smooth_release); active_peaks[match_idx].current_db = _a * m_db + (1.0 - _a) * active_peaks[match_idx].current_db; }
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
                        { ap.fall_streak = (true_mag < ap.current_db) ? ap.fall_streak + 1 : 0; double _a = (true_mag > ap.current_db) ? amp_smooth_attack : ((true_mag < amp_release_fast_drop * ap.current_db || (pitch_shift_semi != 0 && ap.fall_streak >= shift_release_streak_frames)) ? amp_release_fast : amp_smooth_release); ap.current_db = _a * true_mag + (1.0 - _a) * ap.current_db; }
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
                        { ap.fall_streak = (true_mag < ap.current_db) ? ap.fall_streak + 1 : 0; double _a = (true_mag > ap.current_db) ? amp_smooth_attack : ((true_mag < amp_release_fast_drop * ap.current_db || (pitch_shift_semi != 0 && ap.fall_streak >= shift_release_streak_frames)) ? amp_release_fast : amp_smooth_release); ap.current_db = _a * true_mag + (1.0 - _a) * ap.current_db; }
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

                        { ap.fall_streak = (true_mag < ap.current_db) ? ap.fall_streak + 1 : 0; double _a = (true_mag > ap.current_db) ? amp_smooth_attack : ((true_mag < amp_release_fast_drop * ap.current_db || (pitch_shift_semi != 0 && ap.fall_streak >= shift_release_streak_frames)) ? amp_release_fast : amp_smooth_release); ap.current_db = _a * true_mag + (1.0 - _a) * ap.current_db; }
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
            if (!ap.confirmed) {
                // File-start credit (see file_start_confirm_frames): hold the
                // unconfirmed newborn's measurement; it is added to this
                // frame's render list after analysis if the track confirms.
                if (frame_idx < file_start_confirm_frames)
                    file_start_pending.push_back({frame_idx, ap});
                continue;
            }
            frame_info.push_back(ap);
            ever_confirmed_ids.insert(ap.id);
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
    
    
    // File-start credit: release held-back frame-0 measurements of tracks that
    // went on to confirm. Exact lookahead -- analysis is complete here.
    {
        int released = 0;
        for (auto &pr : file_start_pending)
            if (ever_confirmed_ids.count(pr.second.id)) {
                frames_peaks[pr.first].push_back(pr.second);
                released++;
            }
        if (getenv("FILE_START_DEBUG"))
            fprintf(stderr, "file-start credit: %d of %zu held-back births released\n",
                    released, file_start_pending.size());
    }

    // Per-track frequency slope (see chirp_mode). Analysis is complete, so each
    // track's whole trajectory is known: take a centred difference over the
    // frames where the track is present, in Hz per second. A plain centred
    // difference is enough -- the solve still gains 11 dB with the slope 50%
    // wrong (docs 4d.9) -- so no higher-order estimator is warranted.
    if (chirp_mode) {
        unordered_map<int, vector<pair<int,int>>> where;   // track id -> (frame, index)
        for (int fi = 0; fi < (int)frames_peaks.size(); fi++)
            for (int k = 0; k < (int)frames_peaks[fi].size(); k++)
                where[frames_peaks[fi][k].id].push_back({fi, k});
        int chirped = 0, total = 0;
        for (auto &kv : where) {
            vector<pair<int,int>> &occ = kv.second;
            sort(occ.begin(), occ.end());
            for (size_t o = 0; o < occ.size(); o++) {
                size_t a = (o > 0) ? o - 1 : o;
                size_t b = (o + 1 < occ.size()) ? o + 1 : o;
                if (a == b) continue;                       // single sighting: stationary
                int fa = occ[a].first, fb = occ[b].first;
                if (fa >= (int)containsSynthPlacement.size() ||
                    fb >= (int)containsSynthPlacement.size()) continue;
                double dt = (double)(containsSynthPlacement[fb].start -
                                     containsSynthPlacement[fa].start) / (double)sr;
                if (!(dt > 0.0)) continue;
                double df = frames_peaks[fb].at(occ[b].second).freq_hz -
                            frames_peaks[fa].at(occ[a].second).freq_hz;
                double slope = df / dt;
                total++;
                if (fabs(slope) >= chirp_min_slope_hz_s) {
                    frames_peaks[occ[o].first].at(occ[o].second).freq_slope_hz_s = slope;
                    chirped++;
                }
            }
        }
        if (getenv("CHIRP_DEBUG"))
            fprintf(stderr, "chirp: %d of %d track-frames exceed %.0f Hz/s\n",
                    chirped, total, chirp_min_slope_hz_s);
    }

    //==========================================================================
    // JOINT AMP/PHASE, POST-TRACKING (joint_mode == 2)
    //==========================================================================
    // joint_mode 1 solves inside the analysis block, at the frequencies measured
    // there. But tracking then MOVES those frequencies before synthesis — the
    // median-of-3 de-jitter, the sub-200 Hz smoothing, LF replacement — and an
    // amplitude/phase fitted at frequency f is wrong when rendered at f'. Over a
    // 2048-sample frame a 5 Hz discrepancy is ~77 degrees of phase, which alone
    // caps the achievable SRR in the low single digits.
    //
    // Mode 2 therefore re-solves once more at the END, using each frame's FINAL
    // track frequencies — the ones actually about to be rendered — so the basis
    // the amplitudes and phases are fitted to is the basis that gets synthesised.
    // Only long frames carry the 4096 analysis, and only peaks still on the 4096
    // phase reference can be mixed into one solve, so the pass is limited to those.
    // Joint solve under PITCH SHIFT (see joint_shift_mode). At unity the engine
    // re-derives each partial's phase from the analysis every frame, so the jointly
    // fitted (amplitude, phase) PAIR is rendered as fitted and the solve's whole
    // premise holds. Under shift the engine propagates each track's phase
    // independently instead, discarding those relative phases -- and the joint
    // amplitudes were chosen precisely to account for how the partials interfere at
    // the fitted phases. Rendered at different relative phases they are wrong:
    // partials the solve boosted to offset cancellation now add constructively.
    // Measured with the perceptual metric (docs 5c/6b, shifted NMR, lower better):
    // DrumLoop up5 12.9 -> 5.4, Female 9.3 -> -0.7, Happy 8.2 -> 0.4 when the solve
    // is skipped under shift; down-shift gains 5.5-6.9 dB too. This is the
    // "estimate on the basis you synthesise" law of 4c.1, which was only ever
    // checked at unity.
    //   0 = use the joint solve under shift as well (the shipped behaviour)
    //   1 = skip it under shift, keeping the per-peak estimate there
    // Env JOINT_SHIFT.
    int joint_shift_mode = 0;
    if (const char* e = getenv("JOINT_SHIFT")) joint_shift_mode = atoi(e);
    const bool joint_here = (joint_mode == 2) &&
                            !(joint_shift_mode == 1 && pitch_shift_semi != 0);
    if (joint_here) {
        vector<double> seg(ANALYSIS_SIZE);
        // The solve runs AFTER tracking, so it bypasses the amplitude EMA that the
        // tracker applies. Removing the interference bias without replacing that
        // variance control made every per-partial amplitude a raw per-frame
        // measurement again — visible as a jump in trajectory_jitter and env_p2p on
        // 11 files. Re-apply the same asymmetric EMA to the solved amplitudes here,
        // per track id, in frame order. Phase is deliberately NOT smoothed: at unity
        // the engine already re-derives it per frame, so joint phase is no worse in
        // kind, and smoothing it would undo the gain.
        unordered_map<int, double> joint_amp_prev;
        unordered_map<int, pair<double,double>> joint_amp_hist;   // median-of-3 history
        unordered_map<int, int> joint_fall_streak;
        for (int fi = 0; fi < (int)frames_peaks.size() && fi < (int)containsSynthPlacement.size(); fi++) {
            const SynthInformation &si = containsSynthPlacement[fi];
            if (si.size != LONG_SIZE) continue;
            vector<PeakTrack> &fp = frames_peaks[fi];
            vector<int> idx;
            vector<double> f_hz, m_fft, ph, slope;
            for (int i = 0; i < (int)fp.size(); i++) {
                if (fp[i].analysis_fft_size != ANALYSIS_SIZE) continue;
                if (!(fp[i].freq_hz > 0.0) || fp[i].freq_hz >= 0.5 * sr) continue;
                idx.push_back(i);
                f_hz.push_back(fp[i].freq_hz);
                m_fft.push_back(fp[i].current_db);
                ph.push_back(fp[i].phase);
                slope.push_back(chirp_mode ? fp[i].freq_slope_hz_s : 0.0);
            }
            if (idx.empty()) continue;

            int analysis_start = si.start + si.size / 2 - ANALYSIS_SIZE / 2;
            for (int i = 0; i < ANALYSIS_SIZE; i++) {
                int s = analysis_start + i;
                seg[i] = (s >= 0 && s < (int)singleChannelData.size())
                             ? (double)singleChannelData[s] * (double)analysis_hanning[i]
                             : 0.0;
            }
            joint_amp_phase(seg, ANALYSIS_SIZE, sr, f_hz, m_fft, ph,
                            joint_band_bins, joint_reg, joint_iters, 48,
                            chirp_mode ? &slope : nullptr,
                            chirp_mode ? &analysis_hanning : nullptr);
            double frame_max = 0.0;
            if (joint_smooth == 3)
                for (int k = 0; k < (int)idx.size(); k++)
                    if (m_fft[k] > frame_max) frame_max = m_fft[k];
            const double raw_floor = frame_max * pow(10.0, -joint_smooth_raw_db / 20.0);
            for (int k = 0; k < (int)idx.size(); k++) {
                PeakTrack &pk = fp[idx[k]];
                double meas = m_fft[k];
                const bool pass_raw = (joint_smooth == 3 && meas >= raw_floor);
                if (joint_smooth == 2) {
                    // Median-of-3 over the raw solved amplitudes, per track, in
                    // frame order. History holds the two previous RAW values so
                    // the filter never feeds on its own output.
                    int id = pk.id;
                    auto it = joint_amp_hist.find(id);
                    if (it == joint_amp_hist.end()) {
                        joint_amp_hist[id] = {meas, meas};   // first sighting: pass through
                    } else {
                        double a1 = it->second.first, a2 = it->second.second;
                        double lo = min(min(meas, a1), a2), hi = max(max(meas, a1), a2);
                        double med = meas + a1 + a2 - lo - hi;
                        it->second = {meas, a1};
                        meas = med;
                    }
                } else if (joint_smooth && !pass_raw) {
                    int id = pk.id;
                    auto it = joint_amp_prev.find(id);
                    if (it == joint_amp_prev.end()) {
                        joint_amp_prev[id] = meas;      // first sighting: take it as-is
                        joint_fall_streak[id] = 0;
                    } else {
                        double prev = it->second;
                        int &streak = joint_fall_streak[id];
                        streak = (meas < prev) ? streak + 1 : 0;
                        double a = (meas > prev)
                            ? amp_smooth_attack
                            : ((meas < amp_release_fast_drop * prev ||
                                (pitch_shift_semi != 0 && streak >= shift_release_streak_frames))
                                   ? amp_release_fast : amp_smooth_release);
                        meas = a * meas + (1.0 - a) * prev;
                        it->second = meas;
                    }
                }
                pk.current_db = meas;
                pk.phase = ph[k];
            }
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
    if ((pitch_shift_semi != 0 || unity_dedup) && shift_dedup_bins > 0.0) {
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
    // Steady-tone frequency lock (shift mode): even accurate per-frame
    // estimates wiggle a few mHz..0.05 Hz, and propagated phase INTEGRATES
    // the wiggle — adjacent harmonics' relative phases random-walk, so the
    // waveform crest slowly swells and dips (audible slow "volume wobble"
    // on exposed steady tones; every partial's own amplitude is flat).
    // If a track's estimate stays within shift_lock_tol (rel) of its running
    // mean for shift_lock_frames frames, synthesize at the frozen mean:
    // phase then advances perfectly linearly. Vibrato (10-25 cents/frame)
    // blows the tolerance and never locks.
    unordered_map<int, double> shift_lock_mean;
    unordered_map<int, int> shift_lock_count;
    const double shift_lock_tol = 0.0015;  // ~2.6 cents
    const int shift_lock_frames = 8;
    // Harmonic phase coherence: per-member link to its stack root. The offset
    // (member phase minus k*root phase) is captured ONCE when the link forms and
    // then held rigid across frames, so the stack keeps the waveform shape it was
    // measured with instead of drifting out of it.
    struct HarmLink { int root_id = -1; int k = 0; double offset = 0.0; bool captured = false;
                      double theta = 0.0; bool theta_init = false; };
    unordered_map<int, HarmLink> harm_link;

    // Oscillator-bank (synth_mode==1) per-track node lists, captured in the frame
    // loop below and rendered after it by mq_synthesize. The frame loop reuses all
    // of its existing per-frame frequency conditioning (dedup, median, EMA smooth,
    // steady-tone lock) and phase logic; MQ mode just records a node instead of
    // painting a windowed OLA frame. mq_nodes = the (possibly shifted) main model;
    // mq_nodes_unity = the parallel unshifted model for the shift-mode residual.
    unordered_map<int, vector<MQNode>> mq_nodes;
    unordered_map<int, vector<MQNode>> mq_nodes_unity;


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

        // PCOUNT_DEBUG: partials actually rendered per synthesis frame. Env-gated
        // like LF_DEBUG/FTRAJ_DEBUG; completely inert unless PCOUNT_DEBUG is set.
        // Used by the Stage-0 ceiling probe to compare the engine's partial budget
        // against the oracle fit's K on the same material.
        if (getenv("PCOUNT_DEBUG")) {
            double mx = 0.0;
            for (auto &p : frames_peaks[frame_idx]) {
                double a = 4.0 * p.current_db / (double)p.analysis_fft_size;
                if (fabs(a) > mx) mx = fabs(a);
            }
            int n40 = 0, n60 = 0;
            for (auto &p : frames_peaks[frame_idx]) {
                double a = fabs(4.0 * p.current_db / (double)p.analysis_fft_size);
                if (mx > 0 && a >= mx * 0.01)   n40++;   // within 40 dB of loudest
                if (mx > 0 && a >= mx * 0.001)  n60++;   // within 60 dB
            }
            fprintf(stderr, "PCOUNT %d %d %d %d %d\n", frame_idx, frame_size,
                    (int)frames_peaks[frame_idx].size(), n40, n60);
        }
        // PFREQ_DEBUG: dump the per-frame track frequencies + linear amplitudes the
        // engine is about to render, so the ceiling probe can hold the engine's own
        // frequencies fixed and fit only amplitude/phase optimally. Env-gated, inert
        // when unset.
        if (getenv("PFREQ_DEBUG")) {
            fprintf(stderr, "PFRAME %d %d %d %d\n", frame_idx,
                    current_information.start, frame_size,
                    (int)frames_peaks[frame_idx].size());
            for (auto &p : frames_peaks[frame_idx])
                fprintf(stderr, "PF %.4f %.8f %.6f\n", p.freq_hz,
                        4.0 * p.current_db / (double)p.analysis_fft_size, p.phase);
        }

        // Shift-mode: ids alive in THIS frame (rebirth candidates must be dead).
        unordered_set<int> cur_frame_ids;
        if (pitch_shift_semi != 0)
            for (auto &p : frames_peaks[frame_idx]) cur_frame_ids.insert(p.id);

        // HARMONIC PHASE COHERENCE — group this frame's tracks into stacks.
        // Ascending-frequency greedy: a low root claims unassigned tracks within
        // harmonic_lock_tol of k*f_root (k >= 2). Members then derive frequency
        // (exactly k*f_root) and phase (k*root_phase + a captured true offset)
        // from the root, so the stack's relative phases cannot drift apart.
        // k*wrap(theta) is congruent to k*theta mod 2pi for integer k, so wrapped
        // phase bookkeeping stays valid.
        //
        // ordv stays the identity order unless the lock is actually active, because
        // changing the accumulation order of frame_signal changes float rounding —
        // unity must stay byte-identical.
        unordered_map<int, pair<int,int>> member_of;      // member id -> (root id, k)
        struct RootRender { double phase0; double shiftedFreq; double analysis_phase; };
        unordered_map<int, RootRender> root_render;
        vector<int> ordv(frames_peaks[frame_idx].size());
        for (size_t oi = 0; oi < ordv.size(); oi++) ordv[oi] = (int)oi;
        if (pitch_shift_semi != 0 && shift_harmonic_lock) {
            auto &fp = frames_peaks[frame_idx];
            sort(ordv.begin(), ordv.end(), [&fp](int a, int b) {
                return fp[a].freq_hz < fp[b].freq_hz;
            });
            // ROOT VALIDITY. The root must be a plausible FUNDAMENTAL, not merely the
            // lowest-numbered track. Without this the grouping picked whatever sat at
            // the bottom of the list: on 440sawtooth down5 that was sub-audio LF junk
            // at 21.7 Hz, ranked 463rd by amplitude and 56 dB below the loudest
            // partial, and 245-460 real partials were force-quantised onto its 21.7 Hz
            // grid. (Contrast 300hzSaw, whose root is the true 300 Hz fundamental,
            // amplitude rank 1 — which is why the lock worked perfectly there.)
            // The audible result was a lurch at t=0.17 s the moment the steady gate
            // matured, with crest factor jumping 2.24 -> 2.69 into a pinned limiter
            // and stack membership churning 245 -> 127 -> 447 across three frames.
            double frame_max_amp = 0.0;
            for (auto &pk : fp) {
                double A = fabs(4.0 * pk.current_db / (double)pk.analysis_fft_size);
                if (A > frame_max_amp) frame_max_amp = A;
            }
            double root_amp_floor = frame_max_amp *
                                    pow(10.0, -harmonic_root_min_rel_db / 20.0);
            double frame_energy = 0.0;
            for (auto &pk : fp) {
                double A = 4.0 * pk.current_db / (double)pk.analysis_fft_size;
                frame_energy += A * A;
            }

            vector<char> taken(fp.size(), 0);
            for (size_t ii = 0; ii < ordv.size(); ii++) {
                int i = ordv[ii];
                if (taken[i]) continue;
                double fr = fp[i].freq_hz;
                if (fr <= 0.0) continue;
                if (fr > harmonic_root_max_hz) break;   // ascending: no roots left
                if (fr < harmonic_root_min_hz) continue;        // sub-audio / rumble
                double root_amp = fabs(4.0 * fp[i].current_db / (double)fp[i].analysis_fft_size);
                if (root_amp < root_amp_floor) continue;        // too weak to be a fundamental
                // STEADY-ROOT GATE. Members are rendered at exactly k*f_root, so any
                // error in the root's frequency is multiplied by k — on a vibrato
                // source, where the root estimate moves every frame, that is worse
                // than leaving the stack independent (measured on the 440 vibrato
                // saw: shape_consistency 0.914 unlocked -> 0.865..0.885 locked, at
                // every tolerance). Only lock stacks whose root has already passed
                // the steady-tone frequency lock, which a vibrato track never does.
                // shift_lock_count is carried over from the previous frame; use
                // find() so probing it cannot insert a zero entry.
                auto itlk = shift_lock_count.find(fp[i].id);
                if (harmonic_lock_warmup_frames > 0 &&
                    (itlk == shift_lock_count.end() ||
                     itlk->second < harmonic_lock_warmup_frames))
                    continue;
                // Collect candidates first; only commit if the stack looks like a real
                // harmonic series (see the low-k support test below).
                vector<pair<int,int>> cand;             // (index into fp, k)
                bool low_k = false;
                for (size_t jj = ii + 1; jj < ordv.size(); jj++) {
                    int j = ordv[jj];
                    if (taken[j]) continue;
                    double r = fp[j].freq_hz / fr;
                    int hk = (int)floor(r + 0.5);
                    if (hk < 2) continue;
                    // Relative tolerance, ADDITIONALLY capped in absolute Hz. The
                    // relative test alone lets a high partial be dragged a long way:
                    // 0.05% of 20 kHz is 10 Hz, and measured shifts reached 11 Hz.
                    // A genuine harmonic sits far closer than this to k*f0.
                    double dev_hz = fabs(fp[j].freq_hz - (double)hk * fr);
                    if (harmonic_rps_mode) {
                        // RPS absorbs the deviation in theta_k each frame, so only a
                        // relative test is needed and it can be loose.
                        if (dev_hz >= harmonic_lock_tol_rps * fp[j].freq_hz) continue;
                    } else {
                        if (dev_hz >= harmonic_lock_tol * fp[j].freq_hz) continue;
                        if (dev_hz >= harmonic_lock_max_dev_hz) continue;
                    }
                    cand.push_back({j, hk});
                    if (hk <= 3) low_k = true;
                }
                // HARMONIC SUPPORT. A real series has its low harmonics present. A
                // spurious low root instead produces a dense grid that catches
                // unrelated partials by chance at high k, which is exactly how the
                // 21.7 Hz root captured hundreds of them. Requiring k=2 or k=3 to be
                // there costs nothing on a genuine stack.
                // Energy share of this candidate stack within the frame.
                double stack_e = root_amp * root_amp;
                for (auto &c : cand) {
                    double A = 4.0 * fp[c.first].current_db / (double)fp[c.first].analysis_fft_size;
                    stack_e += A * A;
                }
                bool enough_energy = (frame_energy > 0.0) &&
                                     (stack_e >= harmonic_stack_min_energy_frac * frame_energy);

                bool any = false;
                if (low_k && enough_energy && (int)cand.size() >= harmonic_min_members) {
                    for (auto &c : cand) {
                        member_of[fp[c.first].id] = {fp[i].id, c.second};
                        taken[c.first] = 1;
                        any = true;
                    }
                }
                if (any) taken[i] = 1;                  // a root cannot join another stack
            }
        }

        // Ascending order guarantees a stack's root is rendered before its members.
        for (int p_idx : ordv) {
            auto &peak = frames_peaks[frame_idx][p_idx];
            double freq = peak.freq_hz;
            double mag = 4.0 * peak.current_db / (double)peak.analysis_fft_size;
            double phase = peak.phase;

            // Parallel unshifted tonal model for the shift-mode residual:
            // exactly the unity branch (raw freq, analysis-derived phase).
            if (want_unity_model) {
                double inc_u = 2.0 * M_PI * freq / (double)sr;
                int off_u = peak.analysis_fft_size / 2 - frame_size / 2;
                double ph0_u = wrap_phase(phase + inc_u * off_u);
                if (synth_mode == 1) {
                    // Record a unity node at the frame centre (phase carried to
                    // the centre); mq_synthesize renders it with the match-phase rule.
                    double phc = wrap_phase(ph0_u + inc_u * (frame_size / 2.0));
                    mq_nodes_unity[peak.id].push_back(
                        { (double)current_information.start + frame_size / 2.0, freq, mag, phc });
                } else if (is_long) {
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

                // Steady-tone frequency lock (see declaration comment).
                double &lm = shift_lock_mean[peak.id];
                if (lm <= 0.0) lm = freq;
                if (fabs(freq - lm) < shift_lock_tol * lm) {
                    lm += 0.05 * (freq - lm);  // slow mean update
                    if (++shift_lock_count[peak.id] >= shift_lock_frames)
                        freq = lm;
                } else {
                    lm = freq;
                    shift_lock_count[peak.id] = 0;
                }
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

                    // Members override frequency and phase from their root's
                    // just-rendered trajectory. The relative offset is captured
                    // ONCE from the true propagated phases and then held rigid,
                    // so the stack keeps one waveform shape for its lifetime.
                    auto mo = member_of.find(peak.id);
                    if (mo != member_of.end()) {
                        auto rr = root_render.find(mo->second.first);
                        if (rr != root_render.end()) {
                            double hk = (double)mo->second.second;
                            double sf_h = rr->second.shiftedFreq * hk;
                            if (harmonic_rps_mode &&
                                (harmonic_rps_max_k <= 0 || mo->second.second <= harmonic_rps_max_k)) {
                                // theta_k = phi_k - k*phi_root, from THIS frame's analysis
                                // phases (both measured at the same instant), smoothed a
                                // little so a noisy phase read cannot jolt the shape.
                                HarmLink &L = harm_link[peak.id];
                                double th = wrap_phase(phase - hk * rr->second.analysis_phase);
                                if (L.root_id != mo->second.first || L.k != mo->second.second ||
                                    !L.theta_init) {
                                    L.root_id = mo->second.first; L.k = mo->second.second;
                                    L.theta = th; L.theta_init = true; L.captured = true;
                                } else {
                                    // wrapped EMA: move along the shorter arc
                                    double d = wrap_phase(th - L.theta);
                                    L.theta = wrap_phase(L.theta + 0.5 * d);
                                }
                                // Re-anchor the PHASE to the root each frame; leave the
                                // member's own frequency alone (no harmonic quantisation).
                                phase0 = wrap_phase(hk * rr->second.phase0 + L.theta);
                            } else if (sf_h > 0.0 && sf_h < nyquist) {
                                HarmLink &L = harm_link[peak.id];
                                if (L.root_id != mo->second.first ||
                                    L.k != mo->second.second || !L.captured) {
                                    L.root_id = mo->second.first;
                                    L.k = mo->second.second;
                                    L.offset = wrap_phase(phase0 - hk * rr->second.phase0);
                                    L.captured = true;
                                }
                                phase0 = wrap_phase(hk * rr->second.phase0 + L.offset);
                                shiftedFreq = sf_h;
                                phase_inc = 2.0 * M_PI * sf_h / (double)sr;
                                // re-derive the anti-alias gain at the new frequency
                                mag_syn = mag;
                                double aa_lo2 = 0.9 * nyquist;
                                if (sf_h > aa_lo2)
                                    mag_syn = mag * 0.5 * (1.0 + cos(M_PI * (sf_h - aa_lo2) / (nyquist - aa_lo2)));
                            }
                        }
                    }
                }
                if (getenv("HLOCK_STATS")) {
                    static long long locked = 0, total = 0, lastf = -1;
                    if (frame_idx != lastf && lastf >= 0 && frame_idx % 50 == 0)
                        fprintf(stderr, "hlock: %lld of %lld track-frames locked (%.1f%%)\n",
                                locked, total, total ? 100.0 * locked / total : 0.0);
                    lastf = frame_idx; total++;
                    if (harm_link.count(peak.id) && harm_link[peak.id].captured) locked++;
                }
                root_render[peak.id] = {phase0, shiftedFreq, phase};

                if (synth_mode == 1) {
                    // Record a node at the frame centre; mq_synthesize renders the
                    // whole track continuously after the loop. Unity uses the
                    // match-phase rule (shape); shift uses phase propagation.
                    double phc = wrap_phase(phase0 + phase_inc * (frame_size / 2.0));
                    mq_nodes[peak.id].push_back(
                        { (double)current_information.start + frame_size / 2.0,
                          shiftedFreq, mag_syn, phc });
                    continue; // skip OLA render + phase propagation for this frame
                }

                // Linear-FM term (see chirp_mode). tau is measured from the frame
                // CENTRE, where phase_inc's frequency is the measured one, so the
                // quadratic vanishes there and the frame still agrees with the
                // analysis phase. Under shift the rate scales with the frequency.
                double chirp_c = 0.0;
                if (chirp_mode && peak.freq_slope_hz_s != 0.0) {
                    double slope = peak.freq_slope_hz_s * shiftFactor;
                    chirp_c = M_PI * slope / ((double)sr * (double)sr);   // 2pi * slope/2
                }
                const double half = 0.5 * (double)frame_size;
                for (int n = 0; n < frame_size; n++) {
                    double sample_phase = phase0 + phase_inc * n;
                    if (chirp_c != 0.0) {
                        double d = (double)n - half;
                        sample_phase += chirp_c * d * d;
                    }
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
                    double adv = phase_inc * next_delta_samples;
                    if (chirp_c != 0.0) {           // same quadratic, evaluated at the hop
                        double d1 = (double)next_delta_samples - half, d0 = -half;
                        adv += chirp_c * (d1 * d1 - d0 * d0);
                    }
                    synth_phase_by_track[peak.id] = wrap_phase(phase0 + adv);
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

        if (synth_mode == 1) continue; // MQ renders after the loop; no windowing/OLA

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
    
    if (synth_mode == 1) {
        // Oscillator-bank render: one continuous oscillator per track, straight
        // into synthesized_signal (no window_sum normalization — there is no
        // overlap to normalize). Unity pins the waveform shape via measured
        // phase; shift propagates phase. Birth/death ramp = one long hop.
        int mq_fade = LONG_SIZE / 2;
        mq_synthesize(mq_nodes, synthesized_signal,
                      /*match_phase=*/ (pitch_shift_semi == 0), sr, nyquist, mq_fade);
        if (want_unity_model)
            mq_synthesize(mq_nodes_unity, synth_unity,
                          /*match_phase=*/ true, sr, nyquist, mq_fade);
    } else {
        for (size_t i = 0; i < synthesized_signal.size(); i++) {
            if (window_sum[i] > 1.0e-8f) {
                synthesized_signal[i] /= window_sum[i];
                if (want_unity_model) synth_unity[i] /= window_sum[i];
            }
        }
    }
    // Bridge overlap-add coverage notches at the short->long window switch. The
    // short_to_long transition window has a 448-sample zero prefix and the last
    // short frame's Hann tapers to zero at the junction, leaving ~2 samples with
    // window_sum ~ 0 (a COLA hole) -> a 2-sample dropout to zero = an impulsive
    // click on every transient. Linearly interpolate any run of near-zero-coverage
    // samples from the nearest well-covered neighbours (normal audio is untouched).
    // OLA path only: MQ writes no window_sum, so there are no coverage holes.
    if (synth_mode == 0) {
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

    // File-start level cap (see file_start_cap_ms). Compare a short moving RMS of
    // the output against the same measure on the input over the opening window and
    // scale the output down wherever it exceeds the input by more than the
    // headroom. A cap, never a boost: it can only remove energy the input does not
    // support, which is exactly the frame-0 pre-echo and nothing else.
    if (file_start_cap_ms > 0.0 && synth_mode == 0) {
        int CAP = (int)(file_start_cap_ms * 0.001 * sr);
        if (CAP > (int)total_length) CAP = (int)total_length;
        int EW = (int)(0.004 * sr);                       // 4 ms envelope
        if (EW < 16) EW = 16;
        const double lim = pow(10.0, file_start_cap_headroom_db / 20.0);
        if (CAP > EW) {
            vector<double> ci(CAP + 1, 0.0), co(CAP + 1, 0.0);
            for (int i = 0; i < CAP; i++) {
                double xi = (i < (int)singleChannelData.size()) ? (double)singleChannelData[i] : 0.0;
                double yi = synthesized_signal[i];
                ci[i + 1] = ci[i] + xi * xi;
                co[i + 1] = co[i] + yi * yi;
            }
            vector<float> g(CAP, 1.0f);
            for (int i = 0; i < CAP; i++) {
                int a = i - EW / 2, b = a + EW;
                if (a < 0) { a = 0; b = EW; }
                if (b > CAP) { b = CAP; a = CAP - EW; }
                double ei = (ci[b] - ci[a]) / (double)(b - a);
                double eo = (co[b] - co[a]) / (double)(b - a);
                double allowed = ei * lim * lim;
                g[i] = (eo > 1e-20 && eo > allowed) ? (float)sqrt(allowed / eo) : 1.0f;
            }
            // 2 ms box smoothing so the cap cannot buzz, then fade the cap out
            // over the last quarter of the window so it cannot leave a step.
            int SM = (int)(0.002 * sr); if (SM < 4) SM = 4;
            vector<double> cg(CAP + 1, 0.0);
            for (int i = 0; i < CAP; i++) cg[i + 1] = cg[i] + g[i];
            for (int i = 0; i < CAP; i++) {
                int a = i - SM / 2, b = a + SM;
                if (a < 0) { a = 0; b = SM; }
                if (b > CAP) { b = CAP; a = CAP - SM; }
                double gm = (cg[b] - cg[a]) / (double)(b - a);
                double fade = (i > 3 * CAP / 4)
                                  ? (double)(CAP - i) / (double)(CAP - 3 * CAP / 4) : 1.0;
                double gg = 1.0 + (gm - 1.0) * fade;
                synthesized_signal[i] = (float)(synthesized_signal[i] * gg);
                if (want_unity_model) synth_unity[i] = (float)(synth_unity[i] * gg);
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
    // Run when the stochastic fill is enabled, OR when true-phase residual is
    // requested at unity (mode 1 is self-contained: it adds the real residual at
    // residualScale gain and does not need the stochastic makeup gain > 0).
    bool run_true_phase = (residual_mode == 1 && pitch_shift_semi == 0 && settings.residualScale > 0.0);
    if ((residual_noise_gain > 0.0 || run_true_phase) && (pitch_shift_semi == 0 || want_unity_model)) {
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
            if (residual_mode == 1 && pitch_shift_semi == 0) {
                // True-phase residual: add the real model-subtracted spectrum
                // above hp_bin, preserving phase (temporal micro-structure ->
                // natural breath/consonants). fin/fsy carry the same real-FFT
                // packing (re at k, im at RN-k), so subtract componentwise.
                for (int k = hp_bin; k < RN/2; k++) {
                    fout[k]    = fin[k]    - fsy[k];
                    fout[RN-k] = fin[RN-k] - fsy[RN-k];
                }
                if (hp_bin <= RN/2) fout[RN/2] = fin[RN/2] - fsy[RN/2];
            } else {
                // Stochastic residual: white phase on the deficit magnitude.
                // With residual_shift_mode the deficit SPECTRUM is transposed by the
                // same ratio as the tonal content: destination bin j draws its
                // magnitude from source bin j/ratio (linearly interpolated, so no
                // comb gaps when the ratio stretches the spectrum), and amplitudes
                // are scaled by 1/sqrt(ratio) so total noise energy is preserved.
                const bool rshift = (residual_shift_mode == 1 && pitch_shift_semi != 0);
                const double rratio = rshift ? pow(2.0, pitch_shift_semi / 12.0) : 1.0;
                const double rgain_e = rshift ? 1.0 / sqrt(rratio) : 1.0;
                for (int k = 1; k < RN/2; k++) {
                    double mag;
                    if (rshift) {
                        double src = (double)k / rratio;
                        int s0 = (int)floor(src);
                        double fr = src - s0;
                        double a = (s0 >= 1 && s0 < RN/2) ? sm[s0] : 0.0;
                        double b = (s0 + 1 >= 1 && s0 + 1 < RN/2) ? sm[s0 + 1] : 0.0;
                        mag = (a + (b - a) * fr) * rgain_e;
                    } else {
                        mag = sm[k];
                    }
                    double th = 2.0 * M_PI * frand();
                    fout[k]    = (float)(mag * cos(th));
                    fout[RN-k] = (float)(mag * sin(th));
                }
                fout[0] = 0.0f; fout[RN/2] = 0.0f;   // no DC / Nyquist noise
            }
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

        // True-phase mode restores the real residual, so add at unity gain
        // (× residualScale) rather than the stochastic fill's makeup gain.
        double rgain = (residual_mode == 1 && pitch_shift_semi == 0)
                           ? settings.residualScale : residual_noise_gain;

        // Envelope cap (see residual_env_cap). Compare a short moving-RMS of the
        // fill against the same measure on the input and scale the fill down
        // wherever it exceeds it. This is a cap, never a boost, so it can only
        // remove energy the input does not support -- notably the pre-echo the
        // 21 ms residual window spreads ahead of an onset.
        vector<float> env_gain;
        if (residual_env_cap && RL > 0) {
            int EW = (int)(residual_env_ms * 0.001 * sr);
            if (EW < 16) EW = 16;
            if (EW > RL) EW = RL;
            // centred moving mean-square, running sums
            vector<double> ci(RL + 1, 0.0), cr(RL + 1, 0.0);
            for (int i = 0; i < RL; i++) {
                double xi = (i < (int)singleChannelData.size()) ? (double)singleChannelData[i] : 0.0;
                double ri = (res_wsum[i] > wfloor ? res_signal[i] / res_wsum[i]
                                                  : res_signal[i] / wfloor);
                ci[i + 1] = ci[i] + xi * xi;
                cr[i + 1] = cr[i] + ri * ri;
            }
            env_gain.assign(RL, 1.0f);
            for (int i = 0; i < RL; i++) {
                int a = i - EW / 2, b = a + EW;
                if (a < 0) { a = 0; b = EW; }
                if (b > RL) { b = RL; a = RL - EW > 0 ? RL - EW : 0; }
                double ei = (ci[b] - ci[a]) / (double)(b - a);
                double er = (cr[b] - cr[a]) / (double)(b - a);
                // rgain scales the fill on the way in, so compare like with like
                double erg = er * rgain * rgain;
                env_gain[i] = (erg > 1e-20 && erg > ei)
                                  ? (float)sqrt(ei / erg) : 1.0f;
            }
            // 1 ms box smoothing of the gain so the cap cannot buzz
            int SM = (int)(0.001 * sr); if (SM < 4) SM = 4;
            vector<double> cg(RL + 1, 0.0);
            for (int i = 0; i < RL; i++) cg[i + 1] = cg[i] + env_gain[i];
            for (int i = 0; i < RL; i++) {
                int a = i - SM / 2, b = a + SM;
                if (a < 0) { a = 0; b = SM; }
                if (b > RL) { b = RL; a = RL - SM > 0 ? RL - SM : 0; }
                env_gain[i] = (float)((cg[b] - cg[a]) / (double)(b - a));
            }
        }
        if (wfloor > 0.0f) {
            for (int i = 0; i < RL; i++) {
                float denom = res_wsum[i] > wfloor ? res_wsum[i] : wfloor;
                double g = rgain * (env_gain.empty() ? 1.0f : env_gain[i]);
                synthesized_signal[i] += g * res_signal[i] / denom;
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
    
    
    // Pitch-synchronous: un-warp the warped-domain output back to real time so the
    // original vibrato is restored (at the shifted pitch, under a pitch shift).
    if (settings.pitchSync != 0 && !ps_tau.empty()) {
        vector<float> unwarped;
        pitchsync_unwarp(synthesized_signal, ps_tau, unwarped);
        synthesized_signal = std::move(unwarped);
        cout << "[pitch-sync] un-warped to " << synthesized_signal.size() << " samples\n";
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
    cout << "Usage: AdditiveSynthFreqMask <input wav> <output wav> <block size in samples> [pitch shift semitones]" << endl;
    cout << "        <input wav>  Path to the input wav file." << endl;
    cout << "        <output_wav>  Path to the output wav file." << endl;
    cout << "        <block size> block size in samples." << endl;
    cout << "        [pitch shift semitones]  Optional integer; 0 (default) = unity." << endl;
    cout << "        [synth mode]  Optional; 0 (default) = OLA, 1 = oscillator bank (MQ)." << endl;
    cout << "        [residual scale]  Optional; multiplies residual gain (1 default, 0 = off)." << endl;
    cout << "        [unity dedup]  Optional; 1 = run duplicate-partial dedup at unity too." << endl;
    cout << "        [residual mode]  Optional; 0 = stochastic (default), 1 = true-phase (unity)." << endl;
    cout << "        [joint mode]  Optional; 2 = joint least squares at final track freqs (default)," << endl;
    cout << "                      1 = joint solve in analysis (near no-op), 0 = legacy per-peak." << endl;
    cout << "        [residual hp hz]  Optional; override residual high-pass freq (<0 = keep default)." << endl;
    cout << "        [amp mode]  Optional; 0 = parabolic peak (default), 1 = energy-integrated (vibrato)." << endl;
    cout << "        [phase mode]  Optional; 0 = stationary (default), 1 = reassigned (chirp-aware)." << endl;
    cout << "        [pitch sync]  Optional; 1 = auto (default, self-engages on singing), 0 = off, 2 = force on." << endl;
}

bool parseArgs(int argc, const char* argv[], AppSettings& settings)
{
    // 4 args = legacy form (unity, OLA). Optional 5th arg = pitch shift in
    // semitones; optional 6th arg = synth mode (0 OLA, 1 oscillator bank), so a
    // battery can render every pitch condition and both engines from one binary
    // without recompiling. Defaults stay 0, so existing invocations are unchanged.
    if (argc < 4 || argc > 14)
    {
        return false;
    }

    settings.inputWavFilePath = argv[1];
    //settings.inputCSVPath = argv[2];
    settings.outputWavFilePath = argv[2];
    //settings.sampleRate = atof(argv[4]);
    settings.blockSize = atof(argv[3]);
    if (argc >= 5)
        settings.pitchShiftSemi = atoi(argv[4]);
    if (argc >= 6)
        settings.synthMode = atoi(argv[5]);
    if (argc >= 7)
        settings.residualScale = atof(argv[6]);
    if (argc >= 8)
        settings.unityDedup = atoi(argv[7]);
    if (argc >= 9)
        settings.residualMode = atoi(argv[8]);
    if (argc >= 10)
        settings.residualHpHz = atof(argv[9]);
    if (argc >= 11)
        settings.ampMode = atoi(argv[10]);
    if (argc >= 12)
        settings.phaseMode = atoi(argv[11]);
    if (argc >= 13)
        settings.pitchSync = atoi(argv[12]);
    if (argc > 13)
        settings.jointMode = atoi(argv[13]);

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




