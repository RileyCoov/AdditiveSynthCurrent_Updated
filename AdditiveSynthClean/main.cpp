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
vector<int> detect_peaks(vector<double>& magnitude, double threshold) {
    vector<int> peaks;
    for (size_t i = 1; i < magnitude.size()-1; i++) {
        if (magnitude[i] > magnitude[i-1] && magnitude[i] > magnitude[i+1] && magnitude[i] > threshold) {
            peaks.push_back((int) i);
        }
    }
    return peaks;
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
int find_best_match_peak(int peak_bin, const vector<PeakTrack>& active_peaks, double peak_tolerance=2.0) {
    int best_idx = -1;
    double best_diff = numeric_limits<double>::infinity();
    double diff = 0.0;
    for (size_t i = 0; i < active_peaks.size(); i++) {
        if (active_peaks[i].alive) {
            diff = abs((double)active_peaks[i].peak_bin - (double)peak_bin);
            if (diff < peak_tolerance && diff < best_diff) {
                best_diff = diff;
                best_idx = (int)i;
            }
        }
    }
    return best_idx;
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
        MagnitudeFFTVec(mTD.mCurrentFrame);
        for (int j=1; j<halfFFTSize; j++)
        {
            float diff = 20.0f * (log10(mTD.mCurrentFrame[j]) -  log10(mTD.mPrevFrame[j]));
            if (diff >= 0.0f)
                transientList[f] += diff;
        }
        transientList[f] /= halfFFTSize;
        if (transientList[f] > transientThresholdDB && f > 0){
            //To not have a back to back transient detected. We want them to be isolated
            transientList[f] = (transientList[f-1] == 1.0) ? 0.0f : 1.0f;              //1.0f; //Debug the brah brah
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

    // Optional 4th arg: pitch shift in semitones
    if (argc == 5)
        pitch_shift_semi = atoi(argv[4]);

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
    vector<PeakTrack> active_peaks;
    int peak_id_counter = 0;
    //We don't care about the signal once its below 70 db the highest tracked magnitude it obtained
    double threshold_factor = pow(10.0, (-40.0/20));
    
    vector<vector<PeakTrack>> frames_peaks;
    int printCoutner = 1;

    // Store previous frame's phase for instantaneous frequency calculation
    vector<double> prev_phase_spec;
    int prev_frame_size = 0;

    for (int frame_idx = 0; frame_idx < num_frames; frame_idx++) {
        vector<complex<double>> frameGuy = spec[frame_idx];
        frame_size = frameGuy.size()*2;
        bool is_long_window = (frame_size == LONG_SIZE);
        vector<double> mag_spec(frame_size/2, 0.0);
        vector<double> phase_spec(frame_size/2, 0.0);

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

            // Peak detection on the higher-resolution 4096 spectrum
            vector<int> peaks = detect_peaks(analysis_mag, threshold);
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

                        // Smooth blending: trust inst_freq more at low frequencies
                        // (where phase difference is more reliable) and parabolic
                        // more at high frequencies (where phase noise increases).
                        // Sigmoid-like curve: inst_weight goes from ~0.95 at 0 Hz
                        // to ~0.2 at high frequencies, centered around 1000 Hz.
                        double log_f = log10(max(parabolic_freq_hz, 20.0));
                        // Map: 20 Hz (log=1.3) → 0.95, 1000 Hz (log=3.0) → 0.5, 10000 Hz (log=4.0) → 0.15
                        double inst_weight = 0.95 - 0.8 * (log_f - 1.3) / (4.0 - 1.3);
                        inst_weight = max(0.15, min(0.95, inst_weight));
                        freqs_hz[i] = inst_weight * inst_freq + (1.0 - inst_weight) * parabolic_freq_hz;
                    } else {
                        freqs_hz[i] = parabolic_freq_hz;
                    }
                }

                for (size_t i = 0; i < peaks.size(); i++) {
                    double f_hz = freqs_hz[i];
                    double m_db = mags[i];
                    // Convert 4096-domain bin to 2048-domain for peak tracking
                    int p_bin = (int)round((double)peaks[i] / 2.0);
                    double ph = phases[i];

                    double peak_tol = (p_bin >= 5) ? 2.0 : 1.0;
                    int match_idx = find_best_match_peak(p_bin, active_peaks, peak_tol);
                    if (match_idx != -1) {
                        // Always mark as matched so the unmatched-update path
                        // (which uses the 2048 synthesis spectrum) never overwrites
                        // parameters that were measured from the 4096 analysis spectrum.
                        active_peaks[match_idx].edit = true;

                        // Conditional frequency smoothing for low-freq tracked peaks
                        if (f_hz < 200.0 && active_peaks[match_idx].freq_hz > 0.0) {
                            double freq_delta_pct = fabs(f_hz - active_peaks[match_idx].freq_hz) / active_peaks[match_idx].freq_hz;
                            if (freq_delta_pct < 0.05)
                                f_hz = 0.3 * f_hz + 0.7 * active_peaks[match_idx].freq_hz;
                        }

                        // Only update freq/phase/bin from the strongest match
                        // (protects against weaker duplicate matches at nearby bins)
                        if (m_db > active_peaks[match_idx].current_db) {
                            active_peaks[match_idx].freq_hz = f_hz;
                            active_peaks[match_idx].peak_bin = p_bin;
                            active_peaks[match_idx].phase = ph;
                        }

                        // Always update magnitude and death threshold
                        active_peaks[match_idx].current_db = m_db;
                        active_peaks[match_idx].analysis_fft_size = ANALYSIS_SIZE;
                        if (m_db > active_peaks[match_idx].max_db)
                            active_peaks[match_idx].max_db = m_db;
                        double thresholdDB = threshold_factor * active_peaks[match_idx].max_db;
                        if (m_db < thresholdDB)
                            active_peaks[match_idx].alive = false;
                    } else {
                        PeakTrack newPeak(peak_id_counter++, f_hz, m_db, p_bin, ph, ANALYSIS_SIZE);
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
                int long_frame_size = (int)long_spec.size() * 2; // LONG_SIZE
                vector<double> long_mag_spec(long_spec.size(), 0.0);
                vector<double> long_phase_spec(long_spec.size(), 0.0);
                for (int k = 0; k < (int)long_spec.size(); k++) {
                    long_mag_spec[k]   = abs(long_spec[k]);
                    long_phase_spec[k] = arg(long_spec[k]);
                }

                vector<int> peaks = detect_peaks(long_mag_spec, threshold);
                vector<double> freqs, mags;
                if (peaks.size() > 0) {
                    parabolic_interpolation(long_mag_spec, peaks, freqs, mags);
                    for (size_t i = 0; i < peaks.size(); i++) {
                        double f_hz  = freqs[i] * ((double)sr / (double)long_frame_size);
                        double m_db  = mags[i];
                        int    p_bin = peaks[i];
                        double ph    = long_phase_spec[p_bin];

                        double peak_tol = (p_bin >= 5) ? 2.0 : 1.0;
                        int match_idx = find_best_match_peak(p_bin, active_peaks, peak_tol);
                        if (match_idx != -1) {
                            active_peaks[match_idx].freq_hz    = f_hz;
                            active_peaks[match_idx].current_db = m_db;
                            active_peaks[match_idx].peak_bin   = p_bin;
                            active_peaks[match_idx].phase      = ph;
                            active_peaks[match_idx].edit       = true;
                            active_peaks[match_idx].analysis_fft_size = LONG_SIZE;
                            if (m_db > active_peaks[match_idx].max_db)
                                active_peaks[match_idx].max_db = m_db;
                            double thresholdDB = threshold_factor * active_peaks[match_idx].max_db;
                            if (m_db < thresholdDB)
                                active_peaks[match_idx].alive = false;
                        } else {
                            PeakTrack newPeak(peak_id_counter++, f_hz, m_db, p_bin, ph, LONG_SIZE);
                            active_peaks.push_back(newPeak);
                        }
                    }
                }
            }
        }

        for (auto &ap : active_peaks) {
            if (ap.alive && !ap.edit) {
                double scaled_bin;
                // For long frames, use the 4096 analysis spectrum so unmatched peaks
                // stay on the same magnitude/phase scale as matched peaks.
                // For short frames, use the 2048 synthesis spectrum as before.
                if (is_long_window && !analysis_mag.empty()) {
                    // Scale peak_bin from 2048 domain to 4096 domain
                    scaled_bin = ap.peak_bin * 2.0;
                    if (scaled_bin >= 0 && scaled_bin < (int)analysis_mag.size() - 1) {
                        double true_freq, true_mag;
                        single_parabolic_interpolation(analysis_mag, scaled_bin, true_freq, true_mag);
                        ap.current_db = true_mag;
                        ap.phase = interpolate_phase(analysis_phase_store, scaled_bin);
                        if (true_mag > ap.max_db) ap.max_db = true_mag;
                        double thresholdDB = threshold_factor * ap.max_db;
                        if (true_mag < thresholdDB) ap.alive = false;
                    } else {
                        ap.alive = false;
                    }
                } else {
                    if (is_long_window) {
                        scaled_bin = ap.peak_bin;
                    } else {
                        scaled_bin = ap.peak_bin * (double)SHORT_SIZE / LONG_SIZE;
                    }
                    if (scaled_bin >= 0 && scaled_bin < (int)mag_spec.size() - 1) {
                        double true_freq, true_mag;
                        single_parabolic_interpolation(mag_spec, scaled_bin, true_freq, true_mag);
                        ap.current_db = true_mag;
                        ap.phase = interpolate_phase(phase_spec, scaled_bin);
                        if (true_mag > ap.max_db) ap.max_db = true_mag;
                        double thresholdDB = threshold_factor * ap.max_db;
                        if (true_mag < thresholdDB) ap.alive = false;
                    } else {
                        ap.alive = false;
                    }
                }
            }
            ap.edit = false;
        }
        
        //Remove the unalive
        {
            vector<PeakTrack> temp;
            for (auto &p : active_peaks) {
                if (p.alive) temp.push_back(p);
            }
            active_peaks.swap(temp);
        }
        vector<PeakTrack> frame_info;
        for (auto &ap : active_peaks) {
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

    }
    
    
    //==========================================================================
    // SYNTHESIS — Center-Relative Oscillator Bank with 4096 Analysis
    //==========================================================================
    frame_size = LONG_SIZE;
    float total_length = (float)lengthYouNeed;
    vector<float> synthesized_signal(total_length, 0.0f);
    
    vector<float> frame_signal(LONG_SIZE, 0.0f);
    vector<float> frame_signal_short(SHORT_SIZE, 0.0f);

    double nyquist = (double)sr / 2.0;

    vector<int> chordIntervals = {0};

    for (int frame_idx = 0; frame_idx < num_frames; frame_idx++) {
        SynthInformation current_information = containsSynthPlacement[frame_idx];
        int frame_size = current_information.size;
        bool is_long = (frame_size == LONG_SIZE);
        int center = frame_size / 2;

        double t_center_abs = (double)(current_information.start + center) / (double)sr;

        fill(frame_signal.begin(), frame_signal.end(), 0.0f);
        fill(frame_signal_short.begin(), frame_signal_short.end(), 0.0f);

        bool shorter = false;

        for (auto &peak : frames_peaks[frame_idx]) {
            double freq = peak.freq_hz;
            double mag = 4.0 * peak.current_db / (double)peak.analysis_fft_size;
            double phase = peak.phase;

            // Advance phase from FFT start to synthesis frame center.
            double phase_advance_samples = (peak.analysis_fft_size == ANALYSIS_SIZE)
                ? (double)(ANALYSIS_SIZE / 2)
                : (double)center;
            double phase_at_center = phase + 2.0 * M_PI * freq * phase_advance_samples / (double)sr;

            for (int interval : chordIntervals) {
                double shiftFactor = pow(2.0, (interval + pitch_shift_semi) / 12.0);
                double shiftedFreq = freq * shiftFactor;

                if (shiftedFreq >= nyquist || shiftedFreq <= 0.0) continue;

                double shift_phase_offset = 2.0 * M_PI * (shiftedFreq - freq) * t_center_abs;
                shift_phase_offset = fmod(shift_phase_offset, 2.0 * M_PI);
                double correctedPhase = phase_at_center + shift_phase_offset;

                for (int n = 0; n < frame_size; n++) {
                    double t_centered = static_cast<double>(n - center) / sr;
                    double sample = mag * cos(2.0 * M_PI * shiftedFreq * t_centered + correctedPhase);
                    if (is_long) {
                        frame_signal[n] += static_cast<float>(sample);
                    } else {
                        frame_signal_short[n] += static_cast<float>(sample);
                    }
                }
                shorter = !is_long;
            }
        }
        
        // Synthesis window
        if (frame_idx == 0 && !current_information.trans) {
            if (shorter) {
                for (int i = 0; i < frame_size; i++)
                    frame_signal_short[i] *= rect_fade_to_hann_short[i];
            } else {
                for (int i = 0; i < frame_size; i++)
                    frame_signal[i] *= rect_fade_to_hann[i];
            }
        } else {
            vector<float> window_used = current_information.windowApplied;
            if (is_long) {
                for (int i = 0; i < frame_size; i++)
                    frame_signal[i] *= window_used[i];
            } else {
                for (int i = 0; i < frame_size; i++)
                    frame_signal_short[i] *= window_used[i];
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
                synthesized_signal[i] += frame_signal[i - start];
            }
        } else {
            for (int i = start; i < end; i++) {
                synthesized_signal[i] += frame_signal_short[i - start];
            }
        }
    }

    // Output limiting — only scale if clipping
    float max_val = 0.0f;
    for (auto val : synthesized_signal) {
        float abs_val = fabs(val);
        if (abs_val > max_val) max_val = abs_val;
    }
    if (max_val > 1.0f) {
        float scale = 0.95f / max_val;
        for (auto &val : synthesized_signal) {
            val *= scale;
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
    if (argc < 4 || argc > 5)
    {
        return false;
    }

    settings.inputWavFilePath = argv[1];
    settings.outputWavFilePath = argv[2];
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
    const int ChunkSize = 4 + (8 + Subchunk1Size) + (8 + Subchunk2Size);
    const short AudioFormat = 1;
    const short BitsPerSample = 16;
    const short BlockAlign = numChannels * (BitsPerSample / 8);
    const int ByteRate = sampleRate * BlockAlign;
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



