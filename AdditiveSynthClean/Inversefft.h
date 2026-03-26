//
//  Inversefft.h
//  AdditiveSynthChange
//
//  Created by Riley on 3/26/26.
//
#pragma once
#include <vector>
#include <complex>
#include <cmath>

// Full inverse: all bins → N time-domain samples.  O(N²)
inline void InverseRealFFT(float* output, int N,
                           const std::vector<std::complex<double>>& spectrum)
{
    int numBins = (int)spectrum.size();          // expected N/2
    for (int n = 0; n < N; n++) {
        double sum = spectrum[0].real();         // DC (k = 0)
        for (int k = 1; k < numBins; k++) {
            double angle = 2.0 * M_PI * k * n / static_cast<double>(N);
            // Real-signal symmetry: X[N-k] = conj(X[k]), so both halves contribute
            sum += 2.0 * (spectrum[k].real() * cos(angle)
                        - spectrum[k].imag() * sin(angle));
        }
        output[n] = static_cast<float>(sum / N);
    }
}

// Sparse inverse: only evaluates nonzero bins.  O(N * K), K ≪ N/2 typical.
// `masked` has the same size as the full spectrum but most entries are (0,0).
// We skip zero bins entirely, which makes this much faster than the full version
// when only a fraction of bins are active (e.g. peaks ±2 = ~30% of bins).
inline void InverseRealFFT_Sparse(float* output, int N,
                                  const std::vector<std::complex<double>>& masked)
{
    int numBins = (int)masked.size();
    // Collect nonzero bin indices once
    std::vector<int> active;
    active.reserve(numBins);
    for (int k = 0; k < numBins; k++) {
        if (masked[k].real() != 0.0 || masked[k].imag() != 0.0)
            active.push_back(k);
    }

    for (int n = 0; n < N; n++) {
        double sum = 0.0;
        for (int k : active) {
            if (k == 0) {
                sum += masked[0].real();
            } else {
                double angle = 2.0 * M_PI * k * n / static_cast<double>(N);
                sum += 2.0 * (masked[k].real() * cos(angle)
                            - masked[k].imag() * sin(angle));
            }
        }
        output[n] = static_cast<float>(sum / N);
    }
}
