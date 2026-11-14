#include "qcanopy/stimulus.hpp"

#include <cmath>
#include <stdexcept>

namespace qcanopy {
    namespace {
        constexpr double kTwoPi = 6.28318530717958647692;

        // SplitMix64 finalizer.
        uint64_t mix(uint64_t z) {
            z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
            return z ^ (z >> 31);
        }

        struct Rng {
            uint64_t state;
            uint64_t next() { state += UINT64_C(0x9e3779b97f4a7c15); return mix(state); }
            // Uniform on (0, 1]: never 0, so log() below is always finite.
            double unit_open0() { return (static_cast<double>(next() >> 11) + 1.0) * 0x1.0p-53; }
        };

        // Box-Muller: one standard complex Gaussian z = x + iy with x, y ~ N(0, 1/2) i.i.d.
        // (The variance scale is irrelevant since the vector is normalized afterwards.)
        cplx complex_gaussian(Rng& r) {
            const double radius = std::sqrt(-std::log(r.unit_open0()));
            const double theta = kTwoPi * r.unit_open0();
            return std::polar(radius, theta);
        }
    } // namespace

    StateVector haar_random_state(int n, uint64_t seed) {
        if (n < 0 || n > kMaxQubits)
            throw std::invalid_argument("haar_random_state: qubit count out of range");

        StateVector state(n);
        Rng rng{mix(seed)};

        double norm_sq = 0.0;
        for (auto& amp : state.amplitudes) {
            amp = complex_gaussian(rng);
            norm_sq += std::norm(amp);
        }
        // An amplitude is exactly 0 only if its radius draw is exactly 1.0 (probability 2^-53),
        // so norm_sq == 0 is not a practical concern.
        const double inv_norm = 1.0 / std::sqrt(norm_sq);
        for (auto& amp : state.amplitudes) amp *= inv_norm;
        return state;
    }

    std::vector<StateVector> haar_random_states(int n, std::size_t count, uint64_t seed) {
        std::vector<StateVector> states;
        states.reserve(count);
        for (std::size_t k = 0; k < count; ++k) {
            const uint64_t sub_seed = mix(seed ^ mix(k + UINT64_C(0x9e3779b97f4a7c15)));
            states.push_back(haar_random_state(n, sub_seed));
        }
        return states;
    }

} // namespace qcanopy
