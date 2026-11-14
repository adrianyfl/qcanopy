#include "qcanopy/stimulus.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace qcanopy {
    namespace {
        constexpr double kTwoPi = 6.28318530717958647692;

        // Seed domains for batch derivation, so different stimulus kinds never share sub-seeds.
        constexpr uint64_t kDomainHaar = 1;
        constexpr uint64_t kDomainProduct = 2;
        constexpr uint64_t kDomainBasis = 3;

        // SplitMix64 finalizer.
        uint64_t mix(uint64_t z) {
            z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
            return z ^ (z >> 31);
        }

        // Independent sub-seed for item `index` of stimulus kind `domain`.
        uint64_t derive(uint64_t seed, uint64_t domain, uint64_t index) {
            return mix(seed ^ mix(domain) ^ mix(index + UINT64_C(0x9e3779b97f4a7c15)));
        }

        struct Rng {
            uint64_t state;
            uint64_t next() { state += UINT64_C(0x9e3779b97f4a7c15); return mix(state); }
            // Uniform on (0, 1]: never 0, so log() below is always finite.
            double unit_open0() { return (static_cast<double>(next() >> 11) + 1.0) * 0x1.0p-53; }
        };

        // Box-Muller: one standard complex Gaussian z = x + iy with x, y ~ N(0, 1/2) i.i.d.
        // (The variance scale is irrelevant since vectors are normalized afterwards.)
        cplx complex_gaussian(Rng& r) {
            const double radius = std::sqrt(-std::log(r.unit_open0()));
            const double theta = kTwoPi * r.unit_open0();
            return std::polar(radius, theta);
        }

        void check_qubits(int n, const char* what) {
            if (n < 0 || n > kMaxQubits)
                throw std::invalid_argument(std::string(what) + ": qubit count out of range");
        }

        // Fills `count` states from `make(n, sub_seed)`, with sub-seeds drawn from `domain`.
        template <class Make>
        std::vector<StateVector> batch(int n, std::size_t count, uint64_t seed, uint64_t domain, Make make) {
            std::vector<StateVector> states;
            states.reserve(count);
            for (std::size_t k = 0; k < count; ++k) states.push_back(make(n, derive(seed, domain, k)));
            return states;
        }
    } // namespace

    StateVector haar_random_state(int n, uint64_t seed) {
        check_qubits(n, "haar_random_state");

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
        return batch(n, count, seed, kDomainHaar, haar_random_state);
    }

    StateVector random_product_state(int n, uint64_t seed) {
        check_qubits(n, "random_product_state");

        StateVector state(n);
        auto& amps = state.amplitudes;
        Rng rng{mix(seed)};

        // Build the tensor product one qubit at a time. After processing qubits 0..k-1, the
        // first 2^k entries hold their product state. Qubit k is bit k of the index, so
        // tensoring in |v_k> = v0|0> + v1|1> maps entry j to j (bit k = 0) and j + 2^k (bit k = 1).
        amps[0] = {1.0, 0.0};
        for (int k = 0; k < n; ++k) {
            cplx v0 = complex_gaussian(rng), v1 = complex_gaussian(rng);
            const double inv_norm = 1.0 / std::sqrt(std::norm(v0) + std::norm(v1));
            v0 *= inv_norm;
            v1 *= inv_norm;

            const std::size_t half = std::size_t{1} << k;
            for (std::size_t j = 0; j < half; ++j) {
                amps[j + half] = amps[j] * v1;
                amps[j] *= v0;
            }
        }
        return state;
    }

    std::vector<StateVector> random_product_states(int n, std::size_t count, uint64_t seed) {
        return batch(n, count, seed, kDomainProduct, random_product_state);
    }

    StateVector basis_state(int n, uint64_t index) {
        check_qubits(n, "basis_state");
        // n <= kMaxQubits = 32, so 2^n fits in uint64_t.
        if (index >= (uint64_t{1} << n))
            throw std::invalid_argument("basis_state: index out of range for qubit count");

        StateVector state(n);   // |0...0>
        state.amplitudes[0] = {0.0, 0.0};
        state.amplitudes[static_cast<std::size_t>(index)] = {1.0, 0.0};
        return state;
    }

    StateVector random_basis_state(int n, uint64_t seed) {
        check_qubits(n, "random_basis_state");
        Rng rng{mix(seed)};
        // The top n bits of a uniform 64-bit word are uniform on [0, 2^n): no modulo bias.
        const uint64_t index = n == 0 ? 0 : rng.next() >> (64 - n);
        return basis_state(n, index);
    }

    std::vector<StateVector> random_basis_states(int n, std::size_t count, uint64_t seed) {
        return batch(n, count, seed, kDomainBasis, random_basis_state);
    }

    std::vector<StateVector> all_basis_states(int n) {
        check_qubits(n, "all_basis_states");
        const uint64_t dim = uint64_t{1} << n;
        std::vector<StateVector> states;
        states.reserve(static_cast<std::size_t>(dim));
        for (uint64_t i = 0; i < dim; ++i) states.push_back(basis_state(n, i));
        return states;
    }

} // namespace qcanopy
