// GPU vs CPU: runs random circuits on both paths and checks the amplitudes match.
// Exits with 77 (reported by ctest as "Skipped") when the GPU path isn't available.
#include "qcanopy/sim.hpp"
#include "qcanopy/stimulus.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>

using namespace qcanopy;

namespace {
    int failures = 0;
    constexpr double kTol = 1e-10;

    double max_diff(const StateVector& a, const StateVector& b) {
        if (a.amplitudes.size() != b.amplitudes.size()) return INFINITY;
        double d = 0;
        for (std::size_t i = 0; i < a.amplitudes.size(); ++i)
            d = std::max(d, std::abs(a.amplitudes[i] - b.amplitudes[i]));
        return d;
    }

    void expect_close(const StateVector& cpu, const StateVector& gpu, const char* what, int n, uint64_t seed) {
        const double d = max_diff(cpu, gpu);
        if (!(d < kTol)) {
            std::printf("FAIL %s: n=%d seed=%llu max |cpu-gpu| = %g\n",
                        what, n, static_cast<unsigned long long>(seed), d);
            ++failures;
        }
    }

    // Random mix of fixed, parameterized and 2-qubit gates; 2-qubit gates use
    // both orders (q0 < q1 and q0 > q1), since the kernel handles them differently.
    Circuit random_circuit(int n, int num_gates, uint64_t seed) {
        using K = GateKind;
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> angle(-M_PI, M_PI);
        const K one_q[] = {K::H, K::X, K::Y, K::Z, K::S, K::Sdg, K::T, K::SX, K::RX, K::RY, K::RZ, K::P, K::U};
        const K two_q[] = {K::CX, K::CZ, K::SWAP};
        Circuit c(n);
        for (int g = 0; g < num_gates; ++g) {
            const auto q0 = static_cast<uint8_t>(rng() % n);
            if (n > 1 && rng() % 3 == 0) {
                auto q1 = static_cast<uint8_t>(rng() % (n - 1));
                if (q1 >= q0) ++q1;   // any qubit except q0, either side of it
                c.add(make_gate(two_q[rng() % 3], q0, q1));
            } else {
                c.add(make_gate(one_q[rng() % 13], q0, {angle(rng), angle(rng), angle(rng)}));
            }
        }
        c.global_phase = angle(rng);
        return c;
    }
} // namespace

int main() {
    const Simulator cpu(false), gpu(true);
    if (!gpu.uses_gpu()) {
        std::printf("SKIP: GPU path not available (CPU-only build, or no GPU on this node)\n");
        return 77;
    }

    using clock = std::chrono::steady_clock;
    for (int n : {1, 2, 3, 5, 8, 12, 16, 20}) {
        double t_cpu = 0, t_gpu = 0;
        for (uint64_t seed = 1; seed <= 3; ++seed) {
            const Circuit c = random_circuit(n, 300, seed * 1000 + n);

            // Single run from |0...0>.
            auto t0 = clock::now();
            const StateVector want = cpu.execute(c);
            auto t1 = clock::now();
            const StateVector got = gpu.execute(c);
            auto t2 = clock::now();
            t_cpu += std::chrono::duration<double, std::milli>(t1 - t0).count();
            t_gpu += std::chrono::duration<double, std::milli>(t2 - t1).count();
            expect_close(want, got, "execute", n, seed);

            // Batch of random inputs, plus the streaming API.
            if (n <= 12) {
                const Stimulus stim = haar_random_states(n, 4, seed);
                const auto want_b = cpu.execute(c, stim);
                const auto got_b  = gpu.execute(c, stim);
                for (std::size_t i = 0; i < want_b.size(); ++i)
                    expect_close(want_b[i], got_b[i], "execute(stimulus)", n, seed);

                gpu.execute_each(c, stim, [&](std::size_t i, const StateVector& s) {
                    expect_close(want_b[i], s, "execute_each", n, seed);
                    return Control::Continue;
                });
            }
        }
        std::printf("n=%2d  cpu %8.2f ms   gpu %8.2f ms   (3 circuits x 300 gates)\n", n, t_cpu, t_gpu);
    }

    if (failures) { std::printf("%d failure(s)\n", failures); return 1; }
    std::printf("PASS: GPU matches CPU to %g\n", kTol);
    return 0;
}
