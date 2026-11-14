#include "qcanopy/equiv.hpp"

#include <algorithm>
#include <cmath>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace qcanopy {

    Circuit miter(const Circuit& a, const Circuit& b) {
        if (a.n != b.n) throw std::invalid_argument("miter: circuits have different qubit counts");
        Circuit m(a.n);
        m.gates.reserve(a.size() + b.size());
        m.gates = a.gates;
        for (auto it = b.gates.rbegin(); it != b.gates.rend(); ++it) m.add(inverse(*it));
        m.global_phase = a.global_phase - b.global_phase;
        return m;
    }

    EquivalenceChecker::EquivalenceChecker(double tol) : tol_(tol) {
        if (!(tol > 0.0)) throw std::invalid_argument("EquivalenceChecker: tolerance must be positive");
    }

    EquivResult EquivalenceChecker::check(const Circuit& a, const Circuit& b, const Stimulus& stimulus) const {
        if (a.n != b.n)
            throw std::invalid_argument("EquivalenceChecker::check: circuits have different qubit counts");
        if (stimulus.num_qbits != a.n)
            throw std::invalid_argument("EquivalenceChecker::check: stimulus qubit count does not match circuits");
        if (stimulus.empty())
            throw std::invalid_argument("EquivalenceChecker::check: empty stimulus");

        const Circuit m = miter(a, b);
        const std::size_t count = stimulus.size();

        // Per-stimulus measurements, indexed by stimulus. Written and read under `mu`.
        std::vector<cplx> phase(count);          // u_i = o_i / |o_i|
        std::vector<double> residual(count);     // || M|psi_i> - u_i |psi_i> ||
        std::vector<char> done(count, 0);        // char, not bool: distinct bytes per slot
        std::mutex mu;

        // Online early-stop test. It must only return Stop on a genuine failure under the SAME
        // rule as the final verdict below (reference phase = stimulus 0), otherwise cancelled
        // stimuli could hide an unchecked pass. Stimulus 0 is never cancelled (only indices
        // above a stopping index are), so its phase is always available eventually; results
        // that finish before it are re-checked when it arrives.
        auto on_output = [&](std::size_t i, const StateVector& out) -> Control {
            const auto& psi = stimulus.states[i].amplitudes;
            const auto& amp = out.amplitudes;

            cplx overlap{0.0, 0.0};
            for (std::size_t k = 0; k < amp.size(); ++k) overlap += std::conj(psi[k]) * amp[k];
            const double mag = std::abs(overlap);
            const cplx u = mag > 0.0 ? overlap / mag : cplx{1.0, 0.0};

            double res_sq = 0.0;
            for (std::size_t k = 0; k < amp.size(); ++k) res_sq += std::norm(amp[k] - u * psi[k]);
            const double res = std::sqrt(res_sq);

            std::lock_guard<std::mutex> lock(mu);
            phase[i] = u;
            residual[i] = res;
            done[i] = 1;

            if (res > tol_) return Control::Stop;
            if (i == 0) {
                for (std::size_t j = 1; j < count; ++j)
                    if (done[j] && residual[j] <= tol_ && std::abs(phase[j] - u) > tol_) return Control::Stop;
                return Control::Continue;
            }
            if (done[0] && std::abs(u - phase[0]) > tol_) return Control::Stop;
            return Control::Continue;
        };

        sim_.execute_each(m, stimulus, on_output);

        // Final verdict: a deterministic scan in index order over the completed stimuli,
        // so the witness is the lowest-index failure among them.
        EquivResult result;
        result.stimuli_checked = static_cast<std::size_t>(std::count(done.begin(), done.end(), 1));
        const cplx ref = phase[0];
        double max_dev = 0.0;
        for (std::size_t i = 0; i < count; ++i) {
            if (!done[i]) continue;
            const double phase_dev = std::abs(phase[i] - ref);
            Mismatch why = Mismatch::None;
            double dev = 0.0;
            if (residual[i] > tol_) { why = Mismatch::State; dev = residual[i]; }
            else if (phase_dev > tol_) { why = Mismatch::Phase; dev = phase_dev; }

            if (why != Mismatch::None) {
                result.equivalent = false;
                result.mismatch = why;
                result.deviation = dev;
                result.failing_index = i;
                result.failing_kind = stimulus.kinds[i];
                return result;
            }
            max_dev = std::max({max_dev, residual[i], phase_dev});
        }

        // No failure means no Stop was ever returned, so every stimulus ran.
        if (result.stimuli_checked != count)
            throw std::logic_error("EquivalenceChecker::check: stimuli cancelled without a failure");
        result.equivalent = true;
        result.global_phase = ref;
        result.deviation = max_dev;
        return result;
    }

} // namespace qcanopy
