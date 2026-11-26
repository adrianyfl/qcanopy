// `qcanopy sim`: run a circuit on one input state (sim.hpp) and print the output amplitudes.

#include "commands.hpp"
#include "common.hpp"

#include "qcanopy/sim.hpp"
#include "qcanopy/stimulus.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace qcanopy::cli {

    namespace {
        const char* const help =
            "qcanopy sim FILE [options]\n"
            "Simulate a circuit (OpenQASM 3.0, or JSON IR v1 if *.json) and print its output state.\n"
            "\n"
            "  --input SPEC   Input state (default zero):\n"
            "                   zero           |0...0>\n"
            "                   basis:INDEX    basis state |INDEX> (qubit k = bit k)\n"
            "                   haar:SEED      Haar-random state\n"
            "                   product:SEED   random product state\n"
            "  --top K        Show the K most likely basis states (default 16).\n"
            "  --all          Show every basis state with nonzero probability.\n"
            "  --json         Print the shown amplitudes as one JSON object.\n"
            "  --help         Show this help.\n"
            "\n"
            "Basis states are written |q[n-1] ... q[0]>. Probabilities below 1e-12 are hidden.\n"
            "The state vector needs 2^n x 16 bytes (n = 30 -> 16 GiB).\n";

        constexpr double kShowThreshold = 1e-12;

        StateVector input_state(int n, const std::string& spec) {
            if (spec == "zero") return StateVector(n);
            const auto colon = spec.find(':');
            const std::string kind = spec.substr(0, colon);
            if (colon == std::string::npos) throw UsageError("--input expects zero, basis:INDEX, haar:SEED or product:SEED");
            const uint64_t value = parse_uint("--input " + kind, spec.substr(colon + 1));
            if (kind == "basis") {
                if (n < 64 && value >> n) throw UsageError("--input basis index " + std::to_string(value) + " needs more than " + std::to_string(n) + " qubits");
                return basis_state(n, value);
            }
            if (kind == "haar") return haar_random_state(n, value);
            if (kind == "product") return random_product_state(n, value);
            throw UsageError("--input expects zero, basis:INDEX, haar:SEED or product:SEED");
        }

        std::string bits(uint64_t index, int n) {
            std::string s(static_cast<std::size_t>(n), '0');
            for (int q = 0; q < n; ++q) if ((index >> q) & 1) s[static_cast<std::size_t>(n - 1 - q)] = '1';
            return s;
        }
    } // namespace

    int sim_main(int argc, const char* const* argv) {
        return guarded("sim", [&] {
            const Args args(argc, argv, {"--all", "--json"}, {"--input", "--top"});
            if (args.help()) { std::cout << help; return kExitOk; }
            args.expect_positional(1, 1, "one circuit file");
            if (args.has("--all") && args.has("--top")) throw UsageError("--all and --top cannot be combined");
            const std::size_t top = args.has("--all") ? SIZE_MAX
                : static_cast<std::size_t>(parse_uint("--top", args.get("--top").value_or("16")));
            const std::string spec = args.get("--input").value_or("zero");

            const std::string& file = args.positional()[0];
            const Circuit c = load_circuit(file);
            const StateVector out = Simulator{}.execute(c, input_state(c.n, spec));

            // Nonzero entries, most likely first (ties by index).
            std::vector<uint64_t> shown;
            for (uint64_t i = 0; i < out.amplitudes.size(); ++i)
                if (std::norm(out.amplitudes[i]) > kShowThreshold) shown.push_back(i);
            const std::size_t nonzero = shown.size();
            std::stable_sort(shown.begin(), shown.end(), [&](uint64_t x, uint64_t y) {
                return std::norm(out.amplitudes[x]) > std::norm(out.amplitudes[y]);
            });
            if (shown.size() > top) shown.resize(top);

            if (args.has("--json")) {
                std::cout << std::setprecision(17) << "{\"file\":" << json_string(file) << ",\"qubits\":" << c.n
                          << ",\"gates\":" << c.size() << ",\"input\":" << json_string(spec)
                          << ",\"nonzero\":" << nonzero << ",\"amplitudes\":[";
                for (std::size_t k = 0; k < shown.size(); ++k) {
                    const uint64_t i = shown[k];
                    const cplx a = out.amplitudes[i];
                    std::cout << (k ? "," : "") << "{\"index\":" << i << ",\"bits\":\"" << bits(i, c.n)
                              << "\",\"re\":" << a.real() << ",\"im\":" << a.imag() << ",\"p\":" << std::norm(a) << "}";
                }
                std::cout << "]}\n";
                return kExitOk;
            }

            std::cout << c.n << " qubits, " << c.size() << " gates, input " << spec << "\n" << std::fixed;
            for (const uint64_t i : shown) {
                const cplx a = out.amplitudes[i];
                std::cout << "  |" << bits(i, c.n) << ">  " << std::showpos << std::setprecision(6) << a.real() << ' '
                          << a.imag() << "i" << std::noshowpos << "   p=" << std::norm(a) << "\n";
            }
            std::cout << std::defaultfloat;
            if (shown.size() < nonzero)
                std::cout << "  ... " << nonzero - shown.size() << " more (use --top K or --all)\n";
            std::cout << "(" << nonzero << " of " << out.amplitudes.size() << " basis states have p > 1e-12)\n";
            return kExitOk;
        });
    }

} // namespace qcanopy::cli
