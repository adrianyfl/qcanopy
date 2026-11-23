#pragma once
#include "qcanopy/sim.hpp"
namespace qcanopy::cuda {
    bool device_available();
    void evolve(const Circuit& circuit, StateVector& state);
}