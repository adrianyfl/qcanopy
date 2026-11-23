#include "cuda_sim.hpp"
#include <cuda_runtime.h>

namespace qcanopy::cuda {
    bool device_available() {
        int count = 0;
        // cudaGetDeviceCount fails (no driver) or reports 0 devices when there's no GPU.
        return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
    }
    void evolve(const Circuit& circuit, StateVector& state){

    }
}