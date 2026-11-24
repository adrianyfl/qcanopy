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

    __device__ inline double2 operator*(double2 a, double2 b) {
        return make_double2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
    }
    __device__ inline double2 operator+(double2 a, double2 b) {
        return make_double2(a.x + b.x, a.y + b.y);
    }

    __global__ void apply_1q(double2* amp, size_t half_dim, int q, double2* m) {
        size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
        if (i >= half_dim) return;
        size_t bit  = size_t{1} << q;
        size_t idx0 = ((i >> q) << (q + 1)) | (i & (bit - 1));
        size_t idx1 = idx0 | bit;
        // load a0, a1;
        double2 a0 = amp[idx0];
        double2 a1 = amp[idx1];
        //write u00*a0 + u01*a1 and u10*a0 + u11*a1
        amp[idx0] = m[0]*a0 + m[1]*a1;
        amp[idx1] = m[2]*a0 + m[3]*a1;
    }
    __global__ void apply_2q_gate(double2* amp, size_t half_dim, int q0, int q1
                         double2 u00, double2 u01, double2 u10, double2 u11)){
        int q0 = gate.q[0]; // control qubit
        int q1 = gate.q[1]; // target qubit

        // We must insert the zeros from right to left (lowest bit first) 
        // so we don't mess up the offset for the second insertion.
        int q_min = std::min(q0, q1);
        int q_max = std::max(q0, q1);

        size_t mask0 = 1ULL << q0;
        size_t mask1 = 1ULL << q1;

        size_t quarter_dim = state.amplitudes.size() / 4;
        const auto& m = gate.m; // grab gate matrix

        for (size_t i = 0; i < quarter_dim; ++i) {
            
            // Insert 0 at the q_min position
            size_t temp = ((i >> q_min) << (q_min + 1)) | (i & ((1ULL << q_min) - 1));
            
            // Insert 0 at the q_max position
            size_t idx00 = ((temp >> q_max) << (q_max + 1)) | (temp & ((1ULL << q_max) - 1));

            // Produce the other 3 pairs
            size_t idx01 = idx00 | mask0;           // q0 = 1, q1 = 0
            size_t idx10 = idx00 | mask1;           // q0 = 0, q1 = 1
            size_t idx11 = idx00 | mask0 | mask1;   // q0 = 1, q1 = 1

            // Grab amplitudes
            cplx a00 = state.amplitudes[idx00];
            cplx a01 = state.amplitudes[idx01];
            cplx a10 = state.amplitudes[idx10];
            cplx a11 = state.amplitudes[idx11];

            // Apply 4x4 matrix multiplication
            state.amplitudes[idx00] = m[0]*a00  + m[1]*a01  + m[2]*a10  + m[3]*a11;
            state.amplitudes[idx01] = m[4]*a00  + m[5]*a01  + m[6]*a10  + m[7]*a11;
            state.amplitudes[idx10] = m[8]*a00  + m[9]*a01  + m[10]*a10 + m[11]*a11;
            state.amplitudes[idx11] = m[12]*a00 + m[13]*a01 + m[14]*a10 + m[15]*a11;
        }
    }
}