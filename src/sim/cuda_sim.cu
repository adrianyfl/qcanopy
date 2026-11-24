#include "cuda_sim.hpp"
#include <cuda_runtime.h>
#include <algorithm>  
#include <stdexcept>   
#include <string>      

namespace qcanopy::cuda {
    bool device_available() {
        int count = 0;
        // cudaGetDeviceCount fails (no driver) or reports 0 devices when there's no GPU.
        return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
    }
    namespace{
        __device__ inline double2 operator*(double2 a, double2 b) {
            return make_double2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
        }
        __device__ inline double2 operator+(double2 a, double2 b) {
            return make_double2(a.x + b.x, a.y + b.y);
        }
        struct Mat2 { double2 m[4]; };
        __global__ void apply_1q_gate(double2* amp, size_t half_dim, int q, Mat2 u) {
            size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
            if (i >= half_dim) return;
            size_t bit  = size_t{1} << q;
            size_t idx0 = ((i >> q) << (q + 1)) | (i & (bit - 1));
            size_t idx1 = idx0 | bit;
            // load a0, a1;
            double2 a0 = amp[idx0];
            double2 a1 = amp[idx1];
            //write u00*a0 + u01*a1 and u10*a0 + u11*a1
            amp[idx0] = u.m[0]*a0 + u.m[1]*a1;
            amp[idx1] = u.m[2]*a0 + u.m[3]*a1;
        }
        struct Mat4 { double2 m[16]; };
        __global__ void apply_2q_gate(double2* amp, size_t quarter_dim, int q_min, int q_max, size_t mask0, size_t mask1, Mat4 u){
            size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
            if (i >= quarter_dim) return;
                
            // Insert 0 at the q_min position
            size_t temp = ((i >> q_min) << (q_min + 1)) | (i & ((1ULL << q_min) - 1));
            
            // Insert 0 at the q_max position
            size_t idx00 = ((temp >> q_max) << (q_max + 1)) | (temp & ((1ULL << q_max) - 1));
            // Produce the other 3 pairs
            size_t idx01 = idx00 | mask0;           // q0 = 1, q1 = 0
            size_t idx10 = idx00 | mask1;           // q0 = 0, q1 = 1
            size_t idx11 = idx00 | mask0 | mask1;   // q0 = 1, q1 = 1

            // Grab amplitudes
            double2 a00 = amp[idx00];
            double2 a01 = amp[idx01];
            double2 a10 = amp[idx10];
            double2 a11 = amp[idx11];

            // Apply 4x4 matrix multiplication
            amp[idx00] = u.m[0]*a00  + u.m[1]*a01  + u.m[2]*a10  + u.m[3]*a11;
            amp[idx01] = u.m[4]*a00  + u.m[5]*a01  + u.m[6]*a10  + u.m[7]*a11;
            amp[idx10] = u.m[8]*a00  + u.m[9]*a01  + u.m[10]*a10 + u.m[11]*a11;
            amp[idx11] = u.m[12]*a00 + u.m[13]*a01 + u.m[14]*a10 + u.m[15]*a11;
        }

        //------------host--------------
        // Throws with a readable message if a CUDA call failed.
        void check(cudaError_t err, const char* what) {
            if (err != cudaSuccess) {
                throw std::runtime_error(std::string("CUDA ") + what + ": " + cudaGetErrorString(err));
            }
        }

        // Owns a device allocation; frees it on scope exit, including when check() throws.
        struct DeviceBuffer {
            double2* ptr = nullptr;
            explicit DeviceBuffer(size_t n) { check(cudaMalloc(&ptr, n * sizeof(double2)), "cudaMalloc"); }
            ~DeviceBuffer() { cudaFree(ptr); }
            DeviceBuffer(const DeviceBuffer&) = delete;
            DeviceBuffer& operator=(const DeviceBuffer&) = delete;
        };

        constexpr unsigned kThreadsPerBlock = 256;

        // Number of blocks so that blocks * kThreadsPerBlock >= work_items.
        unsigned blocks_for(size_t work_items) {
            return static_cast<unsigned>((work_items + kThreadsPerBlock - 1) / kThreadsPerBlock);
        }

        double2 to_double2(const cplx& c) { return make_double2(c.real(), c.imag()); }
    }
    void evolve(const Circuit& circuit, StateVector& state){
        // std::complex<double> and double2 are both {real, imag} as two doubles,
        // so the amplitudes can be copied byte-for-byte.
        static_assert(sizeof(cplx) == sizeof(double2), "cplx and double2 must have the same layout");

        const size_t dim   = state.amplitudes.size();
        const size_t bytes = dim * sizeof(double2);

        // 1. Make sure the state fits (leave ~10% headroom for the CUDA context etc.).
        size_t free_b = 0, total_b = 0;
        check(cudaMemGetInfo(&free_b, &total_b), "cudaMemGetInfo");
        if (bytes > free_b / 10 * 9) {
            throw std::runtime_error("CUDA evolve: state needs " + std::to_string(bytes >> 20) +
                                        " MiB but only " + std::to_string(free_b >> 20) + " MiB is free on the GPU");
        }

        // 2. Allocate and copy the state to the GPU.
        DeviceBuffer d_amp(dim);
        check(cudaMemcpy(d_amp.ptr, state.amplitudes.data(), bytes, cudaMemcpyHostToDevice), "copy to device");

        // 3. One kernel launch per gate. Launches are queued and run in order on the
        //    default stream, so each gate sees the previous gate's result.
        for (const Gate& gate : circuit.gates) {
            if (gate.nq == 1) {
                Mat2 u;
                for (int k = 0; k < 4; ++k) u.m[k] = to_double2(gate.m[k]);
                const size_t half_dim = dim / 2;
                apply_1q_gate<<<blocks_for(half_dim), kThreadsPerBlock>>>(d_amp.ptr, half_dim, gate.q[0], u);
            } else if (gate.nq == 2) {
                Mat4 u;
                for (int k = 0; k < 16; ++k) u.m[k] = to_double2(gate.m[k]);
                const int q0 = gate.q[0], q1 = gate.q[1];
                const size_t quarter_dim = dim / 4;
                apply_2q_gate<<<blocks_for(quarter_dim), kThreadsPerBlock>>>(
                    d_amp.ptr, quarter_dim,
                    std::min(q0, q1), std::max(q0, q1),     // order for inserting zero bits
                    size_t{1} << q0, size_t{1} << q1,       // masks follow q0/q1 (matrix order)
                    u);
            } else {
                throw std::runtime_error("Simulator only supports 1-qubit and 2-qubit gates.");
            }
        }

        // 4. Launches don't report errors themselves; this catches bad launch configurations.
        check(cudaGetLastError(), "kernel launch");

        // 5. Copy back. Plain cudaMemcpy waits for all queued kernels to finish first,
        //    and also reports any error that happened while they ran.
        check(cudaMemcpy(state.amplitudes.data(), d_amp.ptr, bytes, cudaMemcpyDeviceToHost), "copy to host");

        // 6. Global phase on the CPU, exactly like the CPU path, so results match bit-for-bit here.
        state.apply_global_phase(circuit.global_phase);
    }

}