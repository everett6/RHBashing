#include "riemann.hpp"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#ifdef RH_CUDA
#include <cuda_runtime.h>
static void cuda_check(cudaError_t status) {
    if (status != cudaSuccess) throw std::runtime_error(cudaGetErrorString(status));
}

__global__ void evaluate(const rh::Term* terms, double start, double end,
                         double step, uint64_t offset, uint64_t last, int count, double* out) {
    const int lane = threadIdx.x & 31;
    const int sample = (blockIdx.x * blockDim.x + threadIdx.x) / 32;
    if (sample >= count) return; // Entire warp takes the same branch.
    const uint64_t index = offset + sample;
    // Match host sampling coordinates, without a fused multiply-add.
    const double t = index == last ? end : fmin(end,
        __dadd_rn(start, __dmul_rn(static_cast<double>(index), step)));
    const double a = sqrt(t / (2.0 * rh::pi));
    const int n = static_cast<int>(a);
    const double angle = rh::theta(t);
    double sum = 0.0;
    for (int k = lane + 1; k <= n; k += 32) {
        const rh::Term term = terms[k - 1];
        sum += term.inv_sqrt_n * cos(angle - t * term.log_n);
    }
    for (int shift = 16; shift; shift /= 2)
        sum += __shfl_down_sync(0xffffffff, sum, shift);
    if (lane == 0) out[sample] = 2.0 * sum + rh::correction(a, n);
}

struct Gpu {
    rh::Term* terms = nullptr;
    double* output = nullptr;
    double* host = nullptr;
    ~Gpu() { cudaFree(terms); cudaFree(output); cudaFreeHost(host); }
};
#endif

static double number(const char* text) {
    char* end = nullptr;
    const double value = std::strtod(text, &end);
    if (end == text || *end || !std::isfinite(value))
        throw std::runtime_error("Arguments must be finite numbers");
    return value;
}

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("Usage: rhscan START END STEP");
        const double start = number(argv[1]), end = number(argv[2]), step = number(argv[3]);
        if (start < 200 || end > 1e9 || end <= start || step <= 0 || start + step == start)
            throw std::runtime_error("Require 200 <= START < END <= 1e9 and a positive representable STEP");
        const double intervals = std::ceil((end - start) / step);
        if (!std::isfinite(intervals) || intervals > 1e12)
            throw std::runtime_error("Too many samples (limit: 1e12 intervals); split the range");
        const uint64_t total = static_cast<uint64_t>(intervals) + 1;
        constexpr int batch = 65536;
        const int max_n = static_cast<int>(std::sqrt(end / (2.0 * rh::pi)));
        // Bound work per launch at high t, particularly on desktop display GPUs.
        const int batch_count = std::max(256, std::min(batch, (8 * 1024 * 1024) / max_n));
        std::vector<rh::Term> terms(max_n);
        for (int k = 1; k <= max_n; ++k)
            terms[k - 1] = {std::log(static_cast<double>(k)), 1.0 / std::sqrt(static_cast<double>(k))};

#ifdef RH_CUDA
        Gpu gpu;
        cudaDeviceProp props{};
        cuda_check(cudaGetDeviceProperties(&props, 0));
        std::cerr << "GPU: " << props.name << '\n';
        cuda_check(cudaMalloc(reinterpret_cast<void**>(&gpu.terms), terms.size() * sizeof(rh::Term)));
        cuda_check(cudaMalloc(reinterpret_cast<void**>(&gpu.output), batch * sizeof(double)));
        cuda_check(cudaMallocHost(reinterpret_cast<void**>(&gpu.host), batch * sizeof(double)));
        cuda_check(cudaMemcpy(gpu.terms, terms.data(), terms.size() * sizeof(rh::Term), cudaMemcpyHostToDevice));
        double* values = gpu.host;
#else
        std::vector<double> output(batch);
        double* values = output.data();
#endif
        std::cerr << "Approximate critical-line scan: " << total << " samples, at most "
                  << max_n << " terms/sample. No completeness or RH claim.\n";
        std::cout << "t_left,t_right,z_approx_left,z_approx_right\n" << std::setprecision(17);
        const auto clock_start = std::chrono::steady_clock::now();
        double previous_t = 0, previous_z = 0;
        uint64_t candidates = 0;
        for (uint64_t offset = 0; offset < total; offset += batch_count) {
            const int count = static_cast<int>(std::min<uint64_t>(batch_count, total - offset));
#ifdef RH_CUDA
            evaluate<<<(count + 7) / 8, 256>>>(gpu.terms, start, end, step, offset, total - 1, count, gpu.output);
            cuda_check(cudaGetLastError());
            cuda_check(cudaMemcpy(values, gpu.output, count * sizeof(double), cudaMemcpyDeviceToHost));
#else
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < count; ++j) {
                const double t = offset + j == total - 1 ? end
                    : std::min(end, start + static_cast<double>(offset + j) * step);
                values[j] = rh::z_cpu(t, terms.data());
            }
#endif
            for (int j = 0; j < count; ++j) {
                const double t = offset + j == total - 1 ? end
                    : std::min(end, start + static_cast<double>(offset + j) * step);
                const double z = values[j];
                if (!std::isfinite(z)) throw std::runtime_error("Nonfinite Z value");
                if (offset + j && (z == 0 || previous_z == 0 || std::signbit(z) != std::signbit(previous_z))) {
                    std::cout << previous_t << ',' << t << ',' << previous_z << ',' << z << '\n';
                    ++candidates;
                }
                previous_t = t;
                previous_z = z;
            }
            if (!std::cout) throw std::runtime_error("Failed to write candidate output");
        }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock_start).count();
        std::cout.flush();
        if (!std::cout) throw std::runtime_error("Failed to flush candidate output");
        std::cerr << "Completed " << total << " samples, " << candidates
                  << " candidate brackets in " << seconds << " seconds.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
