#pragma once
#include <cmath>

#ifdef __CUDACC__
#define RH_HD __host__ __device__
#else
#define RH_HD
#endif

namespace rh {
constexpr double pi = 3.141592653589793238462643383279502884;
struct Term { double log_n, inv_sqrt_n; };

RH_HD inline double theta(double t) {
    const double inv = 1.0 / t;
    const double inv2 = inv * inv;
    return 0.5 * t * (log(t / (2.0 * pi)) - 1.0) - pi / 8.0
         + inv * (1.0 / 48.0 + inv2 * (7.0 / 5760.0
         + inv2 * (31.0 / 80640.0)));
}

RH_HD inline double sinc(double x) {
    if (fabs(x) < 1e-4) {
        const double x2 = x * x;
        return 1.0 + x2 * (-1.0 / 6.0 + x2 / 120.0);
    }
    return sin(x) / x;
}

// Rewrite cos(2*pi*(p*p-p-1/16))/cos(2*pi*p) to remove
// the removable singularities at p=1/4 and p=3/4.
RH_HD inline double correction_coefficient(double p) {
    if (p < 0.5) {
        const double d = p - 0.25;
        return (0.75 - p) * sinc(2.0 * pi * d * (p - 0.75))
             / sinc(2.0 * pi * d);
    }
    const double d = p - 0.75;
    return (p - 0.25) * sinc(2.0 * pi * (p - 0.25) * d)
         / sinc(2.0 * pi * d);
}

RH_HD inline double correction(double a, int n) {
    return ((n & 1) ? 1.0 : -1.0)
         * correction_coefficient(a - n) / sqrt(a);
}

inline double z_cpu(double t, const Term* terms) {
    const double a = sqrt(t / (2.0 * pi));
    const int n = static_cast<int>(a);
    const double angle = theta(t);
    double sum = 0.0;
    for (int k = 1; k <= n; ++k)
        sum += terms[k - 1].inv_sqrt_n * cos(angle - t * terms[k - 1].log_n);
    return 2.0 * sum + correction(a, n);
}
} // namespace rh
