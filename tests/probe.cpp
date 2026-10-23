#include "../src/riemann.hpp"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    std::cout << std::setprecision(17);
    for (int j = 1; j < argc; ++j) {
        const double t = std::strtod(argv[j], nullptr);
        const int n = static_cast<int>(std::sqrt(t / (2 * rh::pi)));
        std::vector<rh::Term> terms(n);
        for (int k = 1; k <= n; ++k)
            terms[k - 1] = {std::log(static_cast<double>(k)), 1 / std::sqrt(static_cast<double>(k))};
        std::cout << rh::z_cpu(t, terms.data()) << '\n';
    }
}
