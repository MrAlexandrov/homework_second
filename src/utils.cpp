#include "utils.hpp"

#include <random>

namespace NUtils {

constexpr double EPS = 1e-6;

double GenerateRandomNumber() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    return dis(gen);
}

int GenerateIntNumber(int from, int to) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(from, to);

    return dis(gen);
}

bool Equals(double lhs, double rhs) {
    return std::fabs(lhs - rhs) < EPS;
}

} // namespace NUtils
