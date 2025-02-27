#include "analytical.hpp"

namespace NAnalitycal {

TAnalyticalSolution::TAnalyticalSolution(const TMatrix& P)
    : P_(P)
{
}

std::vector<Type> TAnalyticalSolution::CalculateAndGetDistribution() {
    CalculateDistribution();
    return GetDistribution();
}

std::vector<Type> TAnalyticalSolution::GetDistribution() const {
    return Distribution_;
}

void TAnalyticalSolution::CalculateDistribution() {
    TVector distribution = CalculateDistributionImpl();
    size_t n = distribution.rows();
    Distribution_.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        Distribution_.emplace_back(distribution(i));
    }
}

TVector TAnalyticalSolution::CalculateDistributionImpl() const {
    size_t n = P_.rows();

    TMatrix A = P_.transpose() - TMatrix::Identity(n, n);
    A.row(n - 1).setOnes();

    TVector b = TVector::Zero(n);
    b(n - 1) = Type(1);

    return A.fullPivLu().solve(b);
}

} // namespace NAnalitycal
