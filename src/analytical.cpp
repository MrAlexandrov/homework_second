#include "analytical.hpp"

namespace NAnalitycal {

TAnalyticalSolution::TAnalyticalSolution(const TMatrix& P)
    : P_(P)
    , Distribution_(CalculateStationaryDistribution())
{
}

std::vector<Type> TAnalyticalSolution::GetDistribution() const {
    std::vector<Type> result;
    size_t m = Distribution_.rows();
    result.reserve(m);
    for (int i = 0; i < m; ++i) {
        result.emplace_back(Distribution_(i));
    }
    return result;
}

TVector TAnalyticalSolution::CalculateStationaryDistribution() const {
    size_t n = P_.rows();

    TMatrix A = P_.transpose() - TMatrix::Identity(n, n);
    A.row(n - 1).setOnes();

    TVector b = TVector::Zero(n);
    b(n - 1) = Type(1);

    return A.fullPivLu().solve(b);
}

} // namespace NAnalitycal
