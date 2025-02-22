#pragma once

#include "types.hpp"

#include <vector>

namespace NAnalitycal {

using namespace NTypes;

class TAnalyticalSolution final {
public:
    TAnalyticalSolution(const TMatrix&);
    std::vector<Type> GetDistribution() const;

private:
    TVector CalculateStationaryDistribution() const;

private:
    const TMatrix P_;
    TVector Distribution_;
};

} // namespace NAnalitycal
