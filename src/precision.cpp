#include "precision.hpp"

#include <iomanip>
#include <iostream>

namespace NPrecision {

TPrecision::TPrecision(
    std::ios_base::fmtflags oldFlags,
    std::streamsize oldPrecision,
    size_t newPrecision
)
    : Flags_(oldFlags)
    , Precision_(oldPrecision)
{
    std::cout << std::setprecision(newPrecision) << std::fixed;
}

TPrecision::~TPrecision() {
    std::cout.flags(Flags_);
    std::cout.precision(Precision_);
}

} // namespace NPrecision
