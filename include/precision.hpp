#pragma once

#include <iostream>

namespace NPrecision {

class TPrecision final {
public:
    TPrecision(std::ios_base::fmtflags oldFlags, std::streamsize oldPrecision, size_t newPrecision);

    ~TPrecision();
private:
    std::ios_base::fmtflags Flags_;
    std::streamsize Precision_;
};

} // namespace NPrecision
