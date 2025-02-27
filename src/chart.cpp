#include "chart.hpp"

#include <string>

namespace NChart {

void TChart::SetTitle(const std::string& title) {
    Title_ = title;
}

std::vector<Type> TChart::GetData() const {
    return Data_;
}

std::string TChart::GetTitle() const {
    return Title_;
}

size_t TChart::GetSize() const {
    return Data_.size();
}

} // namespace NChart
