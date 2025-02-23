#pragma once

#include "types.hpp"

#include <string>
#include <vector>

namespace NDrawer {

using namespace NTypes;

class TDrawer{
public:
    TDrawer(const TMatrix& P, std::string filename, const std::vector<Type>& probability);

private:
    void GenerateDotFile() const;
    std::string GenerateNodesPenwidthCommand() const;
    std::string GenerateArrow(int from, int to, Type label, Type penwidth) const;
    void GenerateImage() const;

private:
    TMatrix P_;
    std::string Filename_;
    std::vector<Type> Probability_;
};

} // namespace NDrawer
