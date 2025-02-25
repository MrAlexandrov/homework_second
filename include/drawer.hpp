#pragma once

#include "types.hpp"

#include <string>
#include <vector>

namespace NDrawer {

using namespace NTypes;

class TDrawer{
private:
    TDrawer() = delete;

public:
    static void GenerateAndDrawGraph(const TMatrix& P, std::string filename, const std::vector<Type>& probability = {});

private:
    static void GenerateDotFile(const TMatrix& P, std::string filename, const std::vector<Type>& probability);
    static std::string GenerateNodesPenwidthCommand(const TMatrix& P, const std::vector<Type>& probability);
    static std::string GenerateArrow(int from, int to, Type label, Type penwidth);
    static void GenerateImage(std::string filename);
};

} // namespace NDrawer
