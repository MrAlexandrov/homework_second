#pragma once

#include "types.hpp"

#include <string>

namespace NGrapher {

using namespace NTypes;

class TGrapher{
public:
    TGrapher(const TMatrix& P, std::string filename);

private:
    void GenerateDotFile() const;
    void GenerateImage() const;

private:
    TMatrix P_;
    std::string Filename_;
};

} // namespace NGrapher
