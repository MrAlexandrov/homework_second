#include "grapher.hpp"

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

namespace NGrapher {

constexpr double PenWidthArrowMultiplyer = 5;
constexpr double PenWidthNodeMultiplyer = 20;

TGrapher::TGrapher(
    const TMatrix& P,
    std::string filename = "default",
    const std::vector<Type>& probability = {}
)
    : P_(P)
    , Filename_(filename)
    , Probability_(probability)
{
    GenerateDotFile();
    GenerateImage();
}

void TGrapher::GenerateDotFile() const {
    size_t n = P_.rows();
    std::ofstream dotFile(Filename_ + ".dot");
    if (!dotFile.is_open()) {
        throw std::runtime_error("Failed to create DOT file: " + Filename_);
    }

    dotFile << "digraph MarkovChain {\n";
    dotFile << "    rankdir=LR;\n";
    dotFile << "    node [shape=circle];\n";

    dotFile << GenerateNodesPenwidthCommand();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Type label = P_(i, j);
            Type penwidth = label * PenWidthArrowMultiplyer;
            if (label > 0) {
                dotFile << GenerateArrow(i, j, label, penwidth);
            }
        }
    }

    dotFile << "}\n";
}

std::string TGrapher::GenerateArrow(int from, int to, Type label, Type penwidth) const {
    std::stringstream command;
    command << "    "
            << "S" << from << " -> S" << to
            << " [label=\"" << std::fixed << std::setprecision(2)
            << label << "\", "
            << "penwidth=" << penwidth << "];\n";
    return command.str();
}

std::string TGrapher::GenerateNodesPenwidthCommand() const {
    if (Probability_.empty()) return {};
    std::stringstream command;
    for (size_t i = 0, end = P_.rows(); i < end; ++i) {
        command << "    S" << i 
                << "[label=\"S" 
                << i << "\", penwidth=" << std::max(0.1, Probability_[i] * PenWidthNodeMultiplyer)
                << "];\n";
    }
    return command.str();
}

void TGrapher::GenerateImage() const {
    try {
        std::string command = "dot -Tpng " + Filename_ + ".dot -o " + Filename_ + ".png";
        if (system(command.c_str()) != 0) {
            throw std::runtime_error("Failed to generate PNG image");
        }

        std::cout << "PNG image created successfully: " << Filename_ << ".png" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

} // namespace NGrapher
