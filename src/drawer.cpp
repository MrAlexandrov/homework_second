#include "drawer.hpp"

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

namespace NDrawer {

constexpr double PenWidthArrowMultiplyer = 5;
constexpr double PenWidthNodeMultiplyer = 20;

void TDrawer::GenerateAndDrawGraph(const TMatrix& P, std::string filename, const std::vector<Type>& probability) {
    GenerateDotFile(P, filename, probability);
    GenerateImage(filename);
}

void TDrawer::GenerateDotFile(const TMatrix& P, std::string filename, const std::vector<Type>& probability) {
    size_t n = P.rows();
    std::ofstream dotFile(filename + ".dot");
    if (!dotFile.is_open()) {
        throw std::runtime_error("Failed to create DOT file: " + filename);
    }

    dotFile << "digraph MarkovChain {\n";
    dotFile << "    rankdir=LR;\n";
    dotFile << "    node [shape=circle];\n";

    dotFile << GenerateNodesPenwidthCommand(P, probability);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Type label = P(i, j);
            Type penwidth = label * PenWidthArrowMultiplyer;
            if (label > 0) {
                dotFile << GenerateArrow(i, j, label, penwidth);
            }
        }
    }

    dotFile << "}\n";
}

std::string TDrawer::GenerateArrow(int from, int to, Type label, Type penwidth) {
    std::stringstream command;
    command << "    "
            << "S" << from << " -> S" << to
            << " [label=\"" << std::fixed << std::setprecision(2)
            << label << "\", "
            << "penwidth=" << penwidth << "];\n";
    return command.str();
}

std::string TDrawer::GenerateNodesPenwidthCommand(const TMatrix& P, const std::vector<Type>& probability) {
    if (probability.empty()) return {};
    std::stringstream command;
    for (size_t i = 0, end = P.rows(); i < end; ++i) {
        command << "    S" << i
                << "[label=\"S"
                << i << "\", penwidth=" << std::max(0.1, probability[i] * PenWidthNodeMultiplyer)
                << "];\n";
    }
    return command.str();
}

void TDrawer::GenerateImage(std::string filename) {
    try {
        std::string command = "dot -Tpng " + filename + ".dot -o " + filename + ".png";
        if (system(command.c_str()) != 0) {
            throw std::runtime_error("Failed to generate PNG image");
        }

        std::cout << "PNG image created successfully: " << filename << ".png" << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}

} // namespace NDrawer
