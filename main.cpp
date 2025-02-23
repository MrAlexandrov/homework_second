#include "analytical.hpp"
#include "graph.hpp"
#include "drawer.hpp"
#include "precision.hpp"
#include "simulation.hpp"
#include "utils.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <string_view>
#include <unordered_map>
#include <vector>

using Type = double;
using TMatrix = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;

TMatrix GetTransitionMatrix(int n) {
    TMatrix P(n, n);
    for (int i = 0; i < n; ++i) {
        Type sum{0};
        for (int j = 0; j < n; ++j) {
            Type x;
            std::cin >> x;
            sum += x;
            P(i, j) = x;
        }
        assert(NUtils::Equals(sum, 1));
    }
    return P;
}

void PrintResults(std::string_view text, const std::vector<Type>& results) {
    std::cout << text << "\n";
    for (const auto& i : results) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

std::vector<Type> CountError(const std::vector<Type>& lhs, const std::vector<Type>& rhs) {
    assert(lhs.size() == rhs.size());
    int n = lhs.size();
    std::vector<Type> result(n);
    std::transform(
        lhs.cbegin(), lhs.cend(),
        rhs.cbegin(), result.begin(),
        [](Type lhs, Type rhs) {
            return std::fabs(lhs - rhs);
        }
    );
    return result;
}

void NormalizeMatrix(TMatrix& matrix) {
    int rows = matrix.rows();
    int cols = matrix.cols();
    for (int i = 0; i < rows; ++i) {
        Type sum = 0;
        for (int j = 0; j < cols; ++j) {
            sum += matrix(i, j);
        }
        if (NUtils::Equals(sum, 1)) continue;
        for (int j = 0; j < cols; ++j) {
            if (NUtils::Equals(matrix(i, j), 0)) continue;
            matrix(i, j) /= sum;
        }
        if (NUtils::Equals(sum, 0)) {
            matrix(i, i) = 1;
        }
    }
}

TMatrix GenerateMatrix(const std::vector<std::unordered_map<int, Type>>& other) {
    int n = other.size();

    TMatrix matrix(n, n);
    for (int from = 0; from < n; ++from) {
        for (const auto& [to, weight] : other[from]) {
            matrix(from, to) = weight;
        }
    }
    NormalizeMatrix(matrix);
    return matrix;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage " << argv[0] << " <Imitations> <Iterations>";
        return 1;
    }
    size_t Imitations = std::stoul(argv[1]);
    size_t Iterations = std::stoul(argv[2]);

    size_t N = 5;
    std::cin >> N;
    std::cout << N << std::endl;
    TMatrix P = GetTransitionMatrix(N);
    #ifdef DEBUG
    std::cout << "Determinant: " << P.determinant() << std::endl;
    std::cout << "Determinant (P - E): " << (P - TMatrix::Identity(N, N)).determinant() << std::endl;
    #endif // DEBUG
    NGraph::TGraph graph(P);
    graph.FillStronglyConnectedComponents();
    auto color = graph.GetColor();
    auto stronglyConnectedComponents = graph.GetStronglyConnectedComponents();
    #ifdef DEBUG
    std::cout << "stronglyConnectedComponents:\n";
    for (const auto& i : stronglyConnectedComponents) {
        for (const auto& j : i) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
    #endif // DEBUG
    graph.FillCondensation();
    auto condensatedGraph = graph.GetCondensation();
    #ifdef DEBUG
    std::cout << "condensatedGraph:\n";
    for (int i = 0, end = condensatedGraph.size(); i < end; ++i) {
        std::cout << i << ": ";
        for (const auto& [to, weight] : condensatedGraph[i]) {
            std::cout << to << ", " << weight << "\t\t";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
    #endif // DEBUG
    TMatrix condensatedMatrix = GenerateMatrix(condensatedGraph);

    std::cout << "condensatedMatrix:\n";
    std::cout << condensatedMatrix << "\n\n";

    NAnalitycal::TAnalyticalSolution condensatedSolution(condensatedMatrix);
    auto condensationDistribution = condensatedSolution.GetDistribution();

    std::vector<TMatrix> condensatedComponents;
    for (const auto& component : stronglyConnectedComponents) {
        size_t n = component.size();
        TMatrix comp(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                comp(i, j) = P(component[i], component[j]);
            }
        }
        condensatedComponents.emplace_back(comp);
    }

    std::vector<Type> totalDistribution(N, 0);

    for (size_t currentColor = 0, end = condensatedComponents.size(); currentColor < end; ++currentColor) {
        const auto& currentComponent = condensatedComponents[currentColor];
        int n = currentComponent.rows();
        if (!NUtils::Equals(condensatedMatrix(currentColor, currentColor), 1)) continue;
        #ifdef DEBUG
        std::cout << "currentComponent:\n";
        std::cout << currentComponent << std::endl;
        #endif // DEBUG
        NAnalitycal::TAnalyticalSolution analytic(currentComponent);
        NSimulation::TSimulationSolution imitated(currentComponent, Imitations, Iterations);
        std::vector<Type> analyticDistribution = analytic.GetDistribution();
        std::vector<Type> imitatedDistribution = imitated.GetDistribution();
        std::vector<Type> errors = CountError(analyticDistribution, imitatedDistribution);

        for (int i = 0; i < n; ++i) {
            const auto& globalNumber = stronglyConnectedComponents[currentColor][i];
            // distribution = <probability of being in component> * <probability of being at vertex>
            totalDistribution[globalNumber] = condensatedMatrix(0, currentColor) * analyticDistribution[i];
        }

        {
            std::cout << "Component " << currentColor << ":\n";
            std::cout << "Nodes: ";
            for (const auto& node : stronglyConnectedComponents[currentColor]) {
                std::cout << node << " ";
            }
            std::cout << "\n";
            NPrecision::TPrecision changer(
                std::cout.flags(),
                std::cout.precision(),
                6
            );

            PrintResults("Analytical distribution:", analyticDistribution);
            PrintResults("Imitated distribution:", imitatedDistribution);
            PrintResults("Errors:", errors);
        }
        std::cout << std::endl;
    }
    NDrawer::TDrawer drawing(P, "chain", totalDistribution);
    std::cout << "totalDistribution:\n";
    for (size_t i = 0, end = totalDistribution.size(); i < end; ++i) {
        std::cout << i << ": " << totalDistribution[i] << std::endl;
    }
    return 0;
}
