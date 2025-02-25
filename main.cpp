#include "analytical.hpp"
#include "graph.hpp"
#include "drawer.hpp"
#include "precision.hpp"
#include "simulation.hpp"
#include "utils.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
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

template <typename T>
void PrintResults(std::string_view text, const std::vector<T>& results) {
    std::cout << text << "\n";
    for (const auto& i : results) {
        std::cout << i << " ";
    }
    std::cout << "\n";
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
    std::cout << N << "\n";
    TMatrix P = GetTransitionMatrix(N);
    std::cout << "P:\n" << P << "\n\n";

    NSimulation::TSimulationSolution totalSimulated(P, Imitations, Iterations);
    auto totalSimulatedDistribution = totalSimulated.GetDistribution();
    PrintResults("totalSimulatedDistribution:", totalSimulatedDistribution);

    NGraph::TGraph graph(P);
    graph.FillStronglyConnectedComponents();
    auto color = graph.GetColor();
    auto stronglyConnectedComponents = graph.GetStronglyConnectedComponents();
    graph.FillCondensation();
    auto condensatedGraph = graph.GetCondensation();

    TMatrix condensatedMatrix = GenerateMatrix(condensatedGraph);
    std::cout << "condensatedMatrix:\n" << condensatedMatrix << "\n\n";

    // Distribution, when we start in random node
    NGraph::TGraph condensated(condensatedMatrix);
    auto topologyCondensatedGraph = condensated.GetTopologicalSort();
    std::reverse(topologyCondensatedGraph.begin(), topologyCondensatedGraph.end());
    int amountStronglyConnectedComponents = condensatedMatrix.rows();
    std::vector<Type> probabilityRandomStart(amountStronglyConnectedComponents, 0);
    for (int i = 0; i < amountStronglyConnectedComponents; ++i) {
        probabilityRandomStart[i] = static_cast<double>(stronglyConnectedComponents[i].size()) / N;
    }
    PrintResults("probabilityRandomStart start in component:", probabilityRandomStart);

    for (int i = 0; i < amountStronglyConnectedComponents; ++i) {
        int from = topologyCondensatedGraph[i];
        for (int j = 0; j < amountStronglyConnectedComponents; ++j) {
            if (i == j) continue;
            int to = topologyCondensatedGraph[j];
            Type transition = condensatedMatrix(from, to);
            if (NUtils::Equals(transition, 0)) continue;
            probabilityRandomStart[to] += probabilityRandomStart[from] * transition;
        }
        if (!NUtils::Equals(condensatedMatrix(from, from), 1)) {
            probabilityRandomStart[from] = 0;
        }
    }
    PrintResults("probabilityRandomStart:", probabilityRandomStart);

    NDrawer::TDrawer::GenerateAndDrawGraph(condensatedMatrix, "condensated", probabilityRandomStart);

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

    std::vector<Type> definedStart(N, 0);
    std::vector<Type> randomStart(N, 0);

    for (size_t currentColor = 0, end = condensatedComponents.size(); currentColor < end; ++currentColor) {
        const auto& currentComponent = condensatedComponents[currentColor];
        int n = currentComponent.rows();
        // check if state is irrevocable
        if (!NUtils::Equals(condensatedMatrix(currentColor, currentColor), 1)) continue;
        NAnalitycal::TAnalyticalSolution analytic(currentComponent);
        NSimulation::TSimulationSolution imitated(currentComponent, Imitations, Iterations);
        std::vector<Type> analyticDistribution = analytic.GetDistribution();
        std::vector<Type> imitatedDistribution = imitated.GetDistribution();
        std::vector<Type> errors = CountError(analyticDistribution, imitatedDistribution);

        for (int i = 0; i < n; ++i) {
            const auto& globalNumber = stronglyConnectedComponents[currentColor][i];
            // distribution = <probability of being in component> * <probability of being at vertex>
            definedStart[globalNumber] = condensatedMatrix(0, currentColor) * analyticDistribution[i];
            randomStart[globalNumber] = probabilityRandomStart[currentColor] * analyticDistribution[i];
        }

        {
            std::cout << "Component " << currentColor << ":\n";
            PrintResults("Nodes:", stronglyConnectedComponents[currentColor]);

            std::cout << "currentComponent:\n";
            std::cout << currentComponent << "\n";
            NPrecision::TPrecision changer(
                std::cout.flags(),
                std::cout.precision(),
                6
            );

            PrintResults("Analytical distribution:", analyticDistribution);
            PrintResults("Imitated distribution:", imitatedDistribution);
            PrintResults("Errors:", errors);
        }
        std::cout << "\n";
    }
    NDrawer::TDrawer::GenerateAndDrawGraph(P, "chain", definedStart);

    PrintResults("definedStart:", definedStart);
    Type sumDefinedStart = std::accumulate(definedStart.begin(), definedStart.end(), 0.0);
    assert(NUtils::Equals(sumDefinedStart, 1));

    PrintResults("randomStart:", randomStart);
    Type sumRandomStart = std::accumulate(randomStart.begin(), randomStart.end(), 0.0);
    assert(NUtils::Equals(sumRandomStart, 1));

    auto distributionError = CountError(totalSimulatedDistribution, randomStart);
    PrintResults("distributionError:", distributionError);
    return 0;
}
