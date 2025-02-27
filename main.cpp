#include "analytical.hpp"
#include "graph.hpp"
#include "drawer.hpp"
#include "plotter.hpp"
#include "precision.hpp"
#include "imitation.hpp"
#include "utils.hpp"
#include "types.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <string_view>
#include <unordered_map>
#include <vector>

using namespace NTypes;

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

std::vector<Type> AbsoluteError(const std::vector<Type>& lhs, const std::vector<Type>& rhs) {
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

std::vector<Type> SampleAverage(const std::vector<std::vector<int>>& states) {
    size_t n = states.size();
    assert(n != 0);
    size_t m = states.front().size();
    std::vector<Type> average(m, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            average[j] += states[i][j];
        }
    }
    for (int j = 0; j < m; ++j) {
        average[j] /= n;
    }
    return average;
}

std::vector<Type> StandartDeviation(const std::vector<Type>& average, const std::vector<std::vector<double>>& probabilities) {
    size_t n = probabilities.size();
    size_t m = average.size();
    std::vector<Type> diffs(m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            Type value = average[j] - probabilities[i][j];
            diffs[j] += value * value;
        }
    }
    for (int j = 0; j < m; ++j) {
        diffs[j] /= n - 1;
        diffs[j] = sqrt(diffs[j]);
    }
    return diffs;
}

TMatrix GenerateMatrix(const std::vector<std::unordered_map<int, Type>>& other) {
    int n = other.size();

    TMatrix matrix(n, n);
    for (int from = 0; from < n; ++from) {
        Type sum = 0;
        for (const auto& [to, weight] : other[from]) {
            sum += weight;
            matrix(from, to) = weight;
        }
        for (const auto& [to, weight] : other[from]) {
            matrix(from, to) = weight / sum;
        }
    }
    return matrix;
}

TMatrix RenumberNodes(const TMatrix& input, const std::vector<int>& renumber) {
    size_t n = input.rows();
    assert(n == renumber.size());
    TMatrix P(n, n);
    for (size_t i = 0; i < n; ++i) {
        size_t newI = renumber[i];
        for (size_t j = 0; j < n; ++j) {
            size_t newJ = renumber[j];
            P(i, j) = input(newI, newJ);
        }
    }
    return P;
}

TMatrix ReorderNodes(const TMatrix& input) {
    size_t n = input.rows();
    NGraph::TGraph inputGraph(input);
    std::vector<int> renumberedNodes;
    {
        inputGraph.FillStronglyConnectedComponents();
        auto stronglyConnectedComponents = inputGraph.GetStronglyConnectedComponents();
        renumberedNodes.reserve(n);
        for (const auto& component : stronglyConnectedComponents) {
            for (const auto& node : component) {
                renumberedNodes.push_back(node);
            }
        }
    }
    return RenumberNodes(input, renumberedNodes);
}

std::vector<std::vector<int>> GenerateStates(
    int amountStronglyConnectedComponents
    , const std::vector<std::vector<int>>& stronglyConnectedComponents
    , const TMatrix& condensatedMatrix
    , const NImitation::TImitationSolution& totalImitated
    , int Iterations
) {
    std::vector<std::vector<int>> states;
    for (int currentColor = 0; currentColor < amountStronglyConnectedComponents; ++currentColor) {
        size_t numberNodes = stronglyConnectedComponents[currentColor].size();
        int amount = NUtils::Equals(condensatedMatrix(currentColor, currentColor), 1) ? 2 : 6;
        for (int i = 0; i < amount; ++i) {
            int startNode = NUtils::GenerateIntNumber(0, numberNodes - 1);
            startNode = stronglyConnectedComponents[currentColor][startNode];
            states.emplace_back(totalImitated.ImitationImpl(startNode, Iterations));
        }
    }
    return states;
}

void PlotCharts(int Iterations, const std::vector<std::vector<int>>& states) {
    NPlotter::TPlotter chart("chart");
    std::vector<int> XValues(Iterations);
    std::iota(XValues.begin(), XValues.end(), 0);
    chart.SetXValues(XValues);
    for (int i = 0, end = states.size(); i < end; ++i) {
        chart.EmplaceChart(states[i], std::to_string(i));
    }
    chart.Plot();
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage " << argv[0] << " <Imitations> <Iterations>";
        return 1;
    }
    size_t Imitations = std::stoul(argv[1]);
    size_t Iterations = std::stoul(argv[2]);

    size_t N = 0;
    std::cin >> N;
    TMatrix inputP = GetTransitionMatrix(N);

    TMatrix P = ReorderNodes(inputP);
    std::cout << "P:\n" << P << "\n\n";

    NImitation::TImitationSolution totalImitated(P, Imitations, Iterations);

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

    // probability to reach irrevocable states
    std::vector<Type> probabilityRandomStart(amountStronglyConnectedComponents, 0);
    {    
        for (int currentColor = 0; currentColor < amountStronglyConnectedComponents; ++currentColor) {
            probabilityRandomStart[currentColor] = static_cast<double>(stronglyConnectedComponents[currentColor].size()) / N;
        }
        PrintResults("probabilityRandomStart start in component:", probabilityRandomStart);

        for (int currentColor = 0; currentColor < amountStronglyConnectedComponents; ++currentColor) {
            int from = topologyCondensatedGraph[currentColor];
            for (int j = currentColor + 1; j < amountStronglyConnectedComponents; ++j) {
                int to = topologyCondensatedGraph[j];
                Type transition = condensatedMatrix(from, to);
                if (NUtils::Equals(transition, 0)) continue;
                probabilityRandomStart[to] += probabilityRandomStart[from] * (transition / (1.0 - condensatedMatrix(from, from)));
            }
            if (!NUtils::Equals(condensatedMatrix(from, from), 1)) {
                probabilityRandomStart[from] = 0;
            }
        }
        PrintResults("probabilityRandomStart in irrevocable state:", probabilityRandomStart);
    }
    NDrawer::TDrawer::GenerateAndDrawGraph(condensatedMatrix, "condensated", probabilityRandomStart);

    // creating condensated matrices
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
        int nodesInComponent = currentComponent.rows();
        // check if state is irrevocable
        if (!NUtils::Equals(condensatedMatrix(currentColor, currentColor), 1)) continue;
        NAnalitycal::TAnalyticalSolution analytic(currentComponent);
        NImitation::TImitationSolution imitated(currentComponent, Imitations, Iterations);
        std::vector<Type> analyticDistribution = analytic.CalculateAndGetDistribution();
        std::vector<Type> imitatedDistribution = imitated.ImitateAndGetDistribution();
        std::vector<Type> absoluteError = AbsoluteError(analyticDistribution, imitatedDistribution);

        // counting standart deviation
        auto states = imitated.GetAllStates();
        std::vector<std::vector<double>> probabilities(Imitations, std::vector(Iterations, 0.0));
        for (size_t i = 0; i < Imitations; ++i) {
            for (const auto& currentState : states[i]) {
                ++probabilities[i][currentState];
            }
            for (int j = 0; j < N; ++j) {
                probabilities[i][j] /= Iterations;
            }
        }
        std::vector<Type> standartDeviation = StandartDeviation(imitatedDistribution, probabilities);

        for (int i = 0; i < nodesInComponent; ++i) {
            const auto& globalNumber = stronglyConnectedComponents[currentColor][i];
            // distribution = <probability of being in component> * <probability of being at vertex>
            // definedStart[globalNumber] = (condensatedMatrix(0, currentColor) / (1.0 - condensatedMatrix(0, 0))) * analyticDistribution[i];
            randomStart[globalNumber] = probabilityRandomStart[currentColor] * analyticDistribution[i];
        }

        {
            std::cout << "Component " << currentColor << ":\n";
            PrintResults("Nodes:", stronglyConnectedComponents[currentColor]);

            // std::cout << "currentComponent:\n";
            // std::cout << currentComponent << "\n";
            NPrecision::TPrecision changer(
                std::cout.flags(),
                std::cout.precision(),
                6
            );

            PrintResults("Analytical distribution:", analyticDistribution);
            PrintResults("Imitated distribution:", imitatedDistribution);
            PrintResults("Absolute Error:", absoluteError);
            PrintResults("Standart Deviation:", standartDeviation);
        }
        std::cout << "\n";
    }
    NDrawer::TDrawer::GenerateAndDrawGraph(P, "chain", randomStart);

    {
        NPrecision::TPrecision changer(
            std::cout.flags(),
            std::cout.precision(),
            6
        );
        auto totalImitatedDistribution = totalImitated.ImitateAndGetDistribution();
        PrintResults("Imitated distribution:", totalImitatedDistribution);

        // PrintResults("definedStart:", definedStart);
        // Type sumDefinedStart = std::accumulate(definedStart.begin(), definedStart.end(), 0.0);
        // assert(NUtils::Equals(sumDefinedStart, 1));

        PrintResults("Analytic distribution:", randomStart);
        Type sumRandomStart = std::accumulate(randomStart.begin(), randomStart.end(), 0.0);
        assert(NUtils::Equals(sumRandomStart, 1));

        auto distributionError = AbsoluteError(totalImitatedDistribution, randomStart);
        PrintResults("distributionError:", distributionError);
    }

    std::vector<std::vector<int>> states = GenerateStates(
        amountStronglyConnectedComponents
        , stronglyConnectedComponents
        , condensatedMatrix
        , totalImitated
        , Iterations
    );
    PlotCharts(Iterations, states);
    return 0;
}
