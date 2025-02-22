#include "analytical.hpp"
#include "graph.hpp"
#include "grapher.hpp"
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
    std::cout << "Determinant: " << P.determinant() << std::endl;
    std::cout << "Determinant (P - E): " << (P - TMatrix::Identity(N, N)).determinant() << std::endl;
    NGraph::TGraph graph(P);
    graph.FillStronglyConnectedComponents();
    auto components = graph.GetStronglyConnectedComponents();
    std::cout << "components:\n";
    for (const auto& i : components) {
        for (const auto& j : i) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
    graph.FillCondensation();
    auto condensation = graph.GetCondensation();
    auto color = graph.GetColor();
    std::cout << "condensation:\n";
    for (int i = 0, end = condensation.size(); i < end; ++i) {
        std::cout << i << ": ";
        for (const auto& [to, weight] : condensation[i]) {
            std::cout << to << ", " << weight << "\t\t";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
    TMatrix transitions = GenerateMatrix(condensation);
    std::cout << "transitions:\n";
    std::cout << transitions << std::endl;

    NAnalitycal::TAnalyticalSolution asdf(transitions);
    auto qwer = asdf.GetDistribution();


    std::vector<TMatrix> MatrixComponents;
    for (const auto& component : components) {
        size_t n = component.size();
        TMatrix comp(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                comp(i, j) = P(component[i], component[j]);
            }
        }
        MatrixComponents.emplace_back(comp);
    }

    std::vector<Type> totalDistribution(N, 0);

    for (size_t p = 0, end = MatrixComponents.size(); p < end; ++p) {
        const auto& current = MatrixComponents[p];
        int n = current.rows();
        bool ok = true;
        for (int i = 0; i < n; ++i) {
            Type sum = 0;
            for (int j = 0; j < n; ++j) {
                sum += current(i, j);
            }
            if (!NUtils::Equals(sum, 1)) ok = false;
        }
        if (!ok) continue;
        std::cout << "current:\n";
        std::cout << current << std::endl;

        NAnalitycal::TAnalyticalSolution analytic(current);
        NSimulation::TSimulationSolution simulate(current, Imitations, Iterations);
        std::vector<Type> analyticSolution = analytic.GetDistribution();
        std::vector<Type> simulateSolution = simulate.GetDistribution();
        std::vector<Type> errors = CountError(analyticSolution, simulateSolution);

        for (int i = 0; i < n; ++i) {
            totalDistribution[components[p][i]] = transitions(0, p) * analyticSolution[i];
        }

        {
            NPrecision::TPrecision changer(
                std::cout.flags(),
                std::cout.precision(),
                6
            );

            PrintResults("Analytical distribution:", analyticSolution);
            PrintResults("Imitated distribution:", simulateSolution);
            PrintResults("Errors:", errors);
        }
        std::cout << std::endl;
    }
    NGrapher::TGrapher drawing(P, "chain", totalDistribution);
    std::cout << "totalDistribution:\n";
    for (size_t i = 0, end = totalDistribution.size(); i < end; ++i) {
        std::cout << i << ": " << totalDistribution[i] << std::endl;
    }
    return 0;
}
