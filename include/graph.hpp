#pragma once

#include "types.hpp"

#include <vector>

namespace NGraph {

using namespace NTypes;

class TGraph final {
private:
    struct TEdge {
        int To;
        double Weight;
    };

public:
    TGraph(const TMatrix& P);

    void FillStronglyConnectedComponents();
    std::vector<std::vector<int>> GetStronglyConnectedComponents() const;

    void FillCondensation();
    std::vector<std::unordered_map<int, Type>> GetCondensation() const;

    std::vector<int> GetColor() const;

private:
    void TopologicalSort(int from, std::vector<bool>& used, std::vector<int>& output) const;
    void FillComponent(int from, std::vector<bool>& used);

private:
    int Nodes_;
    std::vector<std::vector<TEdge>> Graph_;
    std::vector<std::vector<TEdge>> ReversedGraph_;
    std::vector<int> Color_;
    std::vector<std::vector<int>> Components_;
    std::vector<std::unordered_map<int, Type>> Condensation_;
};

} // namespace NGraph
