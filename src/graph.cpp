#include "graph.hpp"
#include "utils.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <unordered_map>

namespace NGraph {

TGraph::TGraph(const TMatrix& P)
    : Nodes_(P.rows())
    , Graph_(Nodes_)
    , ReversedGraph_(Nodes_)
    , Color_(Nodes_)
{
    for (int i = 0; i < Nodes_; ++i) {
        for (int j = 0; j < Nodes_; ++j) {
            if (NUtils::Equals(P(i, j), 0)) continue;
            Graph_[i].emplace_back(j, P(i, j));
            ReversedGraph_[j].emplace_back(i, P(i, j));
        }
    }
}

std::vector<std::unordered_map<int, Type>> TGraph::GetCondensation() const {
    return Condensation_;
}

void TGraph::FillCondensation() {
    Condensation_.resize(Components_.size());
    for (int from = 0; from < Nodes_; ++from) {
        for (const auto& [to, weight] : Graph_[from]) {
            auto fromComponent = Color_[from];
            auto toComponent = Color_[to];
            if (NUtils::Equals(weight, 0)) continue;
            Condensation_[fromComponent][toComponent] += weight;
        }
    }
}

void TGraph::FillStronglyConnectedComponents() {
    std::vector<int> output = GetTopologicalSort();
    std::reverse(output.begin(), output.end());

    std::vector<bool> used(Nodes_);
    for (const auto& to : output) {
        if (used[to]) continue;
        Components_.push_back({});
        FillComponent(to, used);
    }
    FillColors();
}

std::vector<std::vector<int>> TGraph::GetStronglyConnectedComponents() const {
    return Components_;
}

std::vector<int> TGraph::GetColor() const {
    return Color_;
}

std::vector<int> TGraph::GetTopologicalSort() const {
    std::vector<bool> used(Nodes_, false);
    std::vector<int> output;

    for (int i = 0; i < Nodes_; ++i) {
        if (used[i]) continue;
        TopologicalSortImpl(i, used, output);
    }
    return output;
}

void TGraph::TopologicalSortImpl(int from, std::vector<bool>& used, std::vector<int>& output) const {
    used[from] = true;
    for (const auto&[to, _] : Graph_[from]) {
        if (used[to]) continue;
        TopologicalSortImpl(to, used, output);
    }
    output.push_back(from);
}

void TGraph::FillComponent(int from, std::vector<bool>& used) {
    used[from] = true;
    Components_.back().push_back(from);
    for (const auto& [to, _] : ReversedGraph_[from]) {
        if (used[to]) continue;
        FillComponent(to, used);
    }
}

void TGraph::FillColors() {
    for (int currentColor = 0, end = Components_.size(); currentColor < end; ++currentColor) {
        for (const auto& element : Components_[currentColor]) {
            Color_[element] = currentColor;
        }
    }
}

} // namespace NGraph
