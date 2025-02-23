#include "graph.hpp"
#include "utils.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <iostream>
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
    #ifdef DEBUG
    std::cout << "Color:\n";
    for (const auto& i : Color_) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    #endif // DEBUG
    return Condensation_;
}

void TGraph::FillCondensation() {
    Condensation_.resize(Components_.size());
    for (int from = 0; from < Nodes_; ++from) {
        #ifdef DEBUG
        std::cout << "from: " << from << std::endl;
        #endif // DEBUG
        for (const auto& [to, weight] : Graph_[from]) {
            auto fromComponent = Color_[from];
            auto toComponent = Color_[to];
            #ifdef DEBUG
            std::cout << "fromComponent, toComponent: " << fromComponent << " " << toComponent << std::endl;
            std::cout << "weight: " << weight << std::endl;
            #endif // DEBUG
            if (NUtils::Equals(weight, 0) || fromComponent == toComponent) continue;
            Condensation_[fromComponent][toComponent] += weight;
        }
    }
}

void TGraph::FillStronglyConnectedComponents() {
    std::vector<bool> used(Nodes_, false);
    std::vector<int> output;

    for (int i = 0; i < Nodes_; ++i) {
        if (used[i]) continue;
        TopologicalSort(i, used, output);
    }
    std::reverse(output.begin(), output.end());

    used.assign(Nodes_, false);
    for (const auto& to : output) {
        if (used[to]) continue;
        Components_.push_back({});
        FillComponent(to, used);
    }
}

std::vector<std::vector<int>> TGraph::GetStronglyConnectedComponents() const {
    return Components_;
}

std::vector<int> TGraph::GetColor() const {
    return Color_;
}

void TGraph::TopologicalSort(int from, std::vector<bool>& used, std::vector<int>& output) const {
    used[from] = true;
    for (const auto&[to, weight] : Graph_[from]) {
        if (used[to]) continue;
        TopologicalSort(to, used, output);
    }
    output.push_back(from);
}

void TGraph::FillComponent(int from, std::vector<bool>& used) {
    used[from] = true;
    Components_.back().push_back(from);
    Color_[from] = static_cast<int>(Components_.size()) - 1;
    for (const auto& [to, weight] : ReversedGraph_[from]) {
        if (used[to]) continue;
        FillComponent(to, used);
    }
}

} // namespace NGraph
