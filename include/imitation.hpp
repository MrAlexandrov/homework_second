#pragma once

#include "types.hpp"

#include <vector>
#include <boost/asio.hpp>

namespace NImitation {

using namespace NTypes;

using TThreadPool = boost::asio::thread_pool;

class TImitationSolution final {
public:
    TImitationSolution(const TMatrix& P, size_t imitations, size_t iterations);

    std::vector<Type> ImitateAndGetDistribution();
    std::vector<std::vector<int>> GetAllStates() const;

private:
    void CalculateDistribution();
    std::vector<Type> GetDistribution() const;

    void ImitateSolution(int imitations, int iterations);
    void Imitation(int startNode, int iterations);
    std::vector<int> ImitationImpl(int startNode, int iterations) const;

    int GetNextState(int currentState) const;

private:
    const TMatrix P_;
    size_t NumberStates_;
    std::vector<std::vector<int>> States_;
    std::vector<Type> Distribution_;
    int Imitations_;
    int Iterations_;
    mutable std::mutex StatesMutex_;
};

} // namespace NImitation
