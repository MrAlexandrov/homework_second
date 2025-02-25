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

    std::vector<Type> GetDistribution() const;

private:
    void ImitateSolution(int imitations, int iterations);

    void Imitation(int iterations);

    int GetNextState(int currentState) const;

private:
    const TMatrix P_;
    size_t NumberStates_;
    std::vector<std::atomic<int>> Count_;
    int Imitations_;
    int Iterations_;
};

} // namespace NImitation
