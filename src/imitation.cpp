#include "imitation.hpp"
#include "utils.hpp"

#include <vector>
#include <boost/asio.hpp>
#include <thread>

namespace NImitation {

using TThreadPool = boost::asio::thread_pool;

TImitationSolution::TImitationSolution(const TMatrix& P, size_t imitations, size_t iterations)
    : P_(P)
    , NumberStates_(P_.rows())
    , Imitations_(imitations)
    , Iterations_(iterations)
{
}

std::vector<std::vector<int>> TImitationSolution::GetAllStates() const {
    return States_;
}

std::vector<Type> TImitationSolution::ImitateAndGetDistribution() {
    ImitateSolution(Imitations_, Iterations_);
    CalculateDistribution();
    return GetDistribution();
}

void TImitationSolution::ImitateSolution(int imitations, int iterations) {
    TThreadPool pool(std::thread::hardware_concurrency());
    for (int startNode = 0; startNode < NumberStates_; ++startNode) {
        for (int _ = 0; _ < imitations; ++_) {
            boost::asio::post(pool,
                [this, startNode, iterations]() {
                    this->Imitation(startNode, iterations);
                }
            );
        }
    }
    pool.join();
}

void TImitationSolution::CalculateDistribution() {
    Distribution_.resize(NumberStates_, 0.0);
    for (const auto& imitation : States_) {
        for (const auto& state : imitation) {
            Distribution_[state] += 1;
        }
    }
    int total = NumberStates_ * Imitations_ * Iterations_;
    for (auto& state : Distribution_){
        state /= total;
    }
}

std::vector<Type> TImitationSolution::GetDistribution() const {
    return Distribution_;
}

void TImitationSolution::Imitation(int startNode, int iterations) {
    std::vector<int> localStates = ImitationImpl(startNode, iterations);

    {
        std::lock_guard<std::mutex> lock(StatesMutex_);
        States_.push_back(std::move(localStates));
    }
}

std::vector<int> TImitationSolution::ImitationImpl(int startNode, int iterations) const {
    std::vector<int> localStates;
    localStates.reserve(iterations);

    int currentState = startNode;
    for (int i = 0; i < iterations; ++i) {
        localStates.push_back(currentState);
        currentState = GetNextState(currentState);
    }

    return localStates;
}

int TImitationSolution::GetNextState(int currentState) const {
    auto randomNumber = NUtils::GenerateRandomNumber();
    Type currentValue = 0;
    for (int i = 0; i < NumberStates_; ++i) {
        currentValue += P_(currentState, i);
        if (randomNumber <= currentValue) {
            return i;
        }
    }
    assert(false);
}

} // namespace NImitation
