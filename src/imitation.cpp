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
    , Count_(NumberStates_)
    , Imitations_(imitations)
    , Iterations_(iterations)
{
    ImitateSolution(Imitations_, Iterations_);
}

std::vector<Type> TImitationSolution::GetDistribution() const {
    std::vector<Type> result;
    result.reserve(Count_.size());
    int total = Imitations_ * Iterations_;
    for (const auto& i : Count_) {
        result.emplace_back(static_cast<Type>(i) / total);
    }
    return result;
}

void TImitationSolution::ImitateSolution(int imitations, int iterations) {
    TThreadPool pool(std::thread::hardware_concurrency());
    for (int i = 0; i < imitations; ++i) {
        boost::asio::post(pool,
            [this, iterations]() {
                this->Imitation(iterations);
            }
        );
    }
    pool.join();
}

void TImitationSolution::Imitation(int iterations) {
    // int currentState = NUtils::GenerateIntNumber(0, NumberStates_ - 1);
    for (
        int i = 0, currentState = NUtils::GenerateIntNumber(0, NumberStates_ - 1);
        i < iterations;
        ++i, currentState = GetNextState(currentState)
    ) {
        Count_[currentState].fetch_add(1);
    }
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
