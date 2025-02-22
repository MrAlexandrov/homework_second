#pragma once

#include <Eigen/Dense>

namespace NTypes {

using Type = double;
using TVector = Eigen::Matrix<Type, Eigen::Dynamic, 1>;
using TMatrix = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;

} // namespace NTypes
