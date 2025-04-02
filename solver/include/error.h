#ifndef INCLUDE_INCLUDE_ERROR_H_
#define INCLUDE_INCLUDE_ERROR_H_
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "common.h"

namespace error {
template <NumericType T> T relative(T &analy_du, T &numerical_du) {
  T denom = abs(analy_du) + std::numeric_limits<T>::epsilon();
  return abs(numerical_du - analy_du) / denom;
}
template <NumericType T> T absolute_diff(T &analy_du, T &numerical_du) {
  return abs(numerical_du - analy_du);
}

} // namespace error

#endif // INCLUDE_INCLUDE_ERROR_H_
