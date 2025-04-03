#ifndef INCLUDE_INCLUDE_COMMON_H_
#define INCLUDE_INCLUDE_COMMON_H_

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

// Define constants and types with template support
using mp_float = boost::multiprecision::cpp_dec_float_50;

template <typename T>
concept NumericType = std::is_same_v<T, mp_float> || std::is_same_v<T, double>;

template <NumericType T> struct MathConstants {
  static T PI() { return boost::math::constants::pi<T>(); }
};

// Template for analytical functions
template <NumericType T> using AnalyticalFunction = T (*)(T, int);
template <NumericType T> using HyperbolicFunction = T (*)(T, T);

// Enum stays the same
enum MethodType { EVEN, ODD };

// Template for test functions
namespace TestFunctions {
template <NumericType T> T func_ex01_u(T x, int k) {
  return exp(T(k) * sin(x));
}

template <NumericType T> T func_ex02_1_u(T x, int k = 0) {
  return cos(T(10) * x);
}

template <NumericType T> T func_ex02_2_u(T x, int k = 0) {
  return cos(x / T(2));
}

template <NumericType T> T func_ex02_3_u(T x, int k = 0) { return x; }

template <NumericType T> T func_ex01_du(T x, int k) {
  return T(k) * cos(x) * exp(T(k) * sin(x));
}

template <NumericType T> T func_ex02_1_du(T x, int k = 0) {
  return -T(10) * sin(T(10) * x);
}

template <NumericType T> T func_ex02_2_du(T x, int k = 0) {
  return -T(0.5) * sin(x / T(2));
}

template <NumericType T> T func_ex02_3_du(T x, int k = 0) { return T(1); }

template <NumericType T> T hyperbolic_exact_u(T x, T t) {
  return exp(sin(x - T(2) * MathConstants<T>::PI() * t));
}

// Initial condition u(x,0) = exp(sin(x))
template <NumericType T> T hyperbolic_initial_u(T x, int k = 0) {
  return exp(sin(x));
}

// Analytical derivative du/dx = cos(x-2πt) * exp(sin(x-2πt))
template <NumericType T> T hyperbolic_analytical_du(T x, T t) {
  T arg = x - T(2) * MathConstants<T>::PI() * t;
  return cos(arg) * exp(sin(arg));
}

// Wrapper for analytical derivative at t=0 (for Differentiator interface)
template <NumericType T> T hyperbolic_du_t0(T x, int k = 0) {
  return hyperbolic_analytical_du(x, T(0));
}

} // namespace TestFunctions

namespace util {}

#endif // INCLUDE_INCLUDE_COMMON_H_
