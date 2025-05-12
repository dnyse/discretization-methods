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

namespace BurgersExact {

template <NumericType T> T phi(T a, T b, T nu, int num_terms = 20) {
  T result = T(0);
  for (int k = -num_terms; k <= num_terms; ++k) {
    T arg = (a - T(2 * k + 1) * MathConstants<T>::PI());
    arg = arg * arg / (T(4) * nu * b);
    result += exp(-arg);
  }
  return result;
}

template <NumericType T> T dphi_dx(T a, T b, T nu, int num_terms = 20) {
  T result = T(0);
  for (int k = -num_terms; k <= num_terms; ++k) {
    T pi_term = T(2 * k + 1) * MathConstants<T>::PI();
    T arg = (a - pi_term) * (a - pi_term) / (T(4) * nu * b);
    // derivative of exp(-arg) with respect to x (through a = x - ct)
    // d/dx[exp(-arg)] = exp(-arg) * (-1) * 2*(a-pi_term)/(4*nu*b) * da/dx
    // da/dx = 1
    T coeff = -(a - pi_term) / (T(2) * nu * b);
    result += coeff * exp(-arg);
  }
  return result;
}

template <NumericType T> T exact_solution(T x, T t, T c = T(4), T nu = T(0.1)) {
  T a = x - c * t;
  T b = t + T(1);
  T phi_val = phi(a, b, nu);
  T dphi_dx_val = dphi_dx(a, b, nu);
  return c - T(2) * nu * dphi_dx_val / phi_val;
}

// Initial condition: exact solution at t=0
template <NumericType T> T initial_condition(T x, int k = 0) {
  return exact_solution(x, T(0));
}

// For the exact derivative with respect to x, we need to be careful with chain
// rule
template <NumericType T>
T exact_derivative(T x, T t, T c = T(4), T nu = T(0.1)) {
  T a = x - c * t;
  T b = t + T(1);

  // We need d/dx[c - 2ν(∂φ/∂x)/φ]
  // = -2ν * d/dx[(∂φ/∂x)/φ]
  // = -2ν * [(∂²φ/∂x²)*φ - (∂φ/∂x)²]/φ²

  T phi_val = phi(a, b, nu);
  T dphi_dx_val = dphi_dx(a, b, nu);

  // Second derivative of phi with respect to x
  T d2phi_dx2 = T(0);
  for (int k = -20; k <= 20; ++k) {
    T pi_term = T(2 * k + 1) * MathConstants<T>::PI();
    T arg = (a - pi_term) * (a - pi_term) / (T(4) * nu * b);

    // d²/dx²[exp(-arg)] requires chain rule twice
    T exp_arg = exp(-arg);
    T first_deriv_coeff = -(a - pi_term) / (T(2) * nu * b);

    // Second derivative term
    T second_deriv_coeff = -T(1) / (T(2) * nu * b);
    T squared_term = first_deriv_coeff * first_deriv_coeff;

    d2phi_dx2 += (second_deriv_coeff + squared_term) * exp_arg;
  }

  // Apply quotient rule: d/dx[(∂φ/∂x)/φ] = [(∂²φ/∂x²)*φ - (∂φ/∂x)²]/φ²
  T numerator = d2phi_dx2 * phi_val - dphi_dx_val * dphi_dx_val;
  T denominator = phi_val * phi_val;

  return -T(2) * nu * numerator / denominator;
}

// Wrapper functions for Differentiator interface
template <NumericType T> T burgers_initial_u(T x, int k = 0) {
  return initial_condition<T>(x);
}

template <NumericType T> T burgers_initial_du(T x, int k = 0) {
  return exact_derivative<T>(x, T(0));
}

} // namespace BurgersExact

#endif // INCLUDE_INCLUDE_COMMON_H_
