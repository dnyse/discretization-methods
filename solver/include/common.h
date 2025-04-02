#ifndef INCLUDE_INCLUDE_COMMON_H_
#define INCLUDE_INCLUDE_COMMON_H_

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

// Define constants and types with template support
using mp_float = boost::multiprecision::cpp_dec_float_50;

template <typename T>
concept NumericType = std::is_same_v<T, mp_float> || std::is_same_v<T, double>;

template <NumericType T>
struct MathConstants {
    static T PI() { return boost::math::constants::pi<T>(); }
};

// Template for analytical functions
template <NumericType T>
using AnalyticalFunction = T (*)(T, int);

// Enum stays the same
enum MethodType { EVEN, ODD };

// Template for test functions
namespace TestFunctions {
    template <NumericType T>
    T func_ex01_u(T x, int k) { return exp(T(k) * sin(x)); }
    
    template <NumericType T>
    T func_ex02_1_u(T x, int k = 0) { return cos(T(10) * x); }
    
    template <NumericType T>
    T func_ex02_2_u(T x, int k = 0) { return cos(x / T(2)); }
    
    template <NumericType T>
    T func_ex02_3_u(T x, int k = 0) { return x; }
    
    template <NumericType T>
    T func_ex01_du(T x, int k) {
        return T(k) * cos(x) * exp(T(k) * sin(x));
    }
    
    template <NumericType T>
    T func_ex02_1_du(T x, int k = 0) { return -T(10) * sin(T(10) * x); }
    
    template <NumericType T>
    T func_ex02_2_du(T x, int k = 0) { return -T(0.5) * sin(x / T(2)); }
    
    template <NumericType T>
    T func_ex02_3_du(T x, int k = 0) { return T(1); }
}

namespace util {}

#endif // INCLUDE_INCLUDE_COMMON_H_
