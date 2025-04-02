#include "common.h"
#include "fourier.h"
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iomanip>
#include <iostream>
#include <vector>

void analyze_convergence(SpectralFourier<double> &solver,
                         const std::vector<int> &N_values,
                         AnalyticalFunction<double> ana_u,
                         AnalyticalFunction<double> ana_du) {
  std::cout << std::setw(10) << "N" << std::setw(20) << "L_inf Error"
            << std::setw(20) << "L_2 Error" << std::setw(15) << "L_inf Ratio"
            << std::setw(15) << "L_2 Ratio" << std::endl;
  std::cout << "---------------------------------------------------------------"
               "----------------"
            << std::endl;

  double L_inf_prev = 0;
  double L_2_prev = 0;

  for (size_t n_idx = 0; n_idx < N_values.size(); n_idx++) {
    int N = N_values[n_idx];

    auto [L_inf, L_2] = solver.compute_errors(N, ana_u, ana_du);

    double L_inf_ratio = 0, L_2_ratio = 0;
    if (n_idx > 0 && L_inf_prev > 0 && L_2_prev > 0) {
      L_inf_ratio = L_inf_prev / L_inf;
      L_2_ratio = L_2_prev / L_2;
    }

    std::cout << std::setw(10) << N << std::scientific << std::setprecision(4)
              << std::setw(20) << L_inf << std::setw(20) << L_2;

    if (n_idx > 0) {
      std::cout << std::fixed << std::setprecision(4) << std::setw(15)
                << L_inf_ratio << std::setw(15) << L_2_ratio;
    } else {
      std::cout << std::setw(15) << "---" << std::setw(15) << "---";
    }
    std::cout << std::endl;

    L_inf_prev = L_inf;
    L_2_prev = L_2;
  }

  std::cout << "---------------------------------------------------------------"
               "----------------"
            << std::endl
            << std::endl;
}

void print_summary(std::vector<int> &k_values,
                   std::vector<std::pair<int, mp_float>> &results) {

  std::cout << "------------------------------------------------------"
            << std::endl;
  std::cout << "k      | Minimum N | Max Relative Error" << std::endl;
  std::cout << "------------------------------------------------------"
            << std::endl;

  for (size_t i = 0; i < k_values.size(); i++) {
    std::cout << std::left << std::setw(7) << k_values[i] << "| "
              << std::setw(9) << results[i].first << " | " << std::scientific
              << std::setprecision(10) << results[i].second << std::endl;
  }
  std::cout << "------------------------------------------------------"
            << std::endl;
}

int main() {
  std::vector<int> k_values = {2, 4, 6, 8, 10, 12};
  mp_float error_threshold = 1e-5;
  int max_N = 100;

  std::cout << "ODD Fourier differentiation matrix accuracy for u(x) = "
               "exp(k*sin(x))\n";
  std::cout << "Error threshold: " << error_threshold << std::endl << std::endl;
  SpectralFourier<mp_float> solver_odd(MethodType::ODD);

  std::vector<std::pair<int, mp_float>> results = solver_odd.find_min_N(
      k_values, TestFunctions::func_ex01_u<mp_float>,
      TestFunctions::func_ex01_du<mp_float>, error_threshold, max_N);
  print_summary(k_values, results);

  std::cout << "EVEN Fourier differentiation matrix accuracy for u(x) = "
               "exp(k*sin(x))\n";
  std::cout << "Error threshold: " << error_threshold << std::endl << std::endl;

  SpectralFourier<mp_float> solver_even(MethodType::EVEN);
  results = solver_even.find_min_N(
      k_values, TestFunctions::func_ex01_u<mp_float>,
      TestFunctions::func_ex01_du<mp_float>, error_threshold, max_N);
  print_summary(k_values, results);

  std::vector<int> N_values = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};

  std::cout << "=================================================" << std::endl;
  std::cout << "Convergence Analysis for Fourier Differentiation" << std::endl;
  std::cout << "=================================================" << std::endl;
  std::cout << "\nMethod: EVEN" << std::endl;

  // This behavior is expected for Fourier methods with smooth periodic
  // functions.
  SpectralFourier<double> solver_even_double(MethodType::EVEN);
  analyze_convergence(solver_even_double, N_values,
                      TestFunctions::func_ex02_1_u<double>,
                      TestFunctions::func_ex02_1_du<double>);
  // For cos(x/2) on a [0,2π] domain, the function is not periodic (one full
  // period would require [0,4π]), which explains the poor performance.
  analyze_convergence(solver_even_double, N_values,
                      TestFunctions::func_ex02_2_u<double>,
                      TestFunctions::func_ex02_2_du<double>);
  // Fourier methods are known to perform poorly for non-periodic functions.
  analyze_convergence(solver_even_double, N_values,
                      TestFunctions::func_ex02_3_u<double>,
                      TestFunctions::func_ex02_3_du<double>);

  return 0;
}
