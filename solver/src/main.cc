#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iomanip>
#include <iostream>
#include <vector>

// Define high precision type with 50 digits of precision
using mp_float = boost::multiprecision::cpp_dec_float_50;

// Constants
const mp_float PI = boost::math::constants::pi<mp_float>();
typedef mp_float (*AnalyticalFunction)(mp_float, int);

enum MethodType { EVEN, ODD };

std::vector<std::vector<mp_float>> build_fourier_D(int N, MethodType method) {
  int grid_pts = method == MethodType::ODD ? N + 1 : N;
  std::vector<std::vector<mp_float>> D(grid_pts,
                                       std::vector<mp_float>(grid_pts, 0));
  if (method == MethodType::ODD) {
    for (int j = 0; j < grid_pts; j++) {
      for (int i = 0; i < grid_pts; i++) {
        if (i != j) {
          D[j][i] = (pow(mp_float(-1), j + i)) /
                    (mp_float(2) * sin((j - i) * PI / (grid_pts)));
        }
      }
    }
  } else {
    for (int j = 0; j < grid_pts; j++) {
      for (int i = 0; i < grid_pts; i++) {
        if (i != j) {
          D[j][i] = (pow(mp_float(-1), j + i)) /
                    (mp_float(2) * tan((j - i) * PI / (grid_pts)));
        }
      }
    }
  }
  return D;
}

mp_float func_ex01_u(mp_float x, int k) { return exp(mp_float(k) * sin(x)); }
mp_float func_ex02_1_u(mp_float x, int k = 0) { return cos(10 * x); }
mp_float func_ex02_2_u(mp_float x, int k = 0) { return cos(x / 2); }
mp_float func_ex02_3_u(mp_float x, int k = 0) { return x; }

mp_float func_ex01_du(mp_float x, int k) {
  return mp_float(k) * cos(x) * exp(mp_float(k) * sin(x));
}
mp_float func_ex02_1_du(mp_float x, int k = 0) { return -10 * sin(10 * x); }
mp_float func_ex02_2_du(mp_float x, int k = 0) { return -0.5 * sin(x / 2); }
mp_float func_ex02_3_du(mp_float x, int k = 0) { return 1; }

mp_float calculate_relative_error(int k, int N, MethodType method,
                                  AnalyticalFunction ana_u,
                                  AnalyticalFunction ana_du) {
  int grid_pts = method == MethodType::ODD ? N + 1 : N;
  std::vector<mp_float> x(grid_pts);
  mp_float step = mp_float(2) * PI / (grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    x[i] = i * step;
  }

  std::vector<mp_float> u(x.size()), analy_du(grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    u[i] = ana_u(x[i], k);
    analy_du[i] = ana_du(x[i], k);
  }

  std::vector<std::vector<mp_float>> D = build_fourier_D(N, method);
  std::vector<mp_float> numerical_du(grid_pts, 0);
  for (int j = 0; j < grid_pts; j++) {
    for (int i = 0; i < grid_pts; i++) {
      numerical_du[j] += D[j][i] * u[i];
    }
  }

  std::vector<mp_float> errors(grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    mp_float denom = abs(analy_du[i]) + std::numeric_limits<double>::epsilon();
    errors[i] = abs(numerical_du[i] - analy_du[i]) / denom;
  }

  return *std::max_element(errors.begin(), errors.end());
}

std::vector<std::pair<int, mp_float>>
find_min_N(const std::vector<int> &k_values, AnalyticalFunction ana_u,
           AnalyticalFunction ana_du, mp_float error_threshold = 1e-5,
           int max_N = 100, MethodType method = MethodType::ODD) {
  std::vector<std::pair<int, mp_float>> results;

  for (int k : k_values) {
    std::cout << "Finding N with max error < 10^-5 for k=" << k << std::endl;
    int best_N = -1;
    mp_float best_error = 1.0;

    for (int N = 4; N <= max_N; N += 2) {
      mp_float max_error =
          calculate_relative_error(k, N, method, ana_u, ana_du);
      // std::cout << "Found best N=" << N << " for k=" << k
      // << " error=" << max_error << std::endl;
      if (max_error < best_error) {
        best_error = max_error;
        best_N = N;
        if (best_error < error_threshold) {
          std::cout << "Found best N=" << N << " for k=" << k
                    << " with error=" << best_error << std::endl
                    << std::endl;
          break;
        }
      }
    }
    results.emplace_back(best_N, best_error);
  }
  return results;
}

std::tuple<mp_float, mp_float> calculate_errors(int N, MethodType method,
                                                AnalyticalFunction ana_u,
                                                AnalyticalFunction ana_du) {
  int grid_pts = method == MethodType::ODD ? N + 1 : N;
  std::vector<mp_float> x(grid_pts);
  mp_float step = mp_float(2) * PI / (grid_pts);

  // Set up grid points
  for (int i = 0; i < grid_pts; i++) {
    x[i] = i * step;
  }

  // Compute analytical solution and its derivative
  std::vector<mp_float> u(grid_pts), analy_du(grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    u[i] = ana_u(x[i], 0);
    analy_du[i] = ana_du(x[i], 0);
  }

  // Build differentiation matrix and compute numerical derivative
  std::vector<std::vector<mp_float>> D = build_fourier_D(N, method);
  std::vector<mp_float> numerical_du(grid_pts, 0);
  for (int j = 0; j < grid_pts; j++) {
    for (int i = 0; i < grid_pts; i++) {
      numerical_du[j] += D[j][i] * u[i];
    }
  }

  // Calculate pointwise relative errors
  std::vector<mp_float> rel_errors(grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    // mp_float denom = abs(analy_du[i]) +
    // std::numeric_limits<double>::epsilon(); rel_errors[i] =
    // abs(numerical_du[i] - analy_du[i]) / denom;
    rel_errors[i] = abs(numerical_du[i] - analy_du[i]);
  }

  // Calculate absolute errors for L2 norm
  std::vector<mp_float> abs_errors(grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    abs_errors[i] = abs(numerical_du[i] - analy_du[i]);
  }

  // Compute L_infinity norm (maximum relative error)
  mp_float L_inf = *std::max_element(rel_errors.begin(), rel_errors.end());

  // Compute L_2 norm (global error)
  mp_float L_2 = 0;
  for (int i = 0; i < grid_pts; i++) {
    L_2 += abs_errors[i] * abs_errors[i];
  }
  L_2 = sqrt(L_2 / grid_pts); // Normalized by number of points

  return std::make_tuple(L_inf, L_2);
}

void analyze_convergence(const std::vector<int> &N_values, MethodType method,
                         AnalyticalFunction ana_u, AnalyticalFunction ana_du) {
  std::cout << std::setw(10) << "N" << std::setw(20) << "L_inf Error"
            << std::setw(20) << "L_2 Error" << std::setw(15) << "L_inf Ratio"
            << std::setw(15) << "L_2 Ratio" << std::endl;
  std::cout << "---------------------------------------------------------------"
               "----------------"
            << std::endl;

  mp_float L_inf_prev = 0;
  mp_float L_2_prev = 0;

  for (size_t n_idx = 0; n_idx < N_values.size(); n_idx++) {
    int N = N_values[n_idx];

    auto [L_inf, L_2] = calculate_errors(N, method, ana_u, ana_du);

    mp_float L_inf_ratio = 0, L_2_ratio = 0;
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

  std::vector<std::pair<int, mp_float>> results =
      find_min_N(k_values, func_ex01_u, func_ex01_du, error_threshold, max_N,
                 MethodType::ODD);
  print_summary(k_values, results);

  std::cout << "EVEN Fourier differentiation matrix accuracy for u(x) = "
               "exp(k*sin(x))\n";
  std::cout << "Error threshold: " << error_threshold << std::endl << std::endl;

  results = find_min_N(k_values, func_ex01_u, func_ex01_du, error_threshold,
                       max_N, MethodType::EVEN);
  print_summary(k_values, results);

  std::vector<int> N_values = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};

  std::cout << "=================================================" << std::endl;
  std::cout << "Convergence Analysis for Fourier Differentiation" << std::endl;
  std::cout << "=================================================" << std::endl;
  std::cout << "\nMethod: EVEN" << std::endl;

  // This behavior is expected for Fourier methods with smooth periodic
  // functions.
  analyze_convergence(N_values, MethodType::EVEN, func_ex02_1_u,
                      func_ex02_1_du);
  // For cos(x/2) on a [0,2π] domain, the function is not periodic (one full
  // period would require [0,4π]), which explains the poor performance.
  analyze_convergence(N_values, MethodType::EVEN, func_ex02_2_u,
                      func_ex02_2_du);
  // Fourier methods are known to perform poorly for non-periodic functions.
  analyze_convergence(N_values, MethodType::EVEN, func_ex02_3_u,
                      func_ex02_3_du);

  return 0;
}
