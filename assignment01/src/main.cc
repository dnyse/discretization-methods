#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

std::vector<std::vector<double>> build_fourier_D(int N) {
  std::vector<std::vector<double>> D(N + 1, std::vector<double>(N + 1, 0.0));

  for (int j = 0; j < N + 1; j++) {
    for (int i = 0; i < N + 1; i++) {
      if (i != j) {
        D[j][i] =
            (std::pow(-1, j + i)) / (2.0 * std::sin((j - i) * M_PI / (N + 1)));
      }
    }
  }

  return D;
}

double analytical_u(double x, int k) { return std::exp(k * std::sin(x)); }

double analytical_du(double x, int k) {
  return k * std::cos(x) * std::exp(k * std::sin(x));
}

double calculate_relative_error(int k, int N) {
  std::vector<double> x(N + 1);
  double step = 2.0 * M_PI / (N + 1);
  for (int i = 0; i < N + 1; i++) {
    x[i] = i * step;
  }

  std::vector<double> u(N + 1);
  std::vector<double> analy_du(N + 1);
  for (int i = 0; i < N + 1; i++) {
    u[i] = analytical_u(x[i], k);
    analy_du[i] = analytical_du(x[i], k);
  }

  std::vector<std::vector<double>> D = build_fourier_D(N);

  std::vector<double> numerical_du(N + 1, 0.0);
  for (int j = 0; j < N + 1; j++) {
    for (int i = 0; i < N + 1; i++) {
      numerical_du[j] += D[j][i] * u[i];
    }
  }

  std::vector<double> errors(N + 1);
  for (int i = 0; i < N + 1; i++) {
    double denom =
        std::abs(analy_du[i]) + std::numeric_limits<double>::epsilon();
    errors[i] = std::abs(numerical_du[i] - analy_du[i]) / denom;
    if (k == 12) {
      std::cout << "[" << N << "] Analy sol: " << analy_du[i]
                << " Num sol:" << numerical_du[i]
                << " diff: " << std::abs(numerical_du[i] - analy_du[i])
                << " denom: " << denom << " error: " << errors[i] << std::endl;
    }
  }

  double max_error = *std::max_element(errors.begin(), errors.end());
  if (k == 12)
    std::cout << "[" << N << "] Max error: " << max_error << std::endl;

  return max_error;
}

std::vector<std::pair<int, double>> find_min_N(const std::vector<int> &k_values,
                                               double error_threshold = 1e-5,
                                               int max_N = 100) {
  std::vector<std::pair<int, double>> results;

  for (int k : k_values) {
    int N = 5;
    double max_error = 1.0;
    int best_N = -1;
    double best_error = max_error;

    for (size_t N = 0; N < max_N + 1; N++) {
      max_error = calculate_relative_error(k, N);
      if (max_error < best_error) {
        best_error = max_error;
        best_N = N;
        if (best_error < error_threshold) {
          break;
        }
      }
    }
    results.push_back({best_N, best_error});
  }

return results;
}

int main() {
  std::vector<int> k_values = {2, 4, 6, 8, 10, 12};
  double error_threshold = 1e-5;
  int max_N = 100;

  std::cout << "Testing Fourier differentiation matrix accuracy for u(x) = "
               "exp(k*sin(x))"
            << std::endl;
  std::cout << "Error threshold: " << error_threshold << std::endl;

  std::vector<std::pair<int, double>> results =
      find_min_N(k_values, error_threshold, max_N);

  std::cout << "\nSummary:" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "k      | Minimum N | Max Relative Error" << std::endl;
  std::cout << "========================================" << std::endl;

  for (size_t i = 0; i < k_values.size(); i++) {
    std::cout << std::left << std::setw(7) << k_values[i] << "| "
              << std::setw(9) << results[i].first << " | " << std::scientific
              << std::setprecision(2) << results[i].second << std::endl;
  }

  std::cout << "========================================" << std::endl;

  return 0;
}
