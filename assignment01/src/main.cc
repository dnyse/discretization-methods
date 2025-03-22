#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <matplot/matplot.h>
#include <numbers>
#include <vector>

using namespace matplot;
std::vector<double> linspace(double min, double max, int n) {
  std::vector<double> result(n);
  if (n == 1) {
    result[0] = min;
    return result;
  }
  double step = (max - min) / (n - 1);
  for (int i = 0; i < n; i++) {
    result[i] = min + i * step;
  }
  return result;
}

double u(double x, int k) { return std::exp(k * std::sin(x)); }
double du(double x, int k) {
  return k * std::cos(x) * std::exp(k * std::sin(x));
}

double D_ij(int i, int j, int N) {
  if (i == j)
    return 0.;

  return -std::pow(-1, j + i) / 2 * 1 /
         std::sin(((j - i) * std::numbers::pi) / (N + 1));
}

void create_fourier_mat(std::vector<std::vector<double>> &mat, int N) {
  // Correctly resize the matrix to be N x N
  mat.resize(N);
  for (size_t i = 0; i < N; i++) {
    mat[i].resize(N); // Resize inner vectors first
    for (size_t j = 0; j < N; j++) {
      mat[i][j] = D_ij(i, j, N); // Directly assign instead of push_back
    }
  }
}

std::vector<double>
compute_numerical_derivative(const std::vector<double> &u_values,
                             const std::vector<std::vector<double>> &D) {
  int N = u_values.size();
  std::vector<double> du_numerical(N, 0.0);

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      du_numerical[j] += D[j][i] * u_values[i];
    }
  }

  return du_numerical;
}

double calculate_max_relative_error(const std::vector<double> &numerical,
                                    const std::vector<double> &analytical,
                                    bool debug = false) {
  double max_error = 0.0;

  for (size_t i = 0; i < numerical.size(); i++) {
    double rel_error = std::abs((numerical[i] - analytical[i]) / analytical[i]);
    max_error = std::max(max_error, rel_error);
    if (debug) {
      double test_error = std::abs(numerical[i] - analytical[i]);
      std::cerr << "DEBUGPRINT[33]: main.cc:67: max_error=" << test_error
                << std::endl;
    }
  }

  return max_error;
}

int find_minimum_N(int k, double target_error = 1e-5, int max_N = 200) {
  for (int N = 5; N <= max_N; N++) {
    // Create uniformly spaced points in [0, 2π]
    std::vector<double> x_points = linspace(0, 2 * std::numbers::pi, N);

    // Compute function values at these points
    std::vector<double> u_values(N);
    std::vector<double> du_analytical_values(N);
    for (int i = 0; i < N; i++) {
      u_values[i] = u(x_points[i], k);
      du_analytical_values[i] = du(x_points[i], k);
    }

    // Create the differentiation matrix and compute numerical derivative
    std::vector<std::vector<double>> D;
    create_fourier_mat(D, N);
    auto du_numerical = compute_numerical_derivative(u_values, D);

    // Calculate the maximum relative error
    double max_error;
    if (k == 8 && N == 200) {
      max_error = calculate_max_relative_error(du_numerical,
                                               du_analytical_values, true);
      auto fig = figure(true);
      fig->backend()->run_command("unset warnings");
      fig->quiet_mode(true);
      hold(on);

      // Plot analytical solution (continuous line)
      plot(x_points, du_analytical_values, "b-")
          ->line_width(2)
          .display_name("Analytical");

      // Plot numerical solution (markers)
      plot(x_points, du_numerical, "ro")
          ->marker_size(3)
          .display_name("Numerical");

      // Add grid, labels, and title
      grid(on);
      xlabel("x");
      ylabel("du/dx");
      title("Comparison of Analytical vs Numerical Derivative for k = " +
            std::to_string(k) + ", N = " + std::to_string(N));
      ::matplot::legend(on);
      save("test.png");

    } else
      max_error =
          calculate_max_relative_error(du_numerical, du_analytical_values);

    if (max_error < target_error) {
      return N;
    }
  }
  return 0;
}

void plot_comparison(int k, int N) {

  // Create a dense grid for plotting the analytical solution smoothly
  int plot_points = 500;
  std::vector<double> x_dense = linspace(0, 2 * std::numbers::pi, plot_points);
  std::vector<double> du_analytical_dense(plot_points);

  for (int i = 0; i < plot_points; i++) {
    du_analytical_dense[i] = du(x_dense[i], k);
  }

  // Create grid for numerical solution with N points
  std::vector<double> x_points = linspace(0, 2 * std::numbers::pi, N);

  // Compute function values at these points
  std::vector<double> u_values(N);
  std::vector<double> du_analytical_values(N);

  for (int i = 0; i < N; i++) {
    u_values[i] = u(x_points[i], k);
    du_analytical_values[i] = du(x_points[i], k);
  }

  // Create the differentiation matrix and compute numerical derivative
  std::vector<std::vector<double>> D;
  create_fourier_mat(D, N);
  auto du_numerical = compute_numerical_derivative(u_values, D);

  // Create the plot
  auto fig = figure(true);
  fig->backend()->run_command("unset warnings");
  fig->quiet_mode(true);
  hold(on);

  // Plot analytical solution (continuous line)
  plot(x_dense, du_analytical_dense, "b-")
      ->line_width(2)
      .display_name("Analytical");

  // Plot numerical solution (markers)
  plot(x_points, du_numerical, "ro")->marker_size(8).display_name("Numerical");

  // Add grid, labels, and title
  grid(on);
  xlabel("x");
  ylabel("du/dx");
  title("Comparison of Analytical vs Numerical Derivative for k = " +
        std::to_string(k) + ", N = " + std::to_string(N));
  ::matplot::legend(on);

  // Compute the error for display
  double max_error =
      calculate_max_relative_error(du_numerical, du_analytical_values);
  text(0.5, 0.05, "Max Relative Error: " + std::to_string(max_error));

  // Save the figure
  save("fourier_comparison_k" + std::to_string(k) + "_N" + std::to_string(N) +
       ".png");
  hold(off);
}

void plot_error_convergence(int k, int max_N = 100) {
  using namespace matplot;

  std::vector<int> N_values;
  std::vector<double> errors;

  // Compute errors for different values of N
  for (int N = 5; N <= max_N; N += 5) {
    // Create uniformly spaced points in [0, 2π]
    std::vector<double> x_points = linspace(0, 2 * std::numbers::pi, N);

    // Compute function values at these points
    std::vector<double> u_values(N);
    std::vector<double> du_analytical_values(N);

    for (int i = 0; i < N; i++) {
      u_values[i] = u(x_points[i], k);
      du_analytical_values[i] = du(x_points[i], k);
    }

    // Create the differentiation matrix and compute numerical derivative
    std::vector<std::vector<double>> D;
    create_fourier_mat(D, N);
    auto du_numerical = compute_numerical_derivative(u_values, D);

    // Calculate the maximum relative error
    double max_error =
        calculate_max_relative_error(du_numerical, du_analytical_values);

    N_values.push_back(N);
    errors.push_back(max_error);
  }

  // Plot error convergence
  figure();
  semilogy(N_values, errors, "b-o")->line_width(2);
  grid(on);
  xlabel("N (grid points)");
  ylabel("Maximum Relative Error");
  title("Error Convergence for k = " + std::to_string(k));

  // Add horizontal line for target error 1e-5
  std::vector<int> x_line = {N_values.front(), N_values.back()};
  std::vector<double> y_line = {1e-5, 1e-5};
  auto horizontalLine = plot(x_line, y_line, "r-");
  horizontalLine->line_width(1.5);
  horizontalLine->display_name("Target Error (1e-5)");
  ::matplot::legend();

  // Save the figure
  save("error_convergence_k" + std::to_string(k) + ".png");
}

int main(int argc, char *argv[]) {
  int arr_k[6] = {2, 4, 6, 8, 10, 12};
  double target_error = 1e-5;

  std::cout << "Testing Fourier differentiation matrix accuracy on u(x) = "
               "exp(k * sin(x))"
            << std::endl;
  std::cout << "Target maximum relative error: " << target_error << std::endl;
  std::cout << "-------------------------------------------------------"
            << std::endl;
  std::cout << "| k value | Minimum N | Maximum Relative Error |" << std::endl;
  std::cout << "-------------------------------------------------------"
            << std::endl;

  // Table output for minimum N values
  for (int k : arr_k) {
    int min_N = find_minimum_N(k, target_error);

    // Compute the actual error for verification
    std::vector<double> x_points = linspace(0, 2 * std::numbers::pi, min_N);

    std::vector<double> u_values(min_N);
    std::vector<double> du_analytical_values(min_N);
    for (int i = 0; i < min_N; i++) {
      u_values[i] = u(x_points[i], k);
      du_analytical_values[i] = du(x_points[i], k);
    }

    std::vector<std::vector<double>> D;
    create_fourier_mat(D, min_N);
    auto du_numerical = compute_numerical_derivative(u_values, D);
    double actual_error =
        calculate_max_relative_error(du_numerical, du_analytical_values);

    std::cout << "| " << std::setw(7) << k << " | " << std::setw(9) << min_N
              << " | " << std::setw(23) << std::scientific
              << std::setprecision(6) << actual_error << " |" << std::endl;

    // Create plots for this k value
    // plot_comparison(k, 200);
    // plot_error_convergence(k);
  }

  std::cout << "-------------------------------------------------------"
            << std::endl;
  std::cout << "Plots have been saved as PNG files in the current directory."
            << std::endl;

  return 0;
}
