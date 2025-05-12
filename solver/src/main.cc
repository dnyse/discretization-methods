#include "../include/burger.h"
#include "../include/common.h"
#include <boost/math/constants/constants.hpp>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <vector>

// Function to run Burgers' equation convergence study using Fourier spectral
// method
void burgers_convergence_study() {
  std::vector<int> N_values = {16, 32, 48, 64, 96, 128, 192, 256};
  double t_final = MathConstants<double>::PI() / 4.0; // t = π/4
  double nu = 0.1;
  double c = 4.0;

  std::vector<double> errors(N_values.size());
  std::vector<double> times(N_values.size());

  std::cout
      << "\nBurgers' Equation Convergence Study (Fourier Spectral Method - ODD)"
      << std::endl;
  std::cout << "Parameters: ν = " << nu << ", c = " << c << ", t_final = π/4"
            << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "L_inf Error"
            << std::setw(20) << "Computation Time" << std::endl;
  std::cout << std::string(48, '-') << std::endl;

  for (size_t i = 0; i < N_values.size(); ++i) {
    int N = N_values[i];
    std::cout << std::setw(8) << N;

    // Create Fourier spectral solver with ODD method
    auto spectral = std::make_shared<SpectralFourier<double>>(MethodType::ODD);
    spectral->build(N);

    // Initialize and solve
    BurgersSolver<double> solver(spectral, nu);
    solver.initialize(N, t_final);
    solver.solve();

    // Get results
    auto [x, u_numerical, u_exact, error, time] = solver.get_results();
    errors[i] = error;
    times[i] = time;

    std::cout << std::scientific << std::setprecision(6) << std::setw(20)
              << error << std::fixed << std::setprecision(4) << std::setw(20)
              << time << "s" << std::endl;
  }

  // Calculate convergence rates
  std::cout << "\nConvergence Rates:" << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "Convergence Rate"
            << std::endl;
  std::cout << std::string(28, '-') << std::endl;

  for (size_t i = 1; i < N_values.size(); ++i) {
    double ratio = double(N_values[i]) / double(N_values[i - 1]);
    double rate = log(errors[i - 1] / errors[i]) / log(ratio);

    std::cout << std::setw(8) << N_values[i] << std::fixed
              << std::setprecision(2) << std::setw(20) << rate << std::endl;
  }
}

// Function to show time evolution using Fourier spectral method
void burgers_time_evolution() {
  int N = 128; // As specified in the exam
  std::vector<double> times = {
      0.0,
      boost::math::constants::pi<double>() / 8.0,  // t = π/8
      boost::math::constants::pi<double>() / 6.0,  // t = π/6
      boost::math::constants::pi<double>() / 4.0}; // t = π/4
  double nu = 0.1;
  double c = 4.0;

  std::cout << "\nBurgers' Equation Time Evolution (N = " << N
            << ", Fourier Spectral - ODD)" << std::endl;
  std::cout << "Parameters: ν = " << nu << ", c = " << c << std::endl;
  std::cout << std::setw(15) << "Time" << std::setw(20) << "L∞ Error"
            << std::setw(20) << "Computation Time" << std::endl;
  std::cout << std::string(55, '-') << std::endl;

  // Create Fourier spectral solver
  auto spectral = std::make_shared<SpectralFourier<double>>(MethodType::ODD);
  spectral->build(N);

  for (double t_final : times) {
    BurgersSolver<double> solver(spectral, nu);
    solver.initialize(N, t_final);
    solver.solve();

    auto [x, u_numerical, u_exact, error, time] = solver.get_results();

    std::cout << std::fixed << std::setprecision(4) << std::setw(15) << t_final
              << std::scientific << std::setprecision(6) << std::setw(20)
              << error << std::fixed << std::setprecision(4) << std::setw(20)
              << time << "s" << std::endl;

    // Could add plotting here if needed
    // plot_solution(x, u_numerical, u_exact, t_final);
  }
}

int main(int argc, char *argv[]) {
  // Print OpenMP status
  std::cout << "\n=== OpenMP Status ===" << std::endl;
  std::cout << "OpenMP Version: " << _OPENMP << std::endl;
  std::cout << "Number of available threads: " << omp_get_max_threads()
            << std::endl;

  std::cout << "\n=== Exercise 3 Part 2: Burgers' Equation ===" << std::endl;
  std::cout << "Using Fourier Spectral Method (ODD)" << std::endl;

  // Part (b): Convergence study
  std::cout << "\nPart (b): Convergence Study" << std::endl;
  burgers_convergence_study();

  // Part (d): Time evolution plots
  std::cout << "\n\nPart (d): Time Evolution" << std::endl;
  burgers_time_evolution();

  return 0;
}
