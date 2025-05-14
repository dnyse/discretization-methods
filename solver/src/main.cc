#include "../include/burger.h"
#include "../include/common.h"
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

// Helper function to check if an argument is present
bool contains_arg(int argc, char *argv[], const std::string &arg) {
  return std::find(argv + 1, argv + argc, arg) != argv + argc;
}

// Print usage information
void print_usage(const char *prog_name) {
  std::cout << "\nUsage: " << prog_name << " [options]\n\n";
  std::cout << "Options:\n";
  std::cout << "  --help       Show this help message\n";
  std::cout << "  --all        Run all exercises\n";
  std::cout << "\nExercise 3 Part 2: Burgers' Equation (Collocation)\n";
  std::cout << "  --ex02b      Run Part (b) - Convergence study\n";
  std::cout << "  --ex02d      Run Part (d) - Time evolution\n";
  std::cout << "\nExercise 3 Part 3: Burgers' Equation (Galerkin)\n";
  std::cout << "  --ex03b      Run Part (b) - CFL determination\n";
  std::cout << "  --ex03c      Run Part (c) - Convergence study\n";
  std::cout << "  --ex03d      Run Part (d) - Comparison and plots\n";
  std::cout << std::endl;
}

// Part 2(b): Burgers convergence study (Collocation)
void burgers_find_max_cfl() {
    std::vector<int> N_values = {16, 32, 48, 64, 96, 128, 192, 256};
    double nu = 0.1;
    
    std::cout << "\nPart 2(b): Finding Maximum CFL Values" << std::endl;
    std::cout << "Parameters: ν = " << nu << std::endl;
    std::cout << std::setw(8) << "N" << std::setw(20) << "Max CFL" << std::endl;
    std::cout << std::string(28, '-') << std::endl;
    
    std::vector<double> max_cfls;
    
    for (int N : N_values) {
        auto spectral = std::make_shared<SpectralFourier<double>>(MethodType::ODD);
        spectral->build(N);
        
        BurgersSolver<double> solver(spectral, nu);
        double max_cfl = solver.find_max_cfl(N);
        max_cfls.push_back(max_cfl);
        
        std::cout << std::setw(8) << N << std::fixed << std::setprecision(3)
                  << std::setw(20) << max_cfl << std::endl;
    }
    
    // Analyze if CFL definition is reasonable
    std::cout << "\nAnalysis: Does the CFL definition seem reasonable?" << std::endl;
    std::cout << "The CFL values should typically decrease as N increases for spectral methods." << std::endl;
    std::cout << "Observed trend: ";
    for (size_t i = 1; i < max_cfls.size(); ++i) {
        if (max_cfls[i] < max_cfls[i-1]) {
            std::cout << "↓ ";
        } else {
            std::cout << "↑ ";
        }
    }
    std::cout << std::endl;
}
// Part 2(d): Burgers time evolution (Collocation)
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

// Part 3(b): CFL determination (Galerkin)
void determine_cfl_values() {
  std::vector<int> N_values = {16, 32, 48, 64, 96, 128, 192, 256};
  double nu = 0.1;

  std::cout << "\nPart 3(b): Determining Maximum CFL Values for "
               "Fourier-Galerkin Method"
            << std::endl;
  std::cout << "Parameters: ν = " << nu << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "Max CFL" << std::endl;
  std::cout << std::string(28, '-') << std::endl;

  std::vector<double> cfl_values;

  for (int N : N_values) {
    auto galerkin = std::make_shared<FourierGalerkin<double>>();
    BurgersGalerkinSolver<double> solver(galerkin, nu);
    solver.initialize(N, 1.0); // Initialize for CFL testing

    double max_cfl = solver.find_cfl_limit();
    cfl_values.push_back(max_cfl);

    std::cout << std::setw(8) << N << std::fixed << std::setprecision(3)
              << std::setw(20) << max_cfl << std::endl;
  }
}

//  Part 3(c): Convergence study (Galerkin)
void galerkin_convergence_study() {
  std::vector<int> N_values = {16, 32, 48, 64, 96, 128, 192, 256};
  double t_final = MathConstants<double>::PI() / 4.0;
  double nu = 0.1;

  std::vector<double> errors(N_values.size());
  std::vector<double> times(N_values.size());

  std::cout << "\nPart 3(c): Fourier-Galerkin Convergence Study at t = π/4"
            << std::endl;
  std::cout << "Parameters: ν = " << nu << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "L∞ Error"
            << std::setw(20) << "Computation Time" << std::endl;
  std::cout << std::string(48, '-') << std::endl;

  for (size_t i = 0; i < N_values.size(); ++i) {
    int N = N_values[i];
    std::cout << std::setw(8) << N;

    auto galerkin = std::make_shared<FourierGalerkin<double>>();
    BurgersGalerkinSolver<double> solver(galerkin, nu);
    solver.initialize(N, t_final);
    solver.solve();

    auto [x, u_numerical, u_exact, error, time, cfl] = solver.get_results();
    errors[i] = error;
    times[i] = time;

    std::cout << std::scientific << std::setprecision(6) << std::setw(20)
              << error << std::fixed << std::setprecision(4) << std::setw(20)
              << time << "s" << std::endl;
  }

  // Calculate convergence rates
  std::cout << "\nConvergence Rates:" << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "Rate" << std::endl;
  std::cout << std::string(28, '-') << std::endl;

  for (size_t i = 1; i < N_values.size(); ++i) {
    double ratio = double(N_values[i]) / double(N_values[i - 1]);
    double rate = log(errors[i - 1] / errors[i]) / log(ratio);

    std::cout << std::setw(8) << N_values[i] << std::fixed
              << std::setprecision(2) << std::setw(20) << rate << std::endl;
  }
}

// Part 3(d): Comparison and plots (Galerkin)
void compare_with_collocation() {
  std::vector<int> N_values = {64, 128};
  double t_final = MathConstants<double>::PI() / 4.0;
  double nu = 0.1;

  std::cout
      << "\nPart 3(d): Comparison of Fourier-Galerkin vs Fourier-Collocation"
      << std::endl;
  std::cout << "Parameters: ν = " << nu << ", t = π/4" << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "Galerkin Error"
            << std::setw(20) << "Collocation Error" << std::endl;
  std::cout << std::string(48, '-') << std::endl;

  for (int N : N_values) {
    // Galerkin method
    auto galerkin = std::make_shared<FourierGalerkin<double>>();
    BurgersGalerkinSolver<double> galerkin_solver(galerkin, nu);
    galerkin_solver.initialize(N, t_final);
    galerkin_solver.solve();
    auto [x_g, u_g, u_exact_g, error_g, time_g, cfl_g] =
        galerkin_solver.get_results();

    // Collocation method (using existing SpectralFourier)
    auto collocation =
        std::make_shared<SpectralFourier<double>>(MethodType::ODD);
    collocation->build(N);
    BurgersSolver<double> collocation_solver(collocation, nu);
    collocation_solver.initialize(N, t_final);
    collocation_solver.solve();
    auto [x_c, u_c, u_exact_c, error_c, time_c] =
        collocation_solver.get_results();

    std::cout << std::setw(8) << N << std::scientific << std::setprecision(6)
              << std::setw(20) << error_g << std::setw(20) << error_c
              << std::endl;
  }
}

int main(int argc, char *argv[]) {
  // Print OpenMP status
  std::cout << "\n=== OpenMP Status ===" << std::endl;
  std::cout << "OpenMP Version: " << _OPENMP << std::endl;
  std::cout << "Number of available threads: " << omp_get_max_threads()
            << std::endl;

  if (argc < 2 || contains_arg(argc, argv, "--help")) {
    print_usage(argv[0]);
    return argc < 2 ? 1 : 0;
  }

  bool run_all = contains_arg(argc, argv, "--all");
  bool run_ex02b = run_all || contains_arg(argc, argv, "--ex02b");
  bool run_ex02d = run_all || contains_arg(argc, argv, "--ex02d");
  bool run_ex03b = run_all || contains_arg(argc, argv, "--ex03b");
  bool run_ex03c = run_all || contains_arg(argc, argv, "--ex03c");
  bool run_ex03d = run_all || contains_arg(argc, argv, "--ex03d");

  if (!run_ex02b && !run_ex02d && !run_ex03b && !run_ex03c && !run_ex03d) {
    std::cerr << "Error: No valid exercise options specified." << std::endl;
    print_usage(argv[0]);
    return 1;
  }

  if (run_ex02b) {
    std::cout << "\n=== Exercise 3 Part 2: Burgers' Equation ===" << std::endl;
    std::cout << "Using Fourier Spectral Method (ODD)" << std::endl;
    std::cout << "\nPart (b): Convergence Study" << std::endl;
    burgers_convergence_study();
  }

  if (run_ex02d) {
    std::cout << "\n\nPart (d): Time Evolution" << std::endl;
    burgers_time_evolution();
  }

  if (run_ex03b) {
    std::cout << "\n=== Exercise 3 Part 3: Burgers' Equation with "
                 "Fourier-Galerkin Method ==="
              << std::endl;
    determine_cfl_values();
  }

  if (run_ex03c) {
    galerkin_convergence_study();
  }

  if (run_ex03d) {
    compare_with_collocation();
  }

  return 0;
}
