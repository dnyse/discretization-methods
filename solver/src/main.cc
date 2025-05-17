#include "../include/burger.h"
#include "../include/common.h"
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <iomanip>
#include <iostream>
#include <matplot/matplot.h> // For plotting
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
  std::cout << "  --ex02b      Run Part (b) - Find maximum CFL values\n";
  std::cout << "  --ex02c      Run Part (c) - Convergence study\n";
  std::cout << "  --ex02d      Run Part (d) - Time evolution\n";
  std::cout << "\nExercise 3 Part 3: Burgers' Equation (Galerkin)\n";
  std::cout << "  --ex03b      Run Part (b) - CFL determination\n";
  std::cout << "  --ex03c      Run Part (c) - Convergence study\n";
  std::cout << "  --ex03d      Run Part (d) - Comparison and plots\n";
  std::cout << std::endl;
}

// Part 2(b): Find maximum CFL values experimentally
std::vector<double> burgers_find_max_cfl() {
  std::vector<int> N_values = {16, 32, 48, 64, 96, 128, 192, 256};
  std::vector<double> max_cfls;
  double nu = 0.1;

  std::cout << "\n=== Part 2(b): Determining Maximum CFL Values ==="
            << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "Max CFL" << std::endl;
  std::cout << std::string(28, '-') << std::endl;

  for (int N : N_values) {
    auto spectral = std::make_shared<SpectralFourier<double>>(MethodType::ODD);
    spectral->build(N);

    BurgersSolver<double> solver(spectral, nu);
    double max_cfl = solver.find_max_cfl(N);
    max_cfls.push_back(max_cfl);

    std::cout << std::setw(8) << N << std::fixed << std::setprecision(4)
              << std::setw(20) << max_cfl << std::endl;
  }
  return max_cfls;
}

// Part 2(c): Convergence study using CFL values from Part 2(b)
void burgers_convergence_study_with_cfl(const std::vector<double> &cfl_values) {
  std::vector<int> N_values = {16, 32, 48, 64, 96, 128, 192, 256};
  double t_final = MathConstants<double>::PI() / 4.0; // t = π/4
  double nu = 0.1;

  std::cout << "\n=== Part 2(c): Convergence Study at t = π/4 ===" << std::endl;
  std::cout << "Using CFL values determined in Part 2(b)" << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(15) << "CFL" << std::setw(20)
            << "L∞ Error" << std::setw(20) << "CPU Time" << std::endl;
  std::cout << std::string(63, '-') << std::endl;

  std::vector<double> errors;

  for (size_t i = 0; i < N_values.size(); ++i) {
    int N = N_values[i];
    double cfl = cfl_values[i] * .2;
    // double cfl = 0.5;

    auto spectral = std::make_shared<SpectralFourier<double>>(MethodType::ODD);
    spectral->build(N);

    BurgersSolver<double> solver(spectral, nu);
    solver.initialize(N, t_final, cfl); // Use the CFL from Part 2(b)
    solver.solve();

    auto [x, u_numerical, u_exact, error, time] = solver.get_results();
    errors.push_back(error);

    using namespace matplot;

    auto fig = figure(true);
    fig->quiet_mode(true);
    fig->backend()->run_command("unset warnings");
    fig->size(1280, 960);
    fig->font_size(18);

    auto p2 = plot(x, u_numerical, "b-o");
    p2->line_width(3);
    p2->display_name("Numerical Solution");

    hold(true);

    auto p1 = plot(x, u_exact, "r");
    p1->line_width(3);
    p1->display_name("Analytical Solution");

    hold(false);

    auto ax = gca();
    ax->font_size(18); // Apply font size to axis-level elements
    ax->xlabel("x");
    ax->ylabel("u(x,t)");
    ax->title("Solution for Burgers Problem at tfinal with N=" +
              std::to_string(N_values[i]) + ", CFL=" + std::to_string(cfl));

    legend()->font_size(16);

    save("burger_tfinal_" + std::to_string(N_values[i]) + ".png");

    std::cout << std::setw(8) << N << std::fixed << std::setprecision(4)
              << std::setw(15) << cfl << std::scientific << std::setprecision(6)
              << std::setw(20) << error << std::fixed << std::setprecision(4)
              << std::setw(20) << time << "s" << std::endl;
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

  std::cout
      << "\nDoes the convergence rate correspond to what you would expect?"
      << std::endl;
  std::cout << "Spectral methods typically show exponential convergence for "
               "smooth periodic solutions."
            << std::endl;
  std::cout << "The Burgers' sawtooth solution has steep gradients, which may "
               "affect convergence."
            << std::endl;
}

// Part 2(d): Burgers time evolution
void burgers_time_evolution(const std::vector<double> &cfl_values) {
  int N = 128; // As specified in the exam
  std::vector<double> times = {
      0.0,
      boost::math::constants::pi<double>() / 8.0,  // t = π/8
      boost::math::constants::pi<double>() / 6.0,  // t = π/6
      boost::math::constants::pi<double>() / 4.0}; // t = π/4
  double nu = 0.1;

  // Find the CFL for N=128 from our results
  std::vector<int> N_values = {16, 32, 48, 64, 96, 128, 192, 256};
  double cfl_128 = 0.5; // default
  for (size_t i = 0; i < N_values.size(); ++i) {
    if (N_values[i] == N) {
      cfl_128 = cfl_values[i];
      break;
    }
  }

  std::cout << "\n=== Part 2(d): Time Evolution for N = " << N
            << " ===" << std::endl;
  std::cout << "Parameters: ν = " << nu << ", CFL = " << cfl_128 << std::endl;
  std::cout << std::setw(15) << "Time" << std::setw(20) << "L∞ Error"
            << std::setw(20) << "CPU Time" << std::endl;
  std::cout << std::string(55, '-') << std::endl;

  // Create Fourier spectral solver
  auto spectral = std::make_shared<SpectralFourier<double>>(MethodType::ODD);
  spectral->build(N);

  for (double t_final : times) {
    BurgersSolver<double> solver(spectral, nu);
    solver.initialize(N, t_final, cfl_128);
    solver.solve();

    auto [x, u_numerical, u_exact, error, time] = solver.get_results();

    std::cout << std::fixed << std::setprecision(4) << std::setw(15) << t_final
              << std::scientific << std::setprecision(6) << std::setw(20)
              << error << std::fixed << std::setprecision(4) << std::setw(20)
              << time << "s" << std::endl;
  }
}

// Part 3(b): CFL determination (Galerkin)
std::vector<double> determine_cfl_values() {
  std::vector<int> N_values = {16, 32, 48, 64, 96, 128, 192, 256};
  std::vector<double> cfl_values;
  double nu = 0.1;

  std::cout << "\n=== Part 3(b): Determining Maximum CFL Values for "
               "Fourier-Galerkin ==="
            << std::endl;
  std::cout << "Parameters: ν = " << nu << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "Max CFL" << std::endl;
  std::cout << std::string(28, '-') << std::endl;

  for (int N : N_values) {
    auto galerkin = std::make_shared<FourierGalerkin<double>>();
    BurgersGalerkinSolver<double> solver(galerkin, nu);
    solver.initialize(N, 1.0); // Initialize for CFL testing

    double max_cfl = solver.find_cfl_limit();
    cfl_values.push_back(max_cfl);

    std::cout << std::setw(8) << N << std::fixed << std::setprecision(3)
              << std::setw(20) << max_cfl << std::endl;
  }

  std::cout << "\nDoes the definition of the time-step restriction in Eq. 4 "
               "seem reasonable?"
            << std::endl;
  std::cout << "Comparing with collocation method CFL values..." << std::endl;

  return cfl_values;
}

// Part 3(c): Convergence study using CFL values from Part 3(b)
void galerkin_convergence_study(const std::vector<double> &cfl_values) {
  std::vector<int> N_values = {16, 32, 48, 64, 96, 128, 192, 256};
  double t_final = MathConstants<double>::PI() / 4.0;
  double nu = 0.1;

  std::cout
      << "\n=== Part 3(c): Fourier-Galerkin Convergence Study at t = π/4 ==="
      << std::endl;
  std::cout << "Using CFL values from Part 3(b)" << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(15) << "CFL" << std::setw(20)
            << "L∞ Error" << std::setw(20) << "CPU Time" << std::endl;
  std::cout << std::string(63, '-') << std::endl;

  std::vector<double> errors;

  for (size_t i = 0; i < N_values.size(); ++i) {
    int N = N_values[i];
    // double cfl = cfl_values[i];
    double cfl = 0.5;

    auto galerkin = std::make_shared<FourierGalerkin<double>>();
    BurgersGalerkinSolver<double> solver(galerkin, nu);
    solver.initialize(N, t_final, cfl); // Use CFL from Part 3(b)
    solver.solve();

    auto [x, u_numerical, u_exact, error, time, actual_cfl] =
        solver.get_results();
    errors.push_back(error);

    std::cout << std::setw(8) << N << std::fixed << std::setprecision(4)
              << std::setw(15) << cfl << std::scientific << std::setprecision(6)
              << std::setw(20) << error << std::fixed << std::setprecision(4)
              << std::setw(20) << time << "s" << std::endl;
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

  std::cout << "\nWhat is the convergence rate observed?" << std::endl;
  std::cout << "Does it correspond to what you would expect?" << std::endl;
  std::cout << "Fourier-Galerkin should show spectral convergence for smooth "
               "solutions."
            << std::endl;
}
// Part 3(d): Comparison between Galerkin and Collocation
void compare_methods(const std::vector<double> &collocation_cfls,
                     const std::vector<double> &galerkin_cfls) {
  std::vector<int> N_values = {64, 128};
  double t_final = MathConstants<double>::PI() / 4.0;
  double nu = 0.1;

  std::cout << "\n=== Part 3(d): Comparison of Fourier-Galerkin vs "
               "Fourier-Collocation ==="
            << std::endl;
  std::cout << "At t = π/4" << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "Collocation CFL"
            << std::setw(20) << "Galerkin CFL" << std::setw(20)
            << "Collocation Error" << std::setw(20) << "Galerkin Error"
            << std::endl;
  std::cout << std::string(88, '-') << std::endl;

  for (int N : N_values) {
    // Find the appropriate CFL values
    double collocation_cfl = 0.5;
    double galerkin_cfl = 0.5;

    std::vector<int> all_N_values = {16, 32, 48, 64, 96, 128, 192, 256};
    for (size_t i = 0; i < all_N_values.size(); ++i) {
      if (all_N_values[i] == N) {
        collocation_cfl = collocation_cfls[i];
        galerkin_cfl = galerkin_cfls[i];
        break;
      }
    }

    // Collocation method
    auto collocation =
        std::make_shared<SpectralFourier<double>>(MethodType::ODD);
    collocation->build(N);
    BurgersSolver<double> collocation_solver(collocation, nu);
    collocation_solver.initialize(N, t_final, collocation_cfl);
    collocation_solver.solve();
    auto [x_c, u_c, u_exact_c, error_c, time_c] =
        collocation_solver.get_results();

    // Galerkin method
    auto galerkin = std::make_shared<FourierGalerkin<double>>();
    BurgersGalerkinSolver<double> galerkin_solver(galerkin, nu);
    galerkin_solver.initialize(N, t_final);
    galerkin_solver.solve();
    auto [x_g, u_g, u_exact_g, error_g, time_g, actual_cfl_g] =
        galerkin_solver.get_results();

    std::cout << std::setw(8) << N << std::fixed << std::setprecision(4)
              << std::setw(20) << collocation_cfl << std::setw(20)
              << galerkin_cfl << std::scientific << std::setprecision(6)
              << std::setw(20) << error_c << std::setw(20) << error_g
              << std::endl;
  }

  std::cout
      << "\nDo you see any significant difference between the two solutions?"
      << std::endl;
  std::cout << "Any reason to use the Galerkin formulation?" << std::endl;
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
  bool run_ex02c = run_all || contains_arg(argc, argv, "--ex02c");
  bool run_ex02d = run_all || contains_arg(argc, argv, "--ex02d");
  bool run_ex03b = run_all || contains_arg(argc, argv, "--ex03b");
  bool run_ex03c = run_all || contains_arg(argc, argv, "--ex03c");
  bool run_ex03d = run_all || contains_arg(argc, argv, "--ex03d");

  if (!run_ex02b && !run_ex02c && !run_ex02d && !run_ex03b && !run_ex03c &&
      !run_ex03d) {
    std::cerr << "Error: No valid exercise options specified." << std::endl;
    print_usage(argv[0]);
    return 1;
  }

  // Store CFL values between parts
  std::vector<double> collocation_cfls;
  std::vector<double> galerkin_cfls;

  // Part 2: Fourier Collocation
  if (run_ex02b || run_ex02c || run_ex02d) {
    std::cout << "\n=== Exercise 3 Part 2: Burgers' Equation (Fourier "
                 "Collocation) ==="
              << std::endl;
    std::cout << "Using odd number of grid points: xⱼ = 2πj/(N+1), j ∈ [0,N]"
              << std::endl;

    if (run_ex02b) {
      collocation_cfls = burgers_find_max_cfl();
    }

    if (run_ex02c) {
      if (collocation_cfls.empty()) {
        std::cout << "\nWarning: Part 2(c) requires CFL values from Part 2(b)."
                  << std::endl;
        std::cout << "Running Part 2(b) first..." << std::endl;
        collocation_cfls = burgers_find_max_cfl();
      }
      burgers_convergence_study_with_cfl(collocation_cfls);
    }

    if (run_ex02d) {
      if (collocation_cfls.empty()) {
        std::cout << "\nWarning: Part 2(d) requires CFL values from Part 2(b)."
                  << std::endl;
        std::cout << "Running Part 2(b) first..." << std::endl;
        collocation_cfls = burgers_find_max_cfl();
      }
      burgers_time_evolution(collocation_cfls);
    }
  }

  // Part 3: Fourier Galerkin
  if (run_ex03b || run_ex03c || run_ex03d) {
    std::cout
        << "\n=== Exercise 3 Part 3: Burgers' Equation (Fourier Galerkin) ==="
        << std::endl;

    if (run_ex03b) {
      galerkin_cfls = determine_cfl_values();
    }

    if (run_ex03c) {
      if (galerkin_cfls.empty()) {
        std::cout << "\nWarning: Part 3(c) requires CFL values from Part 3(b)."
                  << std::endl;
        std::cout << "Running Part 3(b) first..." << std::endl;
        galerkin_cfls = determine_cfl_values();
      }
      galerkin_convergence_study(galerkin_cfls);
    }

    if (run_ex03d) {
      if (collocation_cfls.empty() || galerkin_cfls.empty()) {
        std::cout << "\nWarning: Part 3(d) requires CFL values from Parts 2(b) "
                     "and 3(b)."
                  << std::endl;
        if (collocation_cfls.empty()) {
          std::cout << "Running Part 2(b) first..." << std::endl;
          collocation_cfls = burgers_find_max_cfl();
        }
        if (galerkin_cfls.empty()) {
          std::cout << "Running Part 3(b) first..." << std::endl;
          galerkin_cfls = determine_cfl_values();
        }
      }
      compare_methods(collocation_cfls, galerkin_cfls);
    }
  }

  return 0;
}
