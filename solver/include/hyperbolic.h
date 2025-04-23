#ifndef INCLUDE_INCLUDE_HYPERBOLIC_H_
#define INCLUDE_INCLUDE_HYPERBOLIC_H_

#include "common.h"
#include "diff.h"
#include "finite_diff.h"
#include "fourier.h"
#include "integrate.h"
#include "matplot/freestanding/plot.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <matplot/matplot.h>
#include <memory>
#include <string>

// Using functions from TestFunctions namespace
using namespace TestFunctions;

template <NumericType T> class HyperbolicSolver {
public:
  HyperbolicSolver(std::shared_ptr<Differentiator<T>> differentiator)
      : differentiator_(differentiator),
        integrator_(std::make_shared<RungeKutta4<T>>()) {}

  void initialize(int N, T t_final) {
    N_ = N;
    t_final_ = t_final;

    // Create grid points
    differentiator_->create_grid_pts(N);

    // Set initial condition
    differentiator_->compute_analy_sols(TestFunctions::hyperbolic_initial_u<T>,
                                        TestFunctions::hyperbolic_du_t0<T>);

    // Calculate time step (CFL condition)
    T dx = T(2) * MathConstants<T>::PI() / (N + 1);

    // Set different dt based on the method type
    if (dynamic_cast<SecondOrderFiniteDiff<T> *>(differentiator_.get())) {
      dt_ = T(0.5) * dx / (T(2) * MathConstants<T>::PI());
    } else if (dynamic_cast<FourthOrderFiniteDiff<T> *>(
                   differentiator_.get())) {
      dt_ = T(0.25) * dx / (T(2) * MathConstants<T>::PI());
    } else {
      // Fourier method - more restrictive for stability
      dt_ = T(0.1) * dx / (T(2) * MathConstants<T>::PI());
    }

    // Adjust dt to ensure we hit t_final exactly
    num_steps_ = static_cast<int>(ceil(t_final / dt_));
    dt_ = t_final / T(num_steps_);

    // Initialize solution vector with initial condition
    u_ = std::vector<T>(N + 1);
    for (int i = 0; i <= N; i++) {
      // Make sure to pass 0 as the second argument to avoid ambiguity
      u_[i] = TestFunctions::hyperbolic_initial_u<T>(
          differentiator_->get_x()[i], 0);
    }
  }

  void solve() {
    auto start_time = std::chrono::high_resolution_clock::now();

    T t = T(0);
    for (int step = 0; step < num_steps_; ++step) {
      // Perform time integration step
      // std::cout << "Timestep: " << step << std::endl;
      u_ = integrator_->integrate(
          u_, dt_, *differentiator_); // This will now use the non-const version

      t += dt_;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    computation_time_ = elapsed.count();

    // Compute exact solution for comparison
    u_exact_ = std::vector<T>(N_ + 1);
    for (int i = 0; i <= N_; i++) {
      u_exact_[i] = TestFunctions::hyperbolic_exact_u<T>(
          differentiator_->get_x()[i], t_final_);
    }
  }

  T compute_error() const {
    T max_error = T(0);
    for (int i = 0; i <= N_; i++) {
      T error = abs(u_[i] - u_exact_[i]);
      if (error > max_error) {
        max_error = error;
      }
    }
    return max_error;
  }

  void save_results(const std::string &filename) const {
    std::ofstream file(filename);
    file << "x,u_computed,u_exact,error\n";

    for (int i = 0; i <= N_; i++) {
      file << differentiator_->get_x()[i] << "," << u_[i] << "," << u_exact_[i]
           << "," << abs(u_[i] - u_exact_[i]) << "\n";
    }

    file.close();
  }

  std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, T, double>
  get_results() const {
    return {differentiator_->get_x(), u_, u_exact_, compute_error(),
            computation_time_};
  }

private:
  std::shared_ptr<Differentiator<T>> differentiator_;
  std::shared_ptr<Integrator<T>> integrator_;
  std::vector<T> u_;
  std::vector<T> u_exact_;
  T dt_;
  T t_final_;
  int N_;
  int num_steps_;
  double computation_time_;
};

template <NumericType T>
void plot_func(std::vector<T> &x, std::vector<T> &u_exact,
               std::vector<T> &u_numeric, std::string type, int t) {
  using namespace matplot;

  auto fig = figure(true);
  fig->quiet_mode(true);
  fig->backend()->run_command("unset warnings");
  fig->size(1280, 960);
  fig->font_size(18);

  auto p2 = plot(x, u_numeric, "b-o");
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
  ax->title("Solution for Hyperbolic Problem at t=" + std::to_string(t) +
            " using " + type);

  legend()->font_size(16);

  save("hyperbolic_" + type + "_" + std::to_string(t) + ".png");
}

// Function to run convergence study
template <NumericType T> void convergence_study() {
  std::vector<int> N_values = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
  T t_final = MathConstants<T>::PI();

  std::vector<T> errors_second(N_values.size());
  std::vector<T> errors_fourth(N_values.size());
  std::vector<T> errors_fourier(N_values.size());

  std::vector<double> times_second(N_values.size());
  std::vector<double> times_fourth(N_values.size());
  std::vector<double> times_fourier(N_values.size());

  std::cout << "Performing convergence study at t = π..." << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "Second Order"
            << std::setw(20) << "Fourth Order" << std::setw(20) << "Fourier"
            << std::endl;
  std::cout << std::string(68, '-') << std::endl;

  for (size_t i = 0; i < N_values.size(); ++i) {
    int N = N_values[i];
    std::cout << std::setw(8) << N;

    // Second order method
    auto second_order = std::make_shared<SecondOrderFiniteDiff<T>>();
    HyperbolicSolver<T> solver_second(second_order);
    solver_second.initialize(N, t_final);
    solver_second.solve();
    auto [_, __, ___, error_second, time_second] = solver_second.get_results();
    errors_second[i] = error_second;
    times_second[i] = time_second;
    std::cout << std::scientific << std::setprecision(6) << std::setw(20)
              << error_second;

    // Fourth order method
    auto fourth_order = std::make_shared<FourthOrderFiniteDiff<T>>();
    HyperbolicSolver<T> solver_fourth(fourth_order);
    solver_fourth.initialize(N, t_final);
    solver_fourth.solve();
    auto [____, _____, ______, error_fourth, time_fourth] =
        solver_fourth.get_results();
    errors_fourth[i] = error_fourth;
    times_fourth[i] = time_fourth;
    std::cout << std::setw(20) << error_fourth;

    // Fourier method
    auto fourier = std::make_shared<SpectralFourier<T>>(MethodType::ODD);
    fourier->build(N); // Build Fourier differentiation matrix
    HyperbolicSolver<T> solver_fourier(fourier);
    solver_fourier.initialize(N, t_final);
    solver_fourier.solve();
    auto [_______, ________, _________, error_fourier, time_fourier] =
        solver_fourier.get_results();
    errors_fourier[i] = error_fourier;
    times_fourier[i] = time_fourier;
    std::cout << std::setw(20) << error_fourier << std::endl;
  }

  // Calculate convergence rates
  std::cout << "\nConvergence Rates:" << std::endl;
  std::cout << std::setw(8) << "N" << std::setw(20) << "Second Order"
            << std::setw(20) << "Fourth Order" << std::setw(20) << "Fourier"
            << std::endl;
  std::cout << std::string(68, '-') << std::endl;

  for (size_t i = 1; i < N_values.size(); ++i) {
    T rate_second = log(errors_second[i - 1] / errors_second[i]) / log(T(2));
    T rate_fourth = log(errors_fourth[i - 1] / errors_fourth[i]) / log(T(2));
    T rate_fourier = log(errors_fourier[i - 1] / errors_fourier[i]) / log(T(2));

    std::cout << std::setw(8) << N_values[i] << std::fixed
              << std::setprecision(2) << std::setw(20) << rate_second
              << std::setw(20) << rate_fourth << std::setw(20) << rate_fourier
              << std::endl;
  }

  // Compare N values needed for same accuracy
  T error_second_2048 = errors_second.back();
  std::cout << "\nError for N=2048 using second order: " << std::scientific
            << error_second_2048 << std::endl;

  // Find equivalent N for fourth order
  std::cout << "Equivalent N needed to achieve the same error:" << std::endl;

  for (size_t i = 0; i < errors_fourth.size(); ++i) {
    if (errors_fourth[i] <= error_second_2048) {
      T N_eq_fourth;
      if (i > 0) {
        // Interpolate between values
        N_eq_fourth =
            N_values[i - 1] + (N_values[i] - N_values[i - 1]) *
                                  (error_second_2048 - errors_fourth[i - 1]) /
                                  (errors_fourth[i] - errors_fourth[i - 1]);
      } else {
        N_eq_fourth = N_values[i];
      }
      std::cout << "  Fourth Order: N approx " << N_eq_fourth << std::endl;
      break;
    }
  }

  // Find equivalent N for Fourier
  for (size_t i = 0; i < errors_fourier.size(); ++i) {
    if (errors_fourier[i] <= error_second_2048) {
      T N_eq_fourier;
      if (i > 0) {
        // Interpolate between values
        N_eq_fourier =
            N_values[i - 1] + (N_values[i] - N_values[i - 1]) *
                                  (error_second_2048 - errors_fourier[i - 1]) /
                                  (errors_fourier[i] - errors_fourier[i - 1]);
      } else {
        N_eq_fourier = N_values[i];
      }
      std::cout << "  Fourier: N approx " << N_eq_fourier << std::endl;
      break;
    }
  }
}

// Function to perform long time integration comparison
template <NumericType T> void long_time_integration() {
  std::vector<T> times = {T(0), T(100), T(200)};

  // Second order with N=200 as specified in the exercise
  int N_second = 200;

  // Fourier with N=10 as specified in the exercise
  int N_fourier = 10;

  std::cout << "\nPerforming long time integration comparison..." << std::endl;

  for (T t_final : times) {
    std::cout << "Time t = " << t_final << std::endl;

    // Second order
    auto second_order = std::make_shared<SecondOrderFiniteDiff<T>>();
    HyperbolicSolver<T> solver_second(second_order);
    solver_second.initialize(N_second, t_final);
    solver_second.solve();
    auto [x_second, u_second, u_exact_second, error_second, time_second] =
        solver_second.get_results();

    plot_func(x_second, u_exact_second, u_second, "SecondOrderFiniteDiff",
              int(t_final));

    // Fourier
    auto fourier = std::make_shared<SpectralFourier<T>>(MethodType::ODD);
    fourier->build(N_fourier);
    HyperbolicSolver<T> solver_fourier(fourier);
    solver_fourier.initialize(N_fourier, t_final);
    solver_fourier.solve();
    auto [x_fourier, u_fourier, u_exact_fourier, error_fourier, time_fourier] =
        solver_fourier.get_results();

    plot_func(x_fourier, u_exact_fourier, u_fourier, "SpectralFourier",
              int(t_final));

    std::cout << "  Second order (N=" << N_second << ") error: " << error_second
              << std::endl;
    std::cout << "  Fourier (N=" << N_fourier << ") error: " << error_fourier
              << std::endl;
  }
}

#endif // INCLUDE_INCLUDE_HYPERBOLIC_H_
