#ifndef INCLUDE_INCLUDE_BURGERS_H_
#define INCLUDE_INCLUDE_BURGERS_H_

#include "common.h"
#include "diff.h"
#include "fourier.h"
#include "integrate.h"
#include <chrono>
#include <memory>
#include <vector>

template <NumericType T> class BurgersSolver {
public:
  BurgersSolver(std::shared_ptr<Differentiator<T>> differentiator,
                T nu = T(0.1))
      : differentiator_(differentiator), nu_(nu),
        integrator_(std::make_shared<RungeKutta4<T>>()) {}

  void initialize(int N, T t_final, T cfl) {
    N_ = N;
    t_final_ = t_final;
    cfl_ = cfl;

    // Create grid points (ODD method as specified)
    differentiator_->create_grid_pts(N);

    // Initialize solution vector with initial condition
    u_ = std::vector<T>(N + 1);
    for (int i = 0; i <= N; i++) {
      u_[i] = TestFunctions::burgers_initial_u<T>(differentiator_->get_x()[i]);
    }
  }

  void solve() {
    auto start_time = std::chrono::high_resolution_clock::now();

    T t = T(0);
    T dx = T(2) * MathConstants<T>::PI() / T(N_ + 1);
    T dt_constraint_sec = nu_ / (dx * dx);
    for (int step = 0; step < num_steps_; ++step) {
      // Use the specialized Burgers' equation integrator
      T max_curr = T(0);

      for (int j = 0; j <= N_; j++) {
        T max_pot = abs(u_[j]) / dx + dt_constraint_sec;
        if (max_pot > max_curr) {
          max_curr = max_pot;
        }
      }

      T dt = cfl_ / max_curr;

      u_ = integrator_->integrate_burgers(u_, dt, *differentiator_, nu_);
      t += dt_;

      // Apply periodic boundary conditions
      u_[0] = u_[N_];
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    computation_time_ = elapsed.count();

    // Compute exact solution for comparison
    u_exact_ = std::vector<T>(N_ + 1);
    for (int i = 0; i <= N_; i++) {
      u_exact_[i] = TestFunctions::burgers_exact_u<T>(
          differentiator_->get_x()[i], t_final_, T(4), nu_);
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

  std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, T, double>
  get_results() const {
    return {differentiator_->get_x(), u_, u_exact_, compute_error(),
            computation_time_};
  }
  T find_max_cfl(int N) {
    N_ = N;
    differentiator_->create_grid_pts(N);

    // Binary search for maximum stable CFL
    T cfl_low = T(0.01);
    T cfl_high = T(4.0);
    T cfl_tolerance = T(0.01);

    while (cfl_high - cfl_low > cfl_tolerance) {
      T cfl_test = (cfl_low + cfl_high) / T(2);

      // Test stability with a short run
      if (is_stable_with_cfl(cfl_test, 25)) {
        cfl_low = cfl_test;
      } else {
        cfl_high = cfl_test;
      }
    }

    return cfl_low; // Return maximum stable CFL
  }

  bool is_stable_with_cfl(T cfl, int test_steps) {
    T dx = T(2) * MathConstants<T>::PI() / T(N_ + 1);

    std::vector<T> u_test(N_ + 1);
    for (int i = 0; i <= N_; i++) {
      u_test[i] =
          TestFunctions::burgers_initial_u<T>(differentiator_->get_x()[i]);
    }

    T dt_constraint_sec = nu_ / (dx * dx);

    for (int step = 0; step < test_steps; ++step) {
      T max_curr = T(0);

      for (int j = 0; j <= N_; j++) {
        T max_pot = abs(u_test[j]) / dx + dt_constraint_sec;
        if (max_pot > max_curr) {
          max_curr = max_pot;
        }
      }

      T dt = cfl / max_curr;

      // Take one time step
      u_test =
          integrator_->integrate_burgers(u_test, dt, *differentiator_, nu_);

      // Check for instability (NaN or exponential growth)
      // TODO: (dhub) Check for potential change in test_u prev and test_U curr
      for (const auto &val : u_test) {
        if (std::isnan(val) || std::isinf(val) || abs(val) > T(1e10)) {
          return false;
        }
      }
    }
    return true;
  }

private:
  std::shared_ptr<Differentiator<T>> differentiator_;
  std::shared_ptr<RungeKutta4<T>> integrator_;
  std::vector<T> u_;
  std::vector<T> u_exact_;
  T dt_;
  T t_final_;
  T nu_;
  T cfl_;
  int N_;
  int num_steps_;
  double computation_time_;
};

template <NumericType T> class BurgersGalerkinSolver {
public:
  BurgersGalerkinSolver(std::shared_ptr<FourierGalerkin<T>> galerkin,
                        T nu = T(0.1))
      : galerkin_(galerkin), nu_(nu),
        integrator_(std::make_shared<RungeKutta4<T>>()) {}

  void initialize(int N, T t_final, T cfl = T(0.5)) {
    N_ = N;
    t_final_ = t_final;
    cfl_ = cfl;

    // Create grid points
    galerkin_->create_grid_pts(N);

    // Initialize solution with quadrature formula
    u_ = std::vector<T>(N + 1);
    for (int i = 0; i <= N; i++) {
      u_[i] = TestFunctions::burgers_initial_u<T>(galerkin_->get_x()[i]);
    }

    // Calculate time step using Eq. (4) from exam
    T k_max = T(N) / T(2);

    // Find max|u(xⱼ)| from current solution
    T max_u = T(0);
    for (int j = 0; j <= N; j++) {
      T abs_u = abs(u_[j]);
      if (abs_u > max_u) {
        max_u = abs_u;
      }
    }

    // Apply formula (4): dt ≤ CFL × [max|u(xⱼ)|*k_max + ν*(k_max)²]⁻¹
    T dt_constraint = max_u * k_max + nu_ * k_max * k_max;
    dt_ = cfl_ / dt_constraint;

    // Adjust dt to ensure we hit t_final exactly
    num_steps_ = static_cast<int>(ceil(t_final / dt_));
    dt_ = t_final / T(num_steps_);
  }

  void solve() {
    auto start_time = std::chrono::high_resolution_clock::now();

    T t = T(0);
    for (int step = 0; step < num_steps_; ++step) {
      // Use the existing RungeKutta4 integrator
      u_ = integrator_->integrate_burgers(u_, dt_, *galerkin_, nu_);
      t += dt_;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    computation_time_ = elapsed.count();

    // Compute exact solution for comparison
    u_exact_ = std::vector<T>(N_ + 1);
    for (int i = 0; i <= N_; i++) {
      u_exact_[i] = TestFunctions::burgers_exact_u<T>(galerkin_->get_x()[i],
                                                      t_final_, T(4), nu_);
    }
  }

  // Find experimental CFL for stability
  T find_cfl_limit() {
    N_ = galerkin_->get_x().size() - 1;

    // Initialize solution
    u_ = std::vector<T>(N_ + 1);
    for (int i = 0; i <= N_; i++) {
      u_[i] = TestFunctions::burgers_initial_u<T>(galerkin_->get_x()[i]);
    }

    // Calculate k_max
    T k_max = T(N_) / T(2);

    // Find max|u(xⱼ)|
    T max_u = T(0);
    for (int j = 0; j <= N_; j++) {
      T abs_u = abs(u_[j]);
      if (abs_u > max_u) {
        max_u = abs_u;
      }
    }

    // Calculate the constraint denominator using Eq. (4)
    T dt_constraint = max_u * k_max + nu_ * k_max * k_max;

    // Binary search for maximum stable CFL
    T cfl_low = T(0.01);
    T cfl_high = T(1.0); // Galerkin is typically more stable
    T cfl_tolerance = T(0.01);

    while (cfl_high - cfl_low > cfl_tolerance) {
      T cfl_test = (cfl_low + cfl_high) / T(2);

      // Apply formula (4) with this CFL
      dt_ = cfl_test / dt_constraint;
      cfl_ = cfl_test;

      // Test stability with a short run
      if (test_stability()) {
        cfl_low = cfl_test;
      } else {
        cfl_high = cfl_test;
      }
    }

    return cfl_low;
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

  std::tuple<std::vector<T>, std::vector<T>, std::vector<T>, T, double, T>
  get_results() const {
    return {galerkin_->get_x(), u_,  u_exact_, compute_error(),
            computation_time_,  cfl_};
  }

private:
  bool test_stability() {
    try {
      T t = T(0);
      int test_steps = 10;
      auto u_test = u_;

      for (int step = 0; step < test_steps; ++step) {
        // Use the existing integrator
        u_test = integrator_->integrate_burgers(u_test, dt_, *galerkin_, nu_);

        // Check for NaN or very large values
        for (const auto &val : u_test) {
          if (std::isnan(val) || std::isinf(val) || abs(val) > T(1e10)) {
            return false;
          }
        }
        t += dt_;
      }

      return true;
    } catch (...) {
      return false;
    }
  }

  std::shared_ptr<FourierGalerkin<T>> galerkin_;
  std::shared_ptr<RungeKutta4<T>> integrator_; // Use the existing integrator
  std::vector<T> u_;
  std::vector<T> u_exact_;
  T dt_;
  T t_final_;
  T nu_;
  T cfl_;
  int N_;
  int num_steps_;
  double computation_time_;
};

#endif // INCLUDE_INCLUDE_BURGERS_H_
