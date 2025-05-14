#ifndef INCLUDE_INCLUDE_BURGERS_H_
#define INCLUDE_INCLUDE_BURGERS_H_

#include "common.h"
#include "diff.h"
#include "finite_diff.h"
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

  void initialize(int N, T t_final) {
    N_ = N;
    t_final_ = t_final;

    // Create grid points
    differentiator_->create_grid_pts(N);

    // Calculate time step (CFL condition for Burgers' equation)
    T dx = T(2) * MathConstants<T>::PI() / (N + 1);
    T c_max = T(4); // Maximum expected velocity

    // CFL condition: dt <= min(dx/c_max, dx²/(2ν))
    T dt_adv = T(0.1) * dx / c_max;
    T dt_diff = T(0.1) * dx * dx / (T(2) * nu_);
    dt_ = std::min(dt_adv, dt_diff);

    // Additional safety for different schemes
    if (dynamic_cast<SecondOrderFiniteDiff<T> *>(differentiator_.get())) {
      dt_ *= T(0.5);
    } else if (dynamic_cast<FourthOrderFiniteDiff<T> *>(
                   differentiator_.get())) {
      dt_ *= T(0.25);
    } else if (dynamic_cast<SpectralFourier<T> *>(differentiator_.get())) {
      dt_ *= T(0.1); // Be more conservative with spectral method
    }

    // Adjust dt to ensure we hit t_final exactly
    num_steps_ = static_cast<int>(ceil(t_final / dt_));
    dt_ = t_final / T(num_steps_);

    // Initialize solution vector with initial condition
    u_ = std::vector<T>(N + 1);
    for (int i = 0; i <= N; i++) {
      u_[i] = TestFunctions::burgers_initial_u<T>(differentiator_->get_x()[i]);
    }
  }

  void solve() {
    auto start_time = std::chrono::high_resolution_clock::now();

    T t = T(0);
    for (int step = 0; step < num_steps_; ++step) {
      // Use the specialized Burgers' equation integrator
      u_ = integrator_->integrate_burgers(u_, dt_, *differentiator_, nu_);
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

private:
  std::shared_ptr<Differentiator<T>> differentiator_;
  std::shared_ptr<RungeKutta4<T>> integrator_;
  std::vector<T> u_;
  std::vector<T> u_exact_;
  T dt_;
  T t_final_;
  T nu_;
  int N_;
  int num_steps_;
  double computation_time_;
};

template <NumericType T> class BurgersGalerkinSolver {
public:
  BurgersGalerkinSolver(std::shared_ptr<FourierGalerkin<T>> galerkin,
                        T nu = T(0.1))
      : galerkin_(galerkin), nu_(nu) {}

  void initialize(int N, T t_final) {
    N_ = N;
    t_final_ = t_final;

    // Create grid points
    galerkin_->create_grid_pts(N);

    // Calculate CFL for time step
    T dx = T(2) * MathConstants<T>::PI() / (N + 1);
    T k_max = T(N) / T(2);

    // Time step based on exam formula (Eq. 4)
    // dt ≤ CFL × [max|u(x_j)|*k_max + ν*(k_max)²]^(-1)
    T u_max = T(4); // Estimated maximum velocity for Burgers' sawtooth
    T dt_constraint = u_max * k_max + nu_ * k_max * k_max;

    // We'll determine CFL experimentally, but start conservative
    cfl_ = T(0.5);
    dt_ = cfl_ / dt_constraint;

    // Adjust dt to ensure we hit t_final exactly
    num_steps_ = static_cast<int>(ceil(t_final / dt_));
    dt_ = t_final / T(num_steps_);

    // Initialize solution vector with initial condition
    u_ = std::vector<T>(N + 1);
    for (int i = 0; i <= N; i++) {
      u_[i] = TestFunctions::burgers_initial_u<T>(galerkin_->get_x()[i]);
    }
  }

  void solve() {
    auto start_time = std::chrono::high_resolution_clock::now();

    T t = T(0);
    for (int step = 0; step < num_steps_; ++step) {
      // Use 4th order RK as specified in the exam
      u_ = runge_kutta_4(u_, dt_);
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

  // Implementation of the specific RK4 scheme from the exam
  std::vector<T> runge_kutta_4(const std::vector<T> &u, T dt) {
    size_t N = u.size() - 1;
    std::vector<T> u1(N + 1), u2(N + 1), u3(N + 1), un_p1(N + 1);

    // Convert to spectral space
    auto u_hat = galerkin_->fft(u);

    // Step 1: u1 = u + dt/2 * F(u)
    auto k1_hat = galerkin_->compute_burgers_spectral_rhs(u, nu_);
    std::vector<std::complex<T>> u1_hat(N + 1);
    for (size_t i = 0; i <= N; ++i) {
      u1_hat[i] = u_hat[i] + dt / T(2) * k1_hat[i];
    }
    u1 = galerkin_->ifft(u1_hat);

    // Step 2: u2 = u + dt/2 * F(u1)
    auto k2_hat = galerkin_->compute_burgers_spectral_rhs(u1, nu_);
    std::vector<std::complex<T>> u2_hat(N + 1);
    for (size_t i = 0; i <= N; ++i) {
      u2_hat[i] = u_hat[i] + dt / T(2) * k2_hat[i];
    }
    u2 = galerkin_->ifft(u2_hat);

    // Step 3: u3 = u + dt * F(u2)
    auto k3_hat = galerkin_->compute_burgers_spectral_rhs(u2, nu_);
    std::vector<std::complex<T>> u3_hat(N + 1);
    for (size_t i = 0; i <= N; ++i) {
      u3_hat[i] = u_hat[i] + dt * k3_hat[i];
    }
    u3 = galerkin_->ifft(u3_hat);

    // Step 4: un+1 = (1/3)*(-u + u1 + 2*u2 + u3 + dt/2 * F(u3))
    auto k4_hat = galerkin_->compute_burgers_spectral_rhs(u3, nu_);
    std::vector<std::complex<T>> un_p1_hat(N + 1);
    for (size_t i = 0; i <= N; ++i) {
      un_p1_hat[i] = (T(1) / T(3)) * (-u_hat[i] + u1_hat[i] + T(2) * u2_hat[i] +
                                      u3_hat[i] + dt / T(2) * k4_hat[i]);
    }
    un_p1 = galerkin_->ifft(un_p1_hat);

    return un_p1;
  }

  // Find experimental CFL for stability
  T find_cfl_limit() {
    T test_cfl = T(2.0);
    T cfl_step = T(0.1);
    bool stable = true;

    while (test_cfl > T(0.1) && stable) {
      // Try current CFL value
      cfl_ = test_cfl;
      initialize(N_, T(0.1)); // Short test run

      // Run a few steps and check for stability
      stable = test_stability();

      if (!stable) {
        test_cfl -= cfl_step;
      } else {
        // If stable, this is our CFL limit
        return test_cfl;
      }
    }

    return T(0.1); // Conservative fallback
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
        u_test = runge_kutta_4(u_test, dt_);

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
