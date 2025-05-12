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
      u_exact_[i] =
          TestFunctions::burgers_exact_u<T>(differentiator_->get_x()[i], t_final_, T(4), nu_);
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
#endif // INCLUDE_INCLUDE_BURGERS_H_
