#ifndef INCLUDE_INCLUDE_FINITE_DIFF_H_
#define INCLUDE_INCLUDE_FINITE_DIFF_H_

#include "diff.h"

template <NumericType T>
class SecondOrderFiniteDiff : public Differentiator<T> {
public:
  SecondOrderFiniteDiff() = default;
  void compute_num_sol() override;
  void create_grid_pts(int N) override;

private:
  T dx_;
};

template <NumericType T> class FourthOrderFiniteDiff : public Differentiator<T> {
public:
  FourthOrderFiniteDiff() = default;
  void compute_num_sol() override;
  void create_grid_pts(int N) override;

private:
  T dx_;
};

template <NumericType T> void SecondOrderFiniteDiff<T>::compute_num_sol() {
  size_t N = this->analytical_u_.size() - 1;
  this->numerical_du_ = std::vector<T>(this->analytical_u_.size(), 0.0);

  for (size_t j = 1; j < N; ++j) {
    this->numerical_du_[j] =
        (this->analytical_u_[j + 1] - this->analytical_u_[j - 1]) / (2.0 * dx_);
  }

  // Periodic boundary conditions
  this->numerical_du_[0] =
      (this->analytical_u_[1] - this->analytical_u_[N]) / (2.0 * dx_);
  this->numerical_du_[N] =
      (this->analytical_u_[0] - this->analytical_u_[N - 1]) / (2.0 * dx_);
}

template <NumericType T> void FourthOrderFiniteDiff<T>::compute_num_sol() {
  size_t N = this->analytical_u_.size() - 1;
  this->numerical_du_ = std::vector<T>(this->analytical_u_.size(), 0.0);

  // Interior points
  for (size_t j = 2; j < N - 1; ++j) {
    this->numerical_du_[j] =
        (-this->analytical_u_[j + 2] + 8.0 * this->analytical_u_[j + 1] -
         8.0 * this->analytical_u_[j - 1] + this->analytical_u_[j - 2]) /
        (12.0 * dx_);
  }

  // Handle periodic boundary conditions correctly
  // For j = 0, we need indices -2, -1, 1, 2 which map to N-2, N-1, 1, 2
  this->numerical_du_[0] =
      (-this->analytical_u_[2] + 8.0 * this->analytical_u_[1] -
       8.0 * this->analytical_u_[N] + this->analytical_u_[N - 1]) /
      (12.0 * dx_);

  // For j = 1, we need indices -1, 0, 2, 3 which map to N-1, 0, 2, 3
  this->numerical_du_[1] =
      (-this->analytical_u_[3] + 8.0 * this->analytical_u_[2] -
       8.0 * this->analytical_u_[0] + this->analytical_u_[N]) /
      (12.0 * dx_);

  // For j = N-1, we need indices N-3, N-2, N, N+1 which map to N-3, N-2, 0, 1
  this->numerical_du_[N - 1] =
      (-this->analytical_u_[0] + 8.0 * this->analytical_u_[N] -
       8.0 * this->analytical_u_[N - 2] + this->analytical_u_[N - 3]) /
      (12.0 * dx_);

  // For j = N, we need indices N-2, N-1, N+1, N+2 which map to N-2, N-1, 1, 2
  this->numerical_du_[N] =
      (-this->analytical_u_[1] + 8.0 * this->analytical_u_[0] -
       8.0 * this->analytical_u_[N - 1] + this->analytical_u_[N - 2]) /
      (12.0 * dx_);
}

// TODO: (dhub) Move this to Differentiator class?
template <NumericType T> void SecondOrderFiniteDiff<T>::create_grid_pts(int N) {
  int grid_pts = N + 1;
  this->x_ = std::vector<T>(grid_pts);
  this->dx_ = T(2) * MathConstants<T>::PI() / (grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    this->x_[i] = i * this->dx_;
  }
}
template <NumericType T> void FourthOrderFiniteDiff<T>::create_grid_pts(int N) {
  int grid_pts = N + 1;
  this->x_ = std::vector<T>(grid_pts);
  this->dx_ = T(2) * MathConstants<T>::PI() / (grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    this->x_[i] = i * this->dx_;
  }
}

#endif // INCLUDE_INCLUDE_FINITE_DIFF_H_
