#ifndef INCLUDE_INCLUDE_FINITE_DIFF_H_
#define INCLUDE_INCLUDE_FINITE_DIFF_H_

#include "diff.h"

template <NumericType T>
class SecondOrderFiniteDiff : public Differentiator<T> {
public:
  SecondOrderFiniteDiff() = default;  // Removed <T> from constructor
  void compute_num_sol() override;
  void create_grid_pts(int N) override;

private:
  T dx_;
};

template <NumericType T> class ForthOrderFiniteDiff : public Differentiator<T> {
public:
  ForthOrderFiniteDiff() = default;  // Removed <T> from constructor
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

template <NumericType T> void ForthOrderFiniteDiff<T>::compute_num_sol() {
  size_t N = this->analytical_u_.size() - 1;
  this->numerical_du_ = std::vector<T>(this->analytical_u_.size(), 0.0);

  // Interior points
  for (size_t j = 2; j < N - 1; ++j) {
    this->numerical_du_[j] =
        (-this->analytical_u_[j + 2] + 8.0 * this->analytical_u_[j + 1] -
         8.0 * this->analytical_u_[j - 1] + this->analytical_u_[j - 2]) /
        (12.0 * dx_);
  }

  // Handle periodic boundary conditions
  this->numerical_du_[0] =
      (-this->analytical_u_[2] + 8.0 * this->analytical_u_[1] -
       8.0 * this->analytical_u_[N - 1] + this->analytical_u_[N - 2]) /
      (12.0 * dx_);
  this->numerical_du_[1] =
      (-this->analytical_u_[3] + 8.0 * this->analytical_u_[2] -
       8.0 * this->analytical_u_[0] + this->analytical_u_[N]) /
      (12.0 * dx_);
  this->numerical_du_[N - 1] =
      (-this->analytical_u_[1] + 8.0 * this->analytical_u_[0] -
       8.0 * this->analytical_u_[N - 2] + this->analytical_u_[N - 3]) /
      (12.0 * dx_);
  this->numerical_du_[N] =
      (-this->analytical_u_[0] + 8.0 * this->analytical_u_[N] -
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
template <NumericType T> void ForthOrderFiniteDiff<T>::create_grid_pts(int N) {
  int grid_pts = N + 1;
  this->x_ = std::vector<T>(grid_pts);
  this->dx_ = T(2) * MathConstants<T>::PI() / (grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    this->x_[i] = i * this->dx_;
  }
}

#endif // INCLUDE_INCLUDE_FINITE_DIFF_H_
