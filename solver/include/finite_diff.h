#ifndef INCLUDE_INCLUDE_FINITE_DIFF_H_
#define INCLUDE_INCLUDE_FINITE_DIFF_H_

#include "diff.h"

template <NumericType T> class SecondOrderFiniteDiff : public Differentiator<T> {
public:
  void compute_num_sol() override;

private:
  double dx_;
};

template <NumericType T> class ForthOrderFiniteDiff : public Differentiator<T> {
public:
  void compute_num_sol() override;

private:
  T dx_;
};

template <NumericType T> void SecondOrderFiniteDiff<T>::compute_num_sol() {
  size_t N = this->analytical_u_.size() - 1;
  std::vector<T> dudt(this->analytical_u_.size(), 0.0);

  // Interior points
  for (size_t j = 2; j < N - 1; ++j) {
    dudt[j] = (-this->analytical_u_[j + 2] + 8.0 * this->analytical_u_[j + 1] -
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

#endif // INCLUDE_INCLUDE_FINITE_DIFF_H_
