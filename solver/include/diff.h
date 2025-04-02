#ifndef INCLUDE_INCLUDE_DIFF_H_
#define INCLUDE_INCLUDE_DIFF_H_

#include "common.h"
#include "error.h"

template <NumericType T> class Differentiator {
public:
  Differentiator() = default;
  virtual ~Differentiator() = default;
  virtual void compute_num_sol() = 0;
  virtual void compute_analy_sols(AnalyticalFunction<T> ana_u,
                                  AnalyticalFunction<T> ana_du, int k = 0);
  virtual void create_grid_pts(int N) = 0;
  virtual T L_inf_error();
  virtual T L_2_error();

protected:
  std::vector<T> x_;
  std::vector<T> numerical_du_;
  std::vector<T> analytical_u_;
  std::vector<T> analytical_du_;
};

template <NumericType T>
void Differentiator<T>::compute_analy_sols(AnalyticalFunction<T> ana_u,
                                           AnalyticalFunction<T> ana_du, int k) {
  analytical_u_ = std::vector<T>(x_.size());
  analytical_du_ = std::vector<T>(x_.size());
  for (int i = 0; i < x_.size(); i++) {
    analytical_u_[i] = ana_u(x_[i], k);
    analytical_du_[i] = ana_du(x_[i], k);
  }
}

template <NumericType T> T Differentiator<T>::L_inf_error() {
  std::vector<T> abs_diffs(x_.size());
  for (int i = 0; i < x_.size(); i++) {
    abs_diffs[i] = error::absolute_diff(analytical_du_[i], numerical_du_[i]);
  }
  T L_inf = *std::max_element(abs_diffs.begin(), abs_diffs.end());
  return L_inf;
}

template <NumericType T> T Differentiator<T>::L_2_error() {
  std::vector<T> abs_diffs(x_.size());
  for (int i = 0; i < x_.size(); i++) {
    abs_diffs[i] = error::absolute_diff(analytical_du_[i], numerical_du_[i]);
  }
  T L_2 = 0;
  for (int i = 0; i < x_.size(); i++) {
    L_2 += abs_diffs[i] * abs_diffs[i];
  }
  L_2 = sqrt(L_2 / x_.size()); // Normalized by number of points
  return L_2;
}

#endif // INCLUDE_INCLUDE_DIFF_H_
