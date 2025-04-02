#ifndef INCLUDE_INCLUDE_FOURIER_H_
#define INCLUDE_INCLUDE_FOURIER_H_

#include "diff.h"
#include "error.h"
#include <iostream>

template <NumericType T> class SpectralFourier : public Differentiator<T> {
public:
  SpectralFourier(MethodType method);
  void build(int N);
  int get_num_grd_pts(int N);
  T calc_rel_error(AnalyticalFunction<T> ana_u, AnalyticalFunction<T> ana_du,
                   int N, int k = 0);
  void compute_num_sol() override;
  void create_grid_pts(int N) override;
  std::vector<std::pair<int, T>> find_min_N(const std::vector<int> &k_values,
                                            AnalyticalFunction<T> ana_u,
                                            AnalyticalFunction<T> ana_du,
                                            T error_threshold = T(1e-5),
                                            int max_N = 100);
  std::tuple<T, T> compute_errors(int N, AnalyticalFunction<T> ana_u,
                                  AnalyticalFunction<T> ana_du);

private:
  MethodType method_;
  std::vector<std::vector<T>> D_;
};

template <NumericType T>
SpectralFourier<T>::SpectralFourier(MethodType method) : method_(method) {}

template <NumericType T> void SpectralFourier<T>::build(int N) {
  int grid_pts = get_num_grd_pts(N);
  D_ = std::vector<std::vector<T>>(grid_pts, std::vector<T>(grid_pts, T(0)));
  if (method_ == MethodType::ODD) {
    for (int j = 0; j < grid_pts; j++) {
      for (int i = 0; i < grid_pts; i++) {
        if (i != j) {
          D_[j][i] =
              (pow(T(-1), j + i)) /
              (T(2) * sin((j - i) * MathConstants<T>::PI() / (grid_pts)));
        }
      }
    }
  } else {
    for (int j = 0; j < grid_pts; j++) {
      for (int i = 0; i < grid_pts; i++) {
        if (i != j) {
          D_[j][i] =
              (pow(T(-1), j + i)) /
              (T(2) * tan((j - i) * MathConstants<T>::PI() / (grid_pts)));
        }
      }
    }
  }
}

template <NumericType T> int SpectralFourier<T>::get_num_grd_pts(int N) {
  return method_ == MethodType::ODD ? N + 1 : N;
}

template <NumericType T> void SpectralFourier<T>::create_grid_pts(int N) {
  int grid_pts = get_num_grd_pts(N);
  this->x_ = std::vector<T>(grid_pts);
  T step = T(2) * MathConstants<T>::PI() / (grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    this->x_[i] = i * step;
  }
}

template <NumericType T> void SpectralFourier<T>::compute_num_sol() {
  this->numerical_du_ = std::vector<T>(this->x_.size());
  for (int j = 0; j < this->analytical_u_.size(); j++) {
    for (int i = 0; i < this->analytical_u_.size(); i++) {
      this->numerical_du_[j] += D_[j][i] * this->analytical_u_[i];
    }
  }
}

template <NumericType T>
T SpectralFourier<T>::calc_rel_error(AnalyticalFunction<T> ana_u,
                                     AnalyticalFunction<T> ana_du, int N,
                                     int k) {
  int grid_pts = get_num_grd_pts(N);
  create_grid_pts(N);
  this->compute_analy_sols(ana_u, ana_du, k);
  build(N);
  compute_num_sol();

  std::vector<T> errors(grid_pts);
  for (int i = 0; i < grid_pts; i++) {
    errors[i] =
        error::relative(this->analytical_du_[i], this->numerical_du_[i]);
  }

  return *std::max_element(errors.begin(), errors.end());
}

template <NumericType T>
std::vector<std::pair<int, T>> SpectralFourier<T>::find_min_N(
    const std::vector<int> &k_values, AnalyticalFunction<T> ana_u,
    AnalyticalFunction<T> ana_du, T error_threshold, int max_N) {
  std::vector<std::pair<int, T>> results;
  for (int k : k_values) {
    int best_N = -1;
    T best_error = T(1.0);

    for (int N = 4; N <= max_N; N += 2) {
      T max_error = calc_rel_error(ana_u, ana_du, N, k);
      if (max_error < best_error) {
        best_error = max_error;
        best_N = N;
        if (best_error < error_threshold) {
          std::cout << "Found best N=" << N << " for k=" << k
                    << " with error=" << best_error << std::endl
                    << std::endl;
          break;
        }
      }
    }
    results.emplace_back(best_N, best_error);
  }
  return results;
}

template <NumericType T>
std::tuple<T, T>
SpectralFourier<T>::compute_errors(int N, AnalyticalFunction<T> ana_u,
                                   AnalyticalFunction<T> ana_du) {
  int grid_pts = get_num_grd_pts(N);
  create_grid_pts(N);
  this->compute_analy_sols(ana_u, ana_du);
  build(N);
  compute_num_sol();
  T L_inf = this->L_inf_error();
  T L_2 = this->L_2_error();
  return std::make_tuple(L_inf, L_2);
}

#endif // INCLUDE_INCLUDE_FOURIER_H_
