#ifndef INCLUDE_INCLUDE_FOURIER_H_
#define INCLUDE_INCLUDE_FOURIER_H_

#include "diff.h"
#include "error.h"
#include <iostream>
#include <omp.h>

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
  int size = this->analytical_u_.size();
  this->numerical_du_ = std::vector<T>(size, T(0));

  // For small matrices, use direct matrix-vector multiplication
  if (size <= 1000) {
#pragma omp parallel for
    for (int j = 0; j < size; j++) {
      T sum = T(0);
      for (int i = 0; i < size; i++) {
        sum += D_[j][i] * this->analytical_u_[i];
      }
      this->numerical_du_[j] = sum;
    }
  }
  // For larger matrices, use block-based approach to improve cache efficiency
  else {
    const int block_size = 64; // Adjust based on cache size

#pragma omp parallel for
    for (int j = 0; j < size; j++) {
      T sum = T(0);
      for (int i_block = 0; i_block < size; i_block += block_size) {
        int i_end = std::min(i_block + block_size, size);
        for (int i = i_block; i < i_end; i++) {
          sum += D_[j][i] * this->analytical_u_[i];
        }
      }
      this->numerical_du_[j] = sum;
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

template <NumericType T> class FourierGalerkin : public Differentiator<T> {
public:
  FourierGalerkin() {}

  void create_grid_pts(int N) override {
    N_ = N;
    int num_pts = N + 1;
    this->x_ = std::vector<T>(num_pts);
    T dx = T(2) * MathConstants<T>::PI() / T(num_pts);

    for (int j = 0; j < num_pts; ++j) {
      this->x_[j] = dx * T(j);
    }
  }

  void compute_num_sol() override {
    // Not needed for Galerkin - we work in spectral space
  }

  // Compute Fourier coefficients from physical space values
  std::vector<std::complex<T>> fft(const std::vector<T> &u) {
    int N = u.size() - 1;
    std::vector<std::complex<T>> u_hat(N + 1);

    for (int k = -N / 2; k <= N / 2; ++k) {
      std::complex<T> sum(0, 0);
      for (int j = 0; j <= N; ++j) {
        T arg = -T(2) * MathConstants<T>::PI() * T(k) * T(j) / T(N + 1);
        std::complex<T> exp_term(cos(arg), sin(arg));
        sum += u[j] * exp_term;
      }
      int idx = (k + N / 2 + N + 1) % (N + 1);
      u_hat[idx] = sum / T(N + 1);
    }
    return u_hat;
  }

  // Compute physical space values from Fourier coefficients
  std::vector<T> ifft(const std::vector<std::complex<T>> &u_hat) {
    int N = u_hat.size() - 1;
    std::vector<T> u(N + 1);

    for (int j = 0; j <= N; ++j) {
      std::complex<T> sum(0, 0);
      for (int k = -N / 2; k <= N / 2; ++k) {
        int idx = (k + N / 2 + N + 1) % (N + 1);
        T arg = T(2) * MathConstants<T>::PI() * T(k) * T(j) / T(N + 1);
        std::complex<T> exp_term(cos(arg), sin(arg));
        sum += u_hat[idx] * exp_term;
      }
      u[j] = sum.real();
    }
    return u;
  }

  // Compute spectral derivative
  std::vector<std::complex<T>>
  spectral_derivative(const std::vector<std::complex<T>> &u_hat,
                      int order = 1) {
    int N = u_hat.size() - 1;
    std::vector<std::complex<T>> du_hat(N + 1);

    for (int k = -N / 2; k <= N / 2; ++k) {
      int idx = (k + N / 2 + N + 1) % (N + 1);
      std::complex<T> ik(0, T(k));
      du_hat[idx] = pow(ik, order) * u_hat[idx];
    }
    return du_hat;
  }

  // Compute RHS for Burgers equation in spectral space
  std::vector<std::complex<T>>
  compute_burgers_spectral_rhs(const std::vector<T> &u, T nu) {
    // Transform to spectral space
    auto u_hat = fft(u);

    // Compute spectral derivatives
    auto du_hat = spectral_derivative(u_hat, 1);
    auto d2u_hat = spectral_derivative(u_hat, 2);

    // Transform derivatives back to physical space
    auto du = ifft(du_hat);

    // Compute nonlinear term u * du/dx in physical space
    std::vector<T> nonlinear(u.size());
    for (size_t i = 0; i < u.size(); ++i) {
      nonlinear[i] = u[i] * du[i];
    }

    // Transform nonlinear term to spectral space
    auto nonlinear_hat = fft(nonlinear);

    // Compute full RHS in spectral space: -transform(u*du/dx) + nu*d2u/dx2
    std::vector<std::complex<T>> rhs_hat(u_hat.size());
    for (size_t i = 0; i < u_hat.size(); ++i) {
      rhs_hat[i] = -nonlinear_hat[i] + nu * d2u_hat[i];
    }

    return rhs_hat;
  }

private:
  int N_;
};

#endif // INCLUDE_INCLUDE_FOURIER_H_
