#ifndef INCLUDE_INCLUDE_INTEGRATE_H_
#define INCLUDE_INCLUDE_INTEGRATE_H_

#include "diff.h"

template <NumericType T> class Integrator {
public:
  Integrator() = default;
  virtual ~Integrator() = default;

  virtual std::vector<T> integrate(const std::vector<T> &u, T dt,
                                   Differentiator<T> &differentiator) = 0;

protected:
  std::vector<std::vector<T>> history_numerical_u_ = {};
};

template <NumericType T> class RungeKutta4 : public Integrator<T> {
public:
  std::vector<T> integrate(const std::vector<T> &u, T dt,
                           Differentiator<T> &differentiator) override;
  std::vector<T> integrate_burgers(const std::vector<T> &u, T dt,
                                   Differentiator<T> &differentiator, T nu);
};

template <NumericType T>
std::vector<T> RungeKutta4<T>::integrate(const std::vector<T> &u, T dt,
                                         Differentiator<T> &differentiator) {
  size_t N = u.size() - 1;
  std::vector<T> u1(N + 1), u2(N + 1), u3(N + 1), un_p1(N + 1);

  // u1 = u + dt/2 * F(u)
  std::vector<T> k1 = differentiator.compute_F(const_cast<std::vector<T> &>(u));
  for (size_t i = 0; i <= N; ++i) {
    u1[i] = u[i] + dt / T(2) * k1[i];
  }

  // u2 = u + dt/2 * F(u1)
  std::vector<T> k2 = differentiator.compute_F(u1);
  for (size_t i = 0; i <= N; ++i) {
    u2[i] = u[i] + dt / T(2) * k2[i];
  }

  // u3 = u + dt * F(u2)
  std::vector<T> k3 = differentiator.compute_F(u2);
  for (size_t i = 0; i <= N; ++i) {
    u3[i] = u[i] + dt * k3[i];
  }

  // un+1 = (1/3) * (-u + u1 + 2*u2 + u3 + dt/2 * F(u3))
  std::vector<T> k4 = differentiator.compute_F(u3);
  for (size_t i = 0; i <= N; ++i) {
    un_p1[i] = (T(1) / T(3)) *
               (-u[i] + u1[i] + T(2) * u2[i] + u3[i] + dt / T(2) * k4[i]);
  }

  // Store for history if needed
  this->history_numerical_u_.push_back(un_p1);

  return un_p1;
}

template <NumericType T>
std::vector<T>
RungeKutta4<T>::integrate_burgers(const std::vector<T> &u, T dt,
                                  Differentiator<T> &differentiator, T nu) {
  size_t N = u.size() - 1;
  std::vector<T> u1(N + 1), u2(N + 1), u3(N + 1), un_p1(N + 1);

  // RK4 for Burgers' equation: du/dt = -u*du/dx + ν*d²u/dx²

  // k1 = f(u)
  std::vector<T> u_copy = u;
  std::vector<T> k1 = differentiator.compute_burgers_rhs(u_copy, nu);
  for (size_t i = 0; i <= N; ++i) {
    u1[i] = u[i] + dt / T(2) * k1[i];
  }

  // k2 = f(u + dt/2 * k1)
  std::vector<T> k2 = differentiator.compute_burgers_rhs(u1, nu);
  for (size_t i = 0; i <= N; ++i) {
    u2[i] = u[i] + dt / T(2) * k2[i];
  }

  // k3 = f(u + dt/2 * k2)
  std::vector<T> k3 = differentiator.compute_burgers_rhs(u2, nu);
  for (size_t i = 0; i <= N; ++i) {
    u3[i] = u[i] + dt * k3[i];
  }

  // k4 = f(u + dt * k3)
  std::vector<T> k4 = differentiator.compute_burgers_rhs(u3, nu);


  // un+1 = u + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
  for (size_t i = 0; i <= N; ++i) {
  un_p1[i] = (T(1) / T(3)) *
               (-u[i] + u1[i] + T(2) * u2[i] + u3[i] + dt / T(2) * k4[i]);
  }

  // Store for history if needed
  this->history_numerical_u_.push_back(un_p1);

  return un_p1;
}

#endif // INCLUDE_INCLUDE_INTEGRATE_H_
