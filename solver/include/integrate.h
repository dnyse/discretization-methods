#ifndef INCLUDE_INCLUDE_INTEGRATE_H_
#define INCLUDE_INCLUDE_INTEGRATE_H_

#include "diff.h"

template <NumericType T> class Integrator {
public:
  Integrator() = default;
  virtual ~Integrator() = default;

  virtual std::vector<T> integrate(const std::vector<double> &u, double dt,
                                   const Differentiator<T> &differentiator) = 0;

protected:
  std::vector<std::vector<T>> history_numerical_u_ = {};
};

template <NumericType T> class RungeKutta4 : public Integrator<T> {
public:
  std::vector<T> integrate(const std::vector<double> &u, double dt,
                           const Differentiator<T> &differentiator) override;
};

template <NumericType T>
std::vector<T>
RungeKutta4<T>::integrate(const std::vector<double> &u, double dt,
                          const Differentiator<T> &differentiator) {
  size_t N = u.size() - 1;
  std::vector<double> u1(N + 1), u2(N + 1), u3(N + 1), un_p1(N + 1);

  // u1 = u + dt/2 * F(u)
  std::vector<double> k1 = differentiator.F(u);
  for (size_t i = 0; i <= N; ++i) {
    u1[i] = u[i] + dt / 2.0 * k1[i];
  }

  // u2 = u + dt/2 * F(u1)
  std::vector<double> k2 = differentiator.F(u1);
  for (size_t i = 0; i <= N; ++i) {
    u2[i] = u[i] + dt / 2.0 * k2[i];
  }

  // u3 = u + dt * F(u2)
  std::vector<double> k3 = differentiator.F(u2);
  for (size_t i = 0; i <= N; ++i) {
    u3[i] = u[i] + dt * k3[i];
  }

  // un+1 = (1/3) * (-u + u1 + 2*u2 + u3 + dt/2 * F(u3))
  std::vector<double> k4 = differentiator.F(u3);
  for (size_t i = 0; i <= N; ++i) {
    un_p1[i] =
        (1.0 / 3.0) * (-u[i] + u1[i] + 2.0 * u2[i] + u3[i] + dt / 2.0 * k4[i]);
  }
  this->history_numerical_u_.push_back(un_p1);
  return un_p1;
}

#endif // INCLUDE_INCLUDE_INTEGRATE_H_
