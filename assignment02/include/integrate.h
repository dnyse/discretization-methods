#ifndef INCLUDE_INCLUDE_INTEGRATE_H_
#define INCLUDE_INCLUDE_INTEGRATE_H_

#include "diff.h"

template <NumericType T> class Integrator {
public:
  Integrator() = default;
  virtual ~Integrator() = default;

  virtual std::vector<T> integrate(const std::vector<T> &u, T dt,
                                   Differentiator<T> &differentiator) = 0;  // Changed to non-const reference

protected:
  std::vector<std::vector<T>> history_numerical_u_ = {};
};

template <NumericType T> class RungeKutta4 : public Integrator<T> {
public:
  std::vector<T> integrate(const std::vector<T> &u, T dt,
                           Differentiator<T> &differentiator) override;  // Changed to non-const reference
};

template <NumericType T>
std::vector<T>
RungeKutta4<T>::integrate(const std::vector<T> &u, T dt,
                          Differentiator<T> &differentiator) {  // Changed to non-const reference
  size_t N = u.size() - 1;
  std::vector<T> u1(N + 1), u2(N + 1), u3(N + 1), un_p1(N + 1);

  // u1 = u + dt/2 * F(u)
  std::vector<T> k1 = differentiator.compute_F(const_cast<std::vector<T>&>(u));  // Use non-const version
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

#endif // INCLUDE_INCLUDE_INTEGRATE_H_
