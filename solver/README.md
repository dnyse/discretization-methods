# Numerical Methods for Differentiation

## Overview

This project implements various numerical differentiation methods to analyze their accuracy and performance characteristics. The implementation focuses on three main techniques:

1. Finite difference methods (2nd and 4th order)
2. Spectral Fourier differentiation
3. Application to a scalar hyperbolic PDE

The project is designed to study convergence rates, error behavior, and effectiveness of different numerical differentiation schemes for various functions, with particular attention to periodic boundary conditions.

## Features

- **Multiple Differentiation Methods**:
  - Second-order finite difference
  - Fourth-order finite difference
  - Spectral differentiation using Fourier matrices (EVEN and ODD variants)

- **High-Precision Arithmetic**:
  - Uses Boost's multiprecision library (50 decimal digits)
  - Ensures accurate numerical results for convergence studies

- **Error Analysis**:
  - L∞ (maximum) and L₂ (root mean square) error measurements
  - Convergence rate calculations
  - Relative and absolute error metrics

- **PDE Solver**:
  - Implementation of a scalar hyperbolic equation solver
  - Fourth-order Runge-Kutta time integration
  - Comparative study of different spatial discretization methods

## Requirements

- C++20 compiler (g++-14 recommended)
- CMake 3.16 or newer
- Boost 1.87.0 or newer (specifically Boost.Multiprecision and Boost.Math)

## Building the Project

```bash
# Configure the project
cmake -B build -DCMAKE_CXX_COMPILER=g++-14

# Build the project
cmake --build build
```

## Running the Project

```bash
./build/solver
```

The program will run three exercises that demonstrate different aspects of numerical differentiation:

1. **Exercise 1**: Compares ODD and EVEN Fourier differentiation matrices for the function u(x) = exp(k·sin(x)) with various k values
2. **Exercise 2**: Analyzes convergence rates of Fourier differentiation for three test functions
3. **Exercise 3**: Investigates a scalar hyperbolic PDE using different discretization methods

## Exercise Details

### Exercise 1: Fourier Differentiation (EVEN vs ODD)

This exercise compares two different Fourier differentiation matrix formulations (EVEN and ODD) applied to the function u(x) = exp(k·sin(x)) for k = 2, 4, 6, 8, 10, 12.

For each k value, the program finds the minimum grid size N needed to achieve a relative error below 10⁻⁵.

### Exercise 2: Fourier Differentiation Convergence

This exercise studies the convergence behavior of the EVEN Fourier method for three different functions:

1. f(x) = cos(10x) - A well-behaved periodic function
2. f(x) = cos(x/2) - A function that is not periodic on [0,2π]
3. f(x) = x - A linear non-periodic function

The analysis demonstrates how spectral methods perform differently based on the properties of the function being differentiated, particularly periodicity.

### Exercise 3: Scalar Hyperbolic Problem

This exercise solves the scalar hyperbolic PDE:

∂u(x,t)/∂t = -2π ∂u(x,t)/∂x

with periodic boundary conditions u(0,t) = u(2π,t) and initial condition u(x,0) = exp(sin(x)).

The exercise is divided into two parts:

1. **Part (a)**: Convergence study comparing second-order finite difference, fourth-order finite difference, and Fourier spectral differentiation at time t = π.

2. **Part (b)**: Long-time integration comparison between second-order finite difference with N=200 and Fourier spectral method with N=10, evaluated at times t = 0, 100, and 200.

## Code Structure

- **Common**: Core data types and utility functions
  - `common.h`: Defines numeric types, mathematical constants, and test functions

- **Differentiation**: Base classes and implementations
  - `diff.h`: Base differentiator interface
  - `finite_diff.h`: Finite difference implementations
  - `fourier.h`: Spectral Fourier implementation
  - `error.h`: Error calculation utilities

- **PDE Solver**:
  - `integrate.h`: Time integration methods including RK4
  - `hyperbolic.h`: Hyperbolic PDE solver implementation

- **Application**:
  - `main.cc`: Main program that runs all exercises

## Mathematical Background

### Numerical Differentiation Methods

1. **Second-order Finite Difference**:
   ```
   f'(x_j) ≈ (f(x_{j+1}) - f(x_{j-1}))/(2Δx)
   ```

2. **Fourth-order Finite Difference**:
   ```
   f'(x_j) ≈ (-f(x_{j+2}) + 8f(x_{j+1}) - 8f(x_{j-1}) + f(x_{j-2}))/(12Δx)
   ```

3. **Fourier Differentiation Matrix**:
   - ODD variant: D_{j,i} = (-1)^{j+i} / (2sin((j-i)π/(N+1))), for i≠j
   - EVEN variant: D_{j,i} = (-1)^{j+i} / (2tan((j-i)π/N)), for i≠j

### Hyperbolic PDE

The scalar hyperbolic problem represents a simple wave traveling with constant velocity 2π:
- Analytical solution: u(x,t) = exp(sin(x-2πt))

The 4th-order Runge-Kutta time integration is defined as:
```
u₁ = uⁿ + (Δt/2)F(uⁿ)
u₂ = uⁿ + (Δt/2)F(u₁)
u₃ = uⁿ + ΔtF(u₂)
uⁿ⁺¹ = (1/3)(-uⁿ + u₁ + 2u₂ + u₃ + (Δt/2)F(u₃))
```
where F(u) = -2π ∂u/∂x.

## Results and Analysis

### Exercise 1
Comparing EVEN and ODD methods shows which formulation requires fewer grid points to achieve the desired accuracy for functions with different oscillatory behavior.

### Exercise 2
The convergence analysis demonstrates:
- Spectral accuracy (exponential convergence) for smooth periodic functions
- Reduced effectiveness for non-periodic or discontinuous functions
- The importance of matching the method to the problem characteristics

### Exercise 3
The hyperbolic PDE solver comparison reveals:
- Convergence rates for different spatial discretization methods
- Trade-offs between high-order methods and computational efficiency
- Long-time integration behavior, particularly for spectral vs. finite difference methods
