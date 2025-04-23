# Numerical Differentiation Solver

## Project Overview

This project implements various numerical differentiation methods to study their accuracy, convergence, and performance characteristics when applied to different types of functions and a scalar hyperbolic PDE. The code leverages C++20 features, high-precision arithmetic through Boost's multiprecision library, and parallel computation using OpenMP.

## Features

- **Multiple Differentiation Methods**:
  - Second-order finite difference (central difference)
  - Fourth-order finite difference
  - Spectral differentiation using Fourier matrices (both EVEN and ODD variants)

- **High-Precision Arithmetic**:
  - Uses Boost's `cpp_dec_float_50` type for 50 decimal digits of precision
  - Ensures accurate numerical results for high-precision convergence studies
  - Template-based implementation allowing different numeric types

- **Error Analysis**:
  - L_inf (maximum) error norm
  - L_2 (root mean square) error norm
  - Relative and absolute error metrics
  - Convergence rate calculations

- **Hyperbolic PDE Solver**:
  - Implements a scalar hyperbolic equation solver
  - Uses fourth-order Runge-Kutta time integration
  - Periodic boundary conditions: u(0,t) = u(2π,t)
  - Initial condition: u(x,0) = exp(sin(x))
  - Analytical solution: u(x,t) = exp(sin(x-2πt))

- **Performance Optimization**:
  - OpenMP parallelization for matrix operations
  - Block-based approach for improved cache efficiency with large matrices

- **Visualization**:
  - Integration with Matplot++ for solution visualization
  - Automated figure generation

## Requirements

- C++20 compatible compiler (GCC 10+ recommended)
- CMake 3.16 or newer
- Boost 1.87.0 or newer (specifically Boost.Multiprecision and Boost.Math)
- OpenMP support
- Matplot++ (automatically fetched by CMake during build)

## Building

```bash
# Configure the project
cmake -B build -S . -DCMAKE_CXX_COMPILER=g++-14

# Build the project
cmake --build build
```

## Running

The application provides a command-line interface to run specific exercises:

```bash
# Run all exercises
./build/solver --all

# Run specific exercises
./build/solver --ex01  # EVEN vs ODD Fourier differentiation
./build/solver --ex02  # Fourier differentiation convergence
./build/solver --ex03a  # Hyperbolic PDE convergence study
./build/solver --ex03b  # Hyperbolic PDE long-term integration

# Display help
./build/solver --help
```
## Code Structure
### Core Components

- **Core Utilities and Types**:
  - `common.h`: Defines numeric types, constants, and test functions
  - `error.h`: Error calculation functions

- **Differentiation Framework**:
  - `diff.h`: Base differentiator interface
  - `finite_diff.h`: Second and fourth-order finite difference implementations
  - `fourier.h`: Spectral Fourier differentiation implementation

- **PDE Solver**:
  - `integrate.h`: Time integration methods (RK4)
  - `hyperbolic.h`: Hyperbolic PDE solver implementation

- **Application**:
  - `main.cc`: Main program with command-line interface

### Class Hierarchy

- `Differentiator<T>`: Abstract base class for differentiation methods
  - `SecondOrderFiniteDiff<T>`: Second-order finite difference implementation
  - `FourthOrderFiniteDiff<T>`: Fourth-order finite difference implementation
  - `SpectralFourier<T>`: Fourier spectral differentiation implementation

- `Integrator<T>`: Abstract base class for time integration methods
  - `RungeKutta4<T>`: Fourth-order Runge-Kutta implementation

- `HyperbolicSolver<T>`: Solver for the scalar hyperbolic PDE
