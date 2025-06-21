# Discretization Methods

A comprehensive collection of numerical methods implementations for solving partial differential equations (PDEs) using advanced discretization techniques. This repository demonstrates the implementation and analysis of spectral methods, finite difference schemes, and time integration methods for various classes of PDEs.

## Overview

This project explores the numerical solution of PDEs through three main components:

1. **Assignment 01**: Fourier Differentiation Matrix Analysis
2. **Assignment 02**: Comparative Study of Numerical Differentiation Methods
3. **Exam Project**: Burgers' Equation Spectral Solvers

Each component builds upon fundamental numerical analysis concepts while implementing increasingly sophisticated algorithms for solving different types of PDEs.

## Repository Structure

```
discretization-methods/
├── assignment01/          # Fourier differentiation matrix solver
├── assignment02/          # Multi-method differentiation framework
├── exam/                  # Burgers' equation spectral methods
└── README.md             # This file
```

## Key Features

### **High-Precision Numerical Computing**
- Boost multiprecision arithmetic (50 decimal digits)
- Robust error analysis and convergence studies
- Template-based implementation supporting multiple numeric types

### **Advanced Discretization Methods**
- **Spectral Methods**: Fourier collocation and Galerkin approaches
- **Finite Difference**: Second and fourth-order accurate schemes
- **Matrix-Free Implementations**: Efficient spectral differentiation

### **Comprehensive PDE Solvers**
- **Hyperbolic PDEs**: Advection equations with periodic boundary conditions
- **Burgers' Equation**: Nonlinear viscous conservation laws
- **Time Integration**: Fourth-order Runge-Kutta with stability analysis

### **Performance Optimization**
- OpenMP parallelization for matrix operations
- FFT-based spectral transforms
- Block-based algorithms for cache efficiency

## Quick Start

### Prerequisites

- **C++20 compiler** (GCC 10+ or equivalent)
- **CMake 3.16** or newer
- **Boost 1.87.0** or newer
- **OpenMP** (optional, for parallel execution)

### Building

Each project can be built independently:

```bash
# Assignment 01 - Fourier Differentiation Matrix
cd assignment01
cmake -B build -DCMAKE_CXX_COMPILER=g++-14
cmake --build build
./build/solver

# Assignment 02 - Multi-Method Framework
cd assignment02
cmake -B build -DCMAKE_CXX_COMPILER=g++-14
cmake --build build
./build/solver --all

# Exam Project - Burgers' Equation
cd exam/solver
cmake -B build -DCMAKE_CXX_COMPILER=g++-14
cmake --build build
./build/solver --all
```

## Project Components

### Assignment 01: Fourier Differentiation Matrix Solver

**Focus**: Spectral accuracy analysis for smooth periodic functions

**Key Features**:
- High-precision Fourier differentiation matrix construction
- Error analysis for `u(x) = exp(k sin x)`
- Grid size optimization for specified error tolerances
- Demonstrates spectral convergence rates

**Mathematical Problem**:
```
Find minimum N such that max|Du - u'| < 10⁻⁵
where D is the Fourier differentiation matrix
```

### Assignment 02: Comprehensive Differentiation Framework

**Focus**: Comparative analysis of numerical differentiation methods

**Key Features**:
- Multiple differentiation schemes (finite difference, spectral)
- Hyperbolic PDE solver with time evolution
- Error norm calculations (L∞, L₂)
- Convergence rate analysis across different function classes

**Mathematical Problems**:
1. **Function Differentiation**: Various test functions with known derivatives
2. **Hyperbolic PDE**: `∂u/∂t + 2π ∂u/∂x = 0` with `u(x,0) = exp(sin x)`

### Exam Project: Burgers' Equation Spectral Methods

**Focus**: Advanced spectral methods for nonlinear PDEs

**Key Features**:
- Fourier collocation vs. Fourier Galerkin comparison
- Dealiasing strategies for nonlinear terms
- CFL stability analysis
- Performance benchmarking

**Mathematical Problem**:
```
∂u/∂t + u ∂u/∂x = ν ∂²u/∂x²
```
with analytical solution validation and spectral aliasing control.

## Numerical Methods Implemented

### **Spectral Methods**
- **Fourier Collocation**: Direct PDE satisfaction at collocation points
- **Fourier Galerkin**: Orthogonal projection in spectral space
- **Dealiasing**: 2/3 rule for nonlinear term treatment

### **Finite Difference Methods**
- **Second-Order Central**: Standard finite difference schemes
- **Fourth-Order Compact**: Higher-order accuracy with periodic boundaries

### **Time Integration**
- **Fourth-Order Runge-Kutta**: Specialized scheme for PDE systems
- **CFL Stability Analysis**: Experimental determination of stability limits

## Research Applications

### **Convergence Analysis**
- Spectral convergence rates for smooth periodic functions
- Finite difference convergence for different function classes
- Performance comparison across method types

### **Stability Studies**
- CFL condition determination for different discretizations
- Long-time integration stability
- Method-specific stability characteristics

### **Error Analysis**
- Multiple error norms (L∞, L₂, relative)
- Systematic grid refinement studies
- Method accuracy comparison

## Performance Characteristics

### **Computational Efficiency**
- **Spectral Methods**: Superior for smooth problems, higher setup cost
- **Finite Difference**: Lower memory, consistent performance
- **Parallel Scaling**: OpenMP acceleration for large problems

### **Accuracy Trade-offs**
- **Fourier Methods**: Exponential convergence for periodic problems
- **Finite Difference**: Algebraic convergence, robust for non-smooth data
- **Time Integration**: Balance between stability and accuracy
