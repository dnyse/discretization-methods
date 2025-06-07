# Discretization Methods: Burgers' Equation Solver

## Overview

This project implements spectral methods for solving Burgers' equation, focusing on the comparison between Fourier Collocation and Fourier Galerkin approaches. The implementation demonstrates advanced numerical techniques for nonlinear partial differential equations, including dealiasing strategies to handle spectral aliasing errors.

The project addresses three main components:

1. **Theoretical Analysis**: Wellposedness, consistency, and stability of advection-diffusion equations
2. **Fourier Collocation Method**: Direct spatial discretization with collocation points
3. **Fourier Galerkin Method**: Spectral space evolution with careful aliasing control

## Mathematical Problem

The project solves the viscous Burgers' equation:

```
∂u(x,t)/∂t + u(x,t) ∂u(x,t)/∂x = ν ∂²u(x,t)/∂x²
```

with periodic boundary conditions on x ∈ [0, 2π], using:
- Propagation velocity: c = 4.0
- Viscosity coefficient: ν = 0.1
- Analytical solution for validation and error analysis

## Features

### **Spectral Methods Implementation**
- **Fourier Collocation**: Point-wise PDE satisfaction with odd-numbered grid points
- **Fourier Galerkin**: Spectral coefficient evolution with orthogonal projection
- **Adaptive Dealiasing**: 2/3 rule implementation to prevent aliasing instabilities

### **Advanced Numerical Techniques**
- **4th-order Runge-Kutta**: Specialized time integration scheme from exam specification
- **CFL Stability Analysis**: Experimental determination of stability boundaries
- **Spectral Convergence**: High-order accuracy for smooth periodic solutions

### **Comprehensive Analysis Tools**
- **Error Measurement**: L∞ error computation against analytical solutions
- **Convergence Rate**: Quantitative assessment of spectral convergence
- **Stability Testing**: Systematic CFL boundary detection
- **Performance Comparison**: Direct method-to-method evaluation

### **High-Precision Implementation**
- Boost multiprecision arithmetic support
- Robust FFT operations for spectral transformations
- Careful numerical handling of nonlinear terms

## Requirements

- **C++20 compiler** (g++-14 recommended)
- **CMake 3.16** or newer
- **Boost 1.87.0** or newer (Multiprecision, Math, Constants)
- **OpenMP** (optional, for parallel matrix operations)

## Building the Project

```bash
# Configure the project
cmake -B build -DCMAKE_CXX_COMPILER=g++-14

# Build with optimization
cmake --build build --config Release
```

## Running the Project

```bash
# Run all parts
./build/solver --all

# Run specific parts
./build/solver --ex02b  # Fourier Collocation CFL determination
./build/solver --ex02c  # Fourier Collocation convergence study
./build/solver --ex02d  # Fourier Collocation time evolution
./build/solver --ex03b  # Fourier Galerkin CFL determination
./build/solver --ex03c  # Fourier Galerkin convergence study
./build/solver --ex03d  # Method comparison
```

## Implementation Details

### **Part 2: Fourier Collocation Method**

**Grid Configuration:**
```
xⱼ = 2πj/(N+1), j ∈ [0, N]  (odd method)
```

**Time Step Constraint:**
```
Δt ≤ CFL × [max|u(xⱼ)|/Δx + ν/(Δx)²]⁻¹
```

**Key Features:**
- Direct matrix-vector differentiation using Fourier differentiation matrices
- Point-wise satisfaction of Burgers' equation
- Systematic CFL stability boundary detection

### **Part 3: Fourier Galerkin Method**

**Spectral Expansion:**
```
u(x,t) = Σ û_n(t)e^(inx), n ∈ [-N/2, N/2]
```

**Modified Time Step Constraint:**
```
Δt ≤ CFL × [max|u(xⱼ)|k_max + ν(k_max)²]⁻¹
```
where k_max = N/2

**Dealiasing Strategy:**
- **2/3 Rule**: k_cutoff = 2k_max/3
- Applied to nonlinear term computation
- Prevents spectral aliasing instabilities

### **4th-Order Runge-Kutta Integration**

Specialized scheme as specified in the exam:
```
u₁ = uⁿ + (Δt/2)F(uⁿ)
u₂ = uⁿ + (Δt/2)F(u₁)  
u₃ = uⁿ + ΔtF(u₂)
uⁿ⁺¹ = (1/3)(-uⁿ + u₁ + 2u₂ + u₃ + (Δt/2)F(u₃))
```

where F(u) = -u∂u/∂x + ν∂²u/∂x²

## Code Structure

### **Core Classes**
- **`BurgersSolver`**: Fourier Collocation implementation
- **`BurgersGalerkinSolver`**: Fourier Galerkin implementation  
- **`FourierGalerkin`**: Spectral space operations with FFT
- **`SpectralFourier`**: Differentiation matrix construction
- **`RungeKutta4`**: Time integration with Burgers-specific interface

### **Key Files**
- **`burger.h/cpp`**: Main solver implementations
- **`fourier.h`**: Spectral differentiation and Galerkin methods
- **`integrate.h`**: Time integration schemes
- **`common.h`**: Exact solutions and utility functions
- **`main.cc`**: Command-line interface and test execution

### **Supporting Infrastructure**
- **Error Analysis**: L∞ error computation against analytical solutions
- **CFL Testing**: Stability boundary detection with error spike monitoring
- **Visualization**: Solution plotting and convergence analysis
- **Performance Timing**: Computational efficiency measurement

## Theoretical Background

### **Burgers' Equation Analytics**

The exact solution is given by:
```
u(x,t) = c - 2ν ∂φ(x-ct,t+1)/∂x / φ(x-ct,t+1)

φ(a,b) = Σ exp[-(a-(2k+1)π)²/(4νb)]
```

This represents a sawtooth-like traveling wave with viscous smoothing.

### **Spectral Method Theory**

**Fourier Collocation:**
- Interpolation-based approach
- Direct point-wise PDE satisfaction
- Infinite-order finite difference equivalent

**Fourier Galerkin:**
- Orthogonal projection onto trigonometric polynomials
- Evolution in spectral coefficient space
- Natural conservation properties

### **Aliasing in Nonlinear Terms**

When computing u∂u/∂x:
- Product generates frequencies up to 2k_max
- Exceeds spectral resolution → aliasing
- 2/3 dealiasing rule: retain only |k| ≤ 2k_max/3

## Expected Results

### **Fourier Collocation**
- **CFL Values**: 0.65 - 1.40 (decreasing with grid refinement)
- **Convergence**: Spectral rates (6-8) for smooth solutions
- **Stability**: Consistent, predictable behavior

### **Fourier Galerkin**  
- **CFL Values**: 4.15 - 8.05 (significantly higher stability limits)
- **Convergence**: Mixed behavior due to aliasing sensitivity
- **Performance**: Potential efficiency gains from larger time steps

### **Method Comparison**
- **Accuracy**: Collocation typically more accurate
- **Stability**: Galerkin allows larger CFL values
- **Implementation**: Galerkin requires sophisticated dealiasing
