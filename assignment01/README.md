# Fourier Differentiation Matrix Solver

## Overview

This project implements a Fourier differentiation matrix solver to analyze the accuracy of spectral differentiation for the function \( u(x) = \exp(k \sin x) \). The implementation uses high-precision arithmetic via Boost's multiprecision library to ensure accurate numerical results.

## Features

- Constructs Fourier differentiation matrices for spectral differentiation
- Calculates numerical derivatives using spectral methods
- Compares numerical derivatives with analytical solutions
- Determines the minimum grid size (\( N \)) required to achieve a specified error threshold
- Uses high-precision arithmetic (50 decimal digits) for accurate numerical computations

## Requirements

- C++20 compiler
- CMake 3.16 or newer
- Boost 1.87.0 or newer
  - Specifically requires Boost.Multiprecision and Boost.Math

## Building the Project

```bash
# Configure the project with g++-14 compiler
cmake -B build -DCMAKE_CXX_COMPILER=g++-14

# Build the project
cmake --build build
```

## Usage

Run the executable to calculate the minimum grid size (\( N \)) needed to achieve a relative error of less than \( 10^{-5} \) for different values of \( k \):

```bash
./solver
```

The program will output:
- The minimum \( N \) required for each \( k \) value (\( k = 2, 4, 6, 8, 10, 12 \))
- The maximum relative error achieved with that \( N \)

## Technical Details

### Mathematical Background

The program uses spectral differentiation with Fourier differentiation matrices to approximate derivatives of the function:

\[
u(x) = \exp(k \sin x)\]

The analytical derivative is:

\[
u'(x) = k\cos(x) \exp(k\sin x)\]

### Implementation

- The differentiation matrix \( D \) is constructed using the formula:

  \[
  D_{j,i} = \frac{(-1)^{j+i}}{2 \sin\left(\frac{(j-i)\pi}{N+1}\right)}, \quad \text{for } i \neq j
  \]

- The relative error is calculated as:

  \[
  \text{error} = \frac{|\text{numerical derivative} - \text{analytical derivative}|}{|\text{analytical derivative}|}
  \]

- The program finds the minimum \( N \) for which the maximum relative error falls below \( 10^{-5} \).

## Code Structure

- `main.cc`: Contains the main implementation including:
  - Fourier differentiation matrix construction
  - Analytical and numerical derivative calculations
  - Error analysis and grid size optimization

