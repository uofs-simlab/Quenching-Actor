# Quenching Actor System

A high-performance computational system for studying quenching phenomena using the C++ Actor Framework (CAF). This project implements distributed computational actors to efficiently explore parameter spaces and detect critical boundaries in nonlinear dynamical systems.

## Overview

The Quenching Actor System is designed to study critical phenomena in reaction-diffusion systems, particularly focusing on the FitzHugh-Nagumo model. The system uses actor-based parallelization to distribute computational workloads across multiple nodes and cores.

## Key Features

- **Actor-based Architecture**: Built using the C++ Actor Framework (CAF) for scalable distributed computing
- **Adaptive Refinement**: Dynamic parameter space exploration with boundary detection
- **Multiple Solver Strategies**: Both dynamic and static regridding approaches
- **Performance Optimizations**: Bracket optimization and early out strategies

## Project Structure

```
.
├── cpp_implementation/          # Main C++ implementation
│   ├── include/                # Header files
│   ├── src/                    # Source files  
│   └── README.md              # Detailed implementation docs
├── quenching_with_julia/       # CAF+Julia implementation
├── experimental/               # Julia code
├── job_comparison_analysis.py  # Comparing results
├── plot_data.py               # Visualization utilities
└── waves/                     # Wave data and initial conditions
```

## Dependencies

- **CAF (C++ Actor Framework)**
- **SUNDIALS**
- **Eigen3**
- **GSL**

## Quick Start

### Building the System

```bash
cd cpp_implementation
make all
```