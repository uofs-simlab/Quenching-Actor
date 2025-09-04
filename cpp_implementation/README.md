# Quenching Actor System 

This directory contains the C++ implementation of the quenching study system using the CAF (C++ Actor Framework).

## Architecture Overview

### Core Components

1. **Dynamic Version** (`actor_cpp_dynamic.cc`)
   - Performs adaptive refinement immediately when boundaries are detected
   - Uses spatial hashing for efficient neighbor management
   - Supports both early termination and bracket optimization

2. **Static Version** (`actor_cpp_static.cc`)  
   - Completes all jobs at each refinement level before proceeding to the next
   - Uses synchronized refinement across all gamma (gg) levels
   - More systematic exploration of parameter space
   - Supports both early termination and bracket optimization

### Key Features

#### Bracket Optimization
- **Shared Implementation**: Both versions use the same bracket optimizer module
- **Strategy 1**: Uses completed neighbor results to initialize solver brackets
- **Strategy 2**: Uses cache from higher gamma levels for the same parameters

#### Early Termination
- **Sundials Integration**: Stops solver early when convergence criteria are met
- **Performance Boost**: Significantly reduces computation time 
- **Available in Both Versions**: Can be enabled/disabled independently

## File Structure

```
cpp_implementation/
├── include/
│   ├── bracket_optimizer.h          # Bracket optimization interface
│   ├── dynamic_neighbor_provider.h  # Neighbor provider for dynamic system
│   ├── static_neighbor_provider.h   # Neighbor provider for static system
│   ├── job_structures.h             # Common data structures
│   ├── system_config.h              # Configuration management
│   ├── config.h                     # CAF configuration
│   ├── fitzhugh_nagumo.h           # Physics model
│   ├── bisection_solver.h          # Numerical solver
│   └── wave_loader.h               # Wave data loading
├── src/
│   ├── actor_cpp_dynamic.cc        # Dynamic adaptive refinement system
│   ├── actor_cpp_static.cc         # Static level-by-level refinement system  
│   ├── bracket_optimizer.cpp       # Bracket optimization implementation
│   ├── fitzhugh_nagumo.cpp         # Physics model implementation
│   ├── bisection_solver.cpp        # Numerical solver implementation
│   └── wave_loader.cpp             # Wave data loading implementation
└── Makefile                        # Build configuration
```

## Building

```bash
# Build both versions
make all

# Build individual versions
make dynamic
make static

# Clean build artifacts
make clean
```

### Configuration Options

- `--server-mode`: Run as computation server
- `--host HOST`: Connect to server at HOST (client mode)
- `--port PORT`: Server port (default: 8080)
- `--enable-bracket`: Enable bracket optimization 
- `--enable-early-termination`: Enable early solver termination 