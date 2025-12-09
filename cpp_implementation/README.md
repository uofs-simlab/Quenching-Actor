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

3. **Shared Components**
   - Both versions use the same bracket optimization and solver components
   - Common header files define job structures and system configuration
   - Unified physics model and numerical solver implementations

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
│   ├── box_neighbor_provider.h      # Spatial neighbor detection for both versions
│   ├── job_structures.h             # Common data structures
│   ├── system_config.h              # Configuration management (L=2700 based)
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
├── Makefile.cpp_dynamic            # Build configuration for dynamic version
├── Makefile.cpp_static             # Build configuration for static version  
└── Makefile                        # Unified build configuration (optional)
```

## Building

The project uses separate Makefiles for different versions:

```bash
# Build dynamic version
make -f Makefile.cpp_dynamic clean all

# Build static version  
make -f Makefile.cpp_static clean
make -f Makefile.cpp_static

# Clean build artifacts
make -f Makefile.cpp_dynamic clean
make -f Makefile.cpp_static clean
```

## Running with SLURM

Use the provided SLURM scripts for HPC execution:

```bash
# Submit dynamic version job
sbatch cpp_quench_dynamic.sh

# Submit static version job
sbatch cpp_quench_static.sh
```

### Configuration Options

- `--server-mode`: Run as computation server
- `--host HOST`: Connect to server at HOST (client mode)  
- `--port PORT`: Server port (default: 8080)
- `--enable-bracket`: Enable bracket optimization 
- `--enable-early-termination`: Enable early solver termination

### Script Configuration

The SLURM scripts handle:
- Environment module loading (gcc, openmpi, eigen, gsl, sundials)
- Library path configuration
- Multi-node execution with server/client pattern

#### Dynamic Version (`cpp_quench_dynamic.sh`):
- Builds executable: `actor_cpp_dynamic` 
- Uses early termination and bracket optimization

#### Static Version (`cpp_quench_static.sh`):
- Builds executable: `actor_cpp_static`