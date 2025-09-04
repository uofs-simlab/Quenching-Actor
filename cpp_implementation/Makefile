# Makefile for Quenching Actor System - Updated for HPC Environment
# Builds both dynamic and static versions with shared bracket optimizer

CXX = g++

# Use HPC environment modules and system paths (like your existing Makefiles)
CXXFLAGS = -std=c++17 -O2 -I/globalhome/tus210/HPC/include -I$(EBROOTSUNDIALS)/include -I$(EBROOTGSL)/include -I$(EBROOTEIGEN)/include -Iinclude

# Linker flags using HPC environment
LDFLAGS = -L/globalhome/tus210/HPC/lib64 -L$(EBROOTSUNDIALS)/lib -L$(EBROOTGSL)/lib

# Libraries including scientific computing libraries
LIBS = -lcaf_core -lcaf_io -lsundials_cvode -lsundials_nvecserial -lgsl -lgslcblas -lm

# Set library path for runtime
export LD_LIBRARY_PATH := /globalhome/tus210/HPC/lib64:$(LD_LIBRARY_PATH)

# Set OpenMP environment variable
export OMP_NUM_THREADS := 1

# Source directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

# Common source files (includes new bracket_optimizer and your existing modules)
COMMON_SOURCES = $(SRC_DIR)/bracket_optimizer.cpp \
                 $(SRC_DIR)/fitzhugh_nagumo.cpp \
                 $(SRC_DIR)/bisection_solver.cpp \
                 $(SRC_DIR)/wave_loader.cpp \
                 $(SRC_DIR)/finite_difference.cpp \
                 $(SRC_DIR)/wave_analysis.cpp

# Executable targets
DYNAMIC_TARGET = $(BUILD_DIR)/actor_cpp_dynamic
STATIC_TARGET = $(BUILD_DIR)/actor_cpp_static

# Object files
COMMON_OBJECTS = $(COMMON_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
DYNAMIC_OBJECT = $(BUILD_DIR)/actor_cpp_dynamic.o
STATIC_OBJECT = $(BUILD_DIR)/actor_cpp_static.o

.PHONY: all clean dynamic static install legacy_dynamic legacy_static

all: dynamic static

dynamic: $(DYNAMIC_TARGET)

static: $(STATIC_TARGET)

# Legacy targets that match your existing Makefile naming
legacy_dynamic: actor_cpp_bisection
legacy_static: actor_cpp_static_legacy

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Common object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Dynamic version (new modular)
$(DYNAMIC_TARGET): $(DYNAMIC_OBJECT) $(COMMON_OBJECTS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

# Static version (new modular)
$(STATIC_TARGET): $(STATIC_OBJECT) $(COMMON_OBJECTS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

# Individual object files
$(DYNAMIC_OBJECT): $(SRC_DIR)/actor_cpp_dynamic.cc | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(STATIC_OBJECT): $(SRC_DIR)/actor_cpp_static.cc | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Legacy targets for compatibility with your existing workflow
actor_cpp_bisection: $(SRC_DIR)/actor_cpp_dynamic.cc $(COMMON_SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

actor_cpp_static_legacy: $(SRC_DIR)/actor_cpp_static.cc $(COMMON_SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

# Installation (copy to bin directory if it exists)
install: all
	@if [ -d "bin" ]; then \
		cp $(DYNAMIC_TARGET) bin/; \
		cp $(STATIC_TARGET) bin/; \
		echo "Installed executables to bin/"; \
	else \
		echo "No bin/ directory found, skipping install"; \
	fi

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)
	rm -f actor_cpp_bisection actor_cpp_static_legacy

# Debug build
debug: CXXFLAGS := $(filter-out -O2,$(CXXFLAGS)) -g -O0 -DDEBUG
debug: all

# Show build configuration
info:
	@echo "Compiler: $(CXX)"
	@echo "Flags: $(CXXFLAGS)"
	@echo "Libraries: $(LIBS)"
	@echo "Dynamic target: $(DYNAMIC_TARGET)"
	@echo "Static target: $(STATIC_TARGET)"
	@echo "LD_LIBRARY_PATH: $(LD_LIBRARY_PATH)"

# Dependencies (simplified)
$(DYNAMIC_OBJECT): $(INCLUDE_DIR)/job_structures.h $(INCLUDE_DIR)/bracket_optimizer.h $(INCLUDE_DIR)/dynamic_neighbor_provider.h
$(STATIC_OBJECT): $(INCLUDE_DIR)/job_structures.h $(INCLUDE_DIR)/bracket_optimizer.h $(INCLUDE_DIR)/static_neighbor_provider.h $(INCLUDE_DIR)/system_config.h
$(BUILD_DIR)/bracket_optimizer.o: $(INCLUDE_DIR)/bracket_optimizer.h $(INCLUDE_DIR)/job_structures.h
