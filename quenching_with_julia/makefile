# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++17 -O0 -I/globalhome/tus210/HPC/include

# Linker flags
LDFLAGS = -L/globalhome/tus210/HPC/lib64

# Libraries
LIBS = -lcaf_core -lcaf_io

# # Source files
# SRCS = actor_julia_test_update.cc
# # Target executable
# TARGET = refinement_test


SRCS = actor_exhustive.cc
# Target executable
TARGET = distributed_test_update

# SRCS = baseline.cc
# # Target executable
# TARGET = baseline_test

# Set library path
export LD_LIBRARY_PATH := /globalhome/tus210/HPC/lib64:$(LD_LIBRARY_PATH)

# Set OpenMP environment variable
export OMP_NUM_THREADS := 1

# Build rule
all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(TARGET) $(SRCS) $(LIBS)

# Clean rule
clean:
	rm -f $(TARGET)
