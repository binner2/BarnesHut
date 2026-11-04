# Modern Makefile for Barnes-Hut Simulation
# Alternative to CMake for quick builds

CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -pedantic
OPTFLAGS := -O3 -march=native -mtune=native -ffast-math -funroll-loops
DEBUGFLAGS := -g -O0 -DDEBUG

# Check for OpenMP support
OPENMP := $(shell $(CXX) -fopenmp -E - </dev/null >/dev/null 2>&1 && echo "-fopenmp" || echo "")
ifneq ($(OPENMP),)
    $(info OpenMP support detected)
    CXXFLAGS += $(OPENMP)
else
    $(warning OpenMP not found - building without parallelization)
endif

# Source files
CORE_SOURCES := stdinc.cpp particle.cpp tree.cpp file.cpp
CORE_OBJECTS := $(CORE_SOURCES:.cpp=.o)

# Point-cloud pipeline (pointcloud/)
PC_DIR := pointcloud
PC_SOURCES := $(PC_DIR)/point_cloud.cpp $(PC_DIR)/point_octree.cpp $(PC_DIR)/algorithms.cpp $(PC_DIR)/synthetic.cpp
PC_OBJECTS := $(PC_SOURCES:.cpp=.o)
PNG_CFLAGS := $(shell pkg-config --cflags libpng 2>/dev/null)
PNG_LIBS := $(shell pkg-config --libs libpng 2>/dev/null)
ifneq ($(PNG_LIBS),)
    $(info libpng detected - PNG depth input enabled)
    CXXFLAGS += -DBH_HAVE_PNG=1 $(PNG_CFLAGS)
endif

# Targets
TARGETS := barnes_hut_sim generate_data pointcloud_tool test_pointcloud

.PHONY: all clean release debug install test

# Default: Release build
all: release

# Release build
release: CXXFLAGS += $(OPTFLAGS)
release: $(TARGETS)

# Debug build
debug: CXXFLAGS += $(DEBUGFLAGS)
debug: $(TARGETS)

# Main simulation executable
barnes_hut_sim: BHtreetest.o $(CORE_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	@echo "Built: $@"

# Data generator executable
generate_data: generate_data.o $(CORE_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	@echo "Built: $@"

# Point-cloud tool
pointcloud_tool: $(PC_DIR)/pointcloud_tool.o $(PC_OBJECTS) stdinc.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(PNG_LIBS)
	@echo "Built: $@"

# Point-cloud tests
test_pointcloud: $(PC_DIR)/tests/test_pointcloud.o $(PC_OBJECTS) stdinc.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(PNG_LIBS)
	@echo "Built: $@"

test: test_pointcloud
	./test_pointcloud

# Object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Dependencies (generated automatically)
BHtreetest.o: BHtreetest.cpp file.h tree.h particle.h vektor.h stdinc.h
generate_data.o: generate_data.cpp file.h particle.h vektor.h stdinc.h
stdinc.o: stdinc.cpp stdinc.h
particle.o: particle.cpp particle.h vektor.h stdinc.h
tree.o: tree.cpp tree.h particle.h vektor.h stdinc.h
file.o: file.cpp file.h particle.h vektor.h stdinc.h

# Installation
install: release
	@mkdir -p $(HOME)/bin
	cp $(TARGETS) $(HOME)/bin/
	@echo "Installed to $(HOME)/bin/"

# Clean
clean:
	rm -f *.o $(PC_DIR)/*.o $(PC_DIR)/tests/*.o $(TARGETS)
	@echo "Cleaned build artifacts"

# Help
help:
	@echo "Barnes-Hut Simulation Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build release version (optimized)"
	@echo "  make debug        - Build debug version"
	@echo "  make clean        - Remove build artifacts"
	@echo "  make install      - Install to ~/bin"
	@echo "  make test         - Build and run point-cloud tests"
	@echo ""
	@echo "OpenMP: $(if $(OPENMP),Enabled,Disabled)"
	@echo "Compiler: $(CXX)"
