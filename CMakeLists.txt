cmake_minimum_required(VERSION 3.20)
project(BarnesHutModern VERSION 2.0.0 LANGUAGES CXX)

# Modern C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags for optimization
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Optimization flags
if(MSVC)
    add_compile_options(/W4 /O2 /arch:AVX2)
    if(CMAKE_BUILD_TYPE MATCHES Release)
        add_compile_options(/GL /Oi /Ot)
        add_link_options(/LTCG)
    endif()
else()
    add_compile_options(-Wall -Wextra -pedantic)
    if(CMAKE_BUILD_TYPE MATCHES Release)
        add_compile_options(-O3 -march=native -mtune=native)
        # Enable more aggressive optimizations
        add_compile_options(-ffast-math -funroll-loops)
    endif()
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options(-g -O0)
    endif()
endif()

# Find OpenMP for parallel execution
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found - enabling parallel execution")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(WARNING "OpenMP not found - using serial execution")
endif()

# Source files
set(CORE_SOURCES
    stdinc.cpp
    particle.cpp
    tree.cpp
    file.cpp
)

set(CORE_HEADERS
    stdinc.h
    vektor.h
    particle.h
    tree.h
    file.h
)

# Main simulation executable
add_executable(barnes_hut_sim BHtreetest.cpp ${CORE_SOURCES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(barnes_hut_sim PRIVATE OpenMP::OpenMP_CXX)
endif()

# Data generator executable
add_executable(generate_data generate_data.cpp ${CORE_SOURCES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(generate_data PRIVATE OpenMP::OpenMP_CXX)
endif()

# Quadtree visualizer executable
add_executable(visualize_quadtree visualize_quadtree.cpp quadtree_visualizer.cpp ${CORE_SOURCES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(visualize_quadtree PRIVATE OpenMP::OpenMP_CXX)
endif()

# Installation
install(TARGETS barnes_hut_sim generate_data visualize_quadtree
        RUNTIME DESTINATION bin)

# Print configuration summary
message(STATUS "")
message(STATUS "=== Barnes-Hut Configuration ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "OpenMP: ${OpenMP_CXX_FOUND}")
message(STATUS "=================================")
message(STATUS "")
