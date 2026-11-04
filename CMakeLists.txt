cmake_minimum_required(VERSION 3.20)
project(BarnesHutModern VERSION 2.1.0 LANGUAGES CXX)

# Optional components. The core simulation, the point-cloud pipeline and the
# tests build with a plain C++20 compiler; CUDA/OpenGL are only needed for the
# interactive visualiser.
option(BH_ENABLE_VISUALIZATION "Build the CUDA/OpenGL visualiser (needs CUDA toolkit, GLFW, GLM)" OFF)
option(BH_BUILD_POINTCLOUD "Build the depth-camera point-cloud pipeline" ON)
option(BH_BUILD_TESTS "Build and register tests" ON)

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

# CUDA + OpenGL configuration (visualiser only)
set(ENABLE_VISUALIZATION OFF)
if(BH_ENABLE_VISUALIZATION)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_ARCHITECTURES 75 80 86 89)  # Turing, Ampere, Ada Lovelace
    if(CMAKE_BUILD_TYPE MATCHES Release)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -use_fast_math --expt-relaxed-constexpr")
    else()
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G -g -O0")
    endif()

    find_package(OpenGL REQUIRED)
    find_package(glfw3 QUIET)
    find_package(glm QUIET)
    if(OpenGL_FOUND AND glfw3_FOUND AND glm_FOUND)
        set(ENABLE_VISUALIZATION ON)
        message(STATUS "OpenGL, GLFW, and GLM found - enabling visualization")
    else()
        message(WARNING "Visualization dependencies not found - building without visualization")
        if(NOT glfw3_FOUND)
            message(WARNING "  GLFW3 not found - install with: sudo apt-get install libglfw3-dev")
        endif()
        if(NOT glm_FOUND)
            message(WARNING "  GLM not found - install with: sudo apt-get install libglm-dev")
        endif()
    endif()
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

# Depth-camera point-cloud pipeline
if(BH_BUILD_TESTS)
    enable_testing()
endif()
if(BH_BUILD_POINTCLOUD)
    add_subdirectory(pointcloud)
endif()

# Visualization executable with CUDA and OpenGL
if(ENABLE_VISUALIZATION)
    # GLAD for OpenGL function loading (we'll include it in visualization/)
    set(GLAD_SOURCES visualization/glad/glad.c)

    set(VISUALIZATION_SOURCES
        visualization/cuda_renderer.cu
        visualization/gl_renderer.cpp
        visualization/window_manager.cpp
        visualization/camera.cpp
        ${GLAD_SOURCES}
    )

    set(VISUALIZATION_HEADERS
        visualization/cuda_renderer.cuh
        visualization/gl_renderer.h
        visualization/window_manager.h
        visualization/camera.h
    )

    add_executable(barnes_hut_visual
        visualization/visualization_main.cpp
        ${CORE_SOURCES}
        ${VISUALIZATION_SOURCES}
    )

    target_include_directories(barnes_hut_visual PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/visualization
        ${CMAKE_CURRENT_SOURCE_DIR}/visualization/glad
    )

    target_link_libraries(barnes_hut_visual PRIVATE
        OpenGL::GL
        glfw
        CUDA::cudart
        CUDA::cuda_driver
    )

    if(OpenMP_CXX_FOUND)
        target_link_libraries(barnes_hut_visual PRIVATE OpenMP::OpenMP_CXX)
    endif()

    # Set CUDA properties
    set_target_properties(barnes_hut_visual PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )

    install(TARGETS barnes_hut_visual RUNTIME DESTINATION bin)

    # Copy shaders to build directory
    file(COPY visualization/shaders DESTINATION ${CMAKE_BINARY_DIR})
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "=== Barnes-Hut Configuration ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "OpenMP: ${OpenMP_CXX_FOUND}")
message(STATUS "Visualization (CUDA/OpenGL): ${ENABLE_VISUALIZATION}")
message(STATUS "Point-cloud pipeline: ${BH_BUILD_POINTCLOUD}")
message(STATUS "=================================")
message(STATUS "")
