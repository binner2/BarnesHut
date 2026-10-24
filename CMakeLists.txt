cmake_minimum_required(VERSION 3.24)
project(BarnesHutModern VERSION 2.1.0 LANGUAGES CXX)

option(BH_ENABLE_CUDA "Build the GPU Barnes-Hut library, tests and benchmark" ON)
option(BH_FAST_MATH "Enable -ffast-math / -use_fast_math (breaks NaN/Inf checks in tests)" OFF)
option(BH_CUDA_SANITIZE "Instrument device code for compute-sanitizer memcheck (-fdevice-sanitize=memcheck)" OFF)

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
    # Host-compiler flags only: nvcc does not accept them directly
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wall> $<$<COMPILE_LANGUAGE:CXX>:-Wextra> $<$<COMPILE_LANGUAGE:CXX>:-pedantic>)
    if(CMAKE_BUILD_TYPE MATCHES Release)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-O3> $<$<COMPILE_LANGUAGE:CXX>:-march=native> $<$<COMPILE_LANGUAGE:CXX>:-funroll-loops>)
        if(BH_FAST_MATH)
            add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-ffast-math>)
        endif()
    endif()
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-g> $<$<COMPILE_LANGUAGE:CXX>:-O0>)
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

# CUDA configuration. The device code needs only the CUDA 13.x baseline
# (sm_75+); CUDA 13.4 / CCCL 3.4 is the tested toolchain. CUDA stays
# optional so the CPU code builds on machines without a toolkit.
set(BH_HAVE_CUDA OFF)
if(BH_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        # Default fatbin: SASS for Turing..Blackwell (incl. GB10 = sm_121)
        # plus compute_121 PTX for forward compatibility. Override with
        # -DCMAKE_CUDA_ARCHITECTURES=native (or e.g. "89;121").
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES 75-real 80-real 86-real 89-real 90-real 100-real 120-real 121-real 121-virtual)
        endif()
        enable_language(CUDA)
        find_package(CUDAToolkit 13.0 REQUIRED)
        if(CUDAToolkit_VERSION VERSION_LESS 13.4)
            message(WARNING "CUDA ${CUDAToolkit_VERSION} found; this code is tested with CUDA 13.4 (CCCL 3.4)")
        endif()
        set(CMAKE_CUDA_STANDARD 20)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        set(BH_HAVE_CUDA ON)

        set(BH_CUDA_OPTIONS -lineinfo)
        if(BH_FAST_MATH)
            list(APPEND BH_CUDA_OPTIONS -use_fast_math)
        endif()
        if(BH_CUDA_SANITIZE)
            list(APPEND BH_CUDA_OPTIONS -fdevice-sanitize=memcheck)
        endif()
    else()
        message(WARNING "No CUDA compiler found - building CPU targets only")
    endif()
endif()

# Find OpenGL and related libraries for visualization
find_package(OpenGL QUIET)
find_package(glfw3 QUIET)
find_package(glm QUIET)

if(BH_HAVE_CUDA AND OpenGL_FOUND AND glfw3_FOUND AND glm_FOUND)
    set(ENABLE_VISUALIZATION ON)
    message(STATUS "OpenGL, GLFW, and GLM found - enabling visualization")
else()
    set(ENABLE_VISUALIZATION OFF)
    message(WARNING "Visualization dependencies not found - building without visualization")
    if(NOT OpenGL_FOUND)
        message(WARNING "  OpenGL not found")
    endif()
    if(NOT glfw3_FOUND)
        message(WARNING "  GLFW3 not found - install with: sudo apt-get install libglfw3-dev")
    endif()
    if(NOT glm_FOUND)
        message(WARNING "  GLM not found - install with: sudo apt-get install libglm-dev")
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

# GPU Barnes-Hut library, host-emulation tests and benchmark
if(BH_HAVE_CUDA)
    add_library(bh_gpu STATIC cuda/bh_gpu.cu)
    target_include_directories(bh_gpu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/cuda)
    target_link_libraries(bh_gpu PUBLIC CUDA::cudart)
    target_compile_options(bh_gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${BH_CUDA_OPTIONS}>)

    # Runs the device algorithm on the host: no GPU needed
    add_executable(bh_host_core_test cuda/tests/host_core_test.cu)
    target_include_directories(bh_host_core_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda)

    add_executable(benchmark_cuda
        cuda/test_scenarios.cpp
        cuda/benchmark_framework.cu
        ${CORE_SOURCES}
    )
    target_include_directories(benchmark_cuda PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/cuda
    )
    target_compile_options(benchmark_cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${BH_CUDA_OPTIONS}>)
    target_link_libraries(benchmark_cuda PRIVATE bh_gpu)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(benchmark_cuda PRIVATE OpenMP::OpenMP_CXX)
        target_compile_options(benchmark_cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>)
    endif()

    install(TARGETS benchmark_cuda RUNTIME DESTINATION bin)

    enable_testing()
    add_test(NAME host_core COMMAND bh_host_core_test)
    # CPU tests always run; GPU tests are skipped when no device is present
    add_test(NAME correctness COMMAND benchmark_cuda --correctness-only)
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
message(STATUS "CUDA: ${BH_HAVE_CUDA} (${CUDAToolkit_VERSION})")
message(STATUS "CUDA Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "Visualization: ${ENABLE_VISUALIZATION}")
message(STATUS "GPU benchmark/tests: ${BH_HAVE_CUDA} (benchmark_cuda, bh_host_core_test)")
message(STATUS "=================================")
message(STATUS "")
