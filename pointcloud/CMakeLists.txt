# Point-cloud pipeline (depth camera -> octree -> Barnes-Hut kernel sums)
option(BH_POINTCLOUD_PNG "Read 16-bit PNG depth images with libpng" ON)

add_library(bh_pointcloud STATIC
    point_cloud.cpp
    point_octree.cpp
    algorithms.cpp
    synthetic.cpp
    ../stdinc.cpp
)
target_include_directories(bh_pointcloud PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_compile_features(bh_pointcloud PUBLIC cxx_std_20)
if(OpenMP_CXX_FOUND)
    target_link_libraries(bh_pointcloud PUBLIC OpenMP::OpenMP_CXX)
endif()

if(BH_POINTCLOUD_PNG)
    find_package(PNG QUIET)
    if(PNG_FOUND)
        target_link_libraries(bh_pointcloud PUBLIC PNG::PNG)
        target_compile_definitions(bh_pointcloud PUBLIC BH_HAVE_PNG=1)
        message(STATUS "pointcloud: libpng found - PNG depth input enabled")
    else()
        message(STATUS "pointcloud: libpng not found - use PGM/raw depth or tools/depth_png_to_pgm.py")
    endif()
endif()

add_executable(pointcloud_tool pointcloud_tool.cpp)
target_link_libraries(pointcloud_tool PRIVATE bh_pointcloud)

add_executable(test_pointcloud tests/test_pointcloud.cpp)
target_link_libraries(test_pointcloud PRIVATE bh_pointcloud)
add_test(NAME pointcloud COMMAND test_pointcloud)

install(TARGETS pointcloud_tool RUNTIME DESTINATION bin)
