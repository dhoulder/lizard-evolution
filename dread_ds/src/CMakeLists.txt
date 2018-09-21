cmake_minimum_required (VERSION 3.0)
project(dread_ds)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(dread_ds dread_ds.cpp)
target_link_libraries(dread_ds Eigen3::Eigen)

add_executable(test_dd test.cpp)
target_link_libraries(test_dd dread_ds)
