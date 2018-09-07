cmake_minimum_required (VERSION 3.0)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_library(dread_ds dread_ds.cpp)
target_link_libraries (dread_ds Eigen3::Eigen)
