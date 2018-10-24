cmake_minimum_required (VERSION 3.0)
project(dread_ds)
find_package (Boost 1.65 COMPONENTS program_options REQUIRED)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(dread_ds STATIC model-config.cpp load-env.cpp species.cpp model.cpp simulation.cpp)
target_link_libraries(dread_ds ${Boost_LIBRARIES})

add_executable(test_dd test.cpp)
target_link_libraries(test_dd dread_ds)
