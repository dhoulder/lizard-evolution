cmake_minimum_required (VERSION 3.0)
project(dread_ds)
find_package (Boost 1.65 COMPONENTS program_options REQUIRED)
find_package (GDAL REQUIRED)

include_directories(${Boost_INCLUDE_DIR} ${GDAL_INCLUDE_DIR})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(dread_ds STATIC
  model-config.cpp deme.cpp environment.cpp species.cpp
  model.cpp simulation.cpp)
target_link_libraries(dread_ds ${Boost_LIBRARIES} ${GDAL_LIBRARY})

add_executable(test_dd test.cpp)
target_link_libraries(test_dd dread_ds)
