cmake_minimum_required (VERSION 3.8)
project(dreadds)
find_package (Boost 1.64
  COMPONENTS program_options filesystem REQUIRED)
find_package (GDAL REQUIRED)

include_directories(${Boost_INCLUDE_DIR} ${GDAL_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
enable_testing()
add_subdirectory(test)

add_library(dreadds STATIC
  model-config.cpp deme.cpp environment.cpp env-params.cpp
  output-file.cpp model.cpp species.cpp
)
target_link_libraries(dreadds
  ${GDAL_LIBRARY}
  Boost::program_options
  Boost::filesystem)

add_executable(dreadds-bin dreadds.cpp)
set_target_properties(dreadds-bin PROPERTIES OUTPUT_NAME dreadds)
target_link_libraries(dreadds-bin
  dreadds)

add_executable(grids2tsv grids2tsv.cpp)
target_link_libraries(grids2tsv
  dreadds)

add_executable(csv2grid csv2grid.cpp)
target_link_libraries(csv2grid
    ${GDAL_LIBRARY})

set_property(TARGET csv2grid grids2tsv dreadds dreadds-bin
             PROPERTY CXX_STANDARD 11)
