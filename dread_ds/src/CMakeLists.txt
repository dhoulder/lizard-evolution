cmake_minimum_required (VERSION 3.0)
project(dreadds)
find_package (Boost 1.64
  COMPONENTS program_options filesystem REQUIRED)
find_package (GDAL REQUIRED)

include_directories(${Boost_INCLUDE_DIR} ${GDAL_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(dreadds STATIC
  model-config.cpp deme.cpp environment.cpp species.cpp
  output-file.cpp model.cpp simulation.cpp)
target_link_libraries(dreadds
  ${GDAL_LIBRARY}
  Boost::program_options
  Boost::filesystem)

add_executable(dreadds-bin dreadds.cpp)
set_target_properties(dreadds-bin PROPERTIES OUTPUT_NAME dreadds)
target_link_libraries(dreadds-bin
  dreadds
  ${GDAL_LIBRARY}
  Boost::program_options
  Boost::filesystem)

add_executable(grids2tsv grids2tsv.cpp)
target_link_libraries(grids2tsv
  dreadds
  ${GDAL_LIBRARY})

set_property(TARGET grids2tsv dreadds dreadds-bin PROPERTY CXX_STANDARD 11)
