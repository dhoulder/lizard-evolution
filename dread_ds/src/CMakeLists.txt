cmake_minimum_required (VERSION 3.0)
project(dread_ds)
find_package (Boost 1.64 COMPONENTS program_options REQUIRED)
find_package (GDAL REQUIRED)

include_directories(${Boost_INCLUDE_DIR} ${GDAL_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(dread_ds STATIC
  model-config.cpp deme.cpp environment.cpp species.cpp
  model.cpp simulation.cpp)
target_link_libraries(dread_ds ${GDAL_LIBRARY} Boost::program_options)

add_executable(test_dd test.cpp)
target_link_libraries(test_dd
	dread_ds
	${GDAL_LIBRARY}
	Boost::program_options)
set_property(TARGET dread_ds test_dd PROPERTY CXX_STANDARD 11)
