# DREaD_ds C++ implementation

This directory contains C++ source code that implements the DREaD_ds
evolution model and provides an API in the form of a runtime library.

See ./simulation.h for API usage.

### Prerequisites

* C++11 compiler. For gcc, version 4.8 or later

* Boost C++ library. See https://www.boost.org/
  See `./CMakeLists.txt` for version information.

* `cmake`
  See https://cmake.org/

* `GDAL`
  See  https://www.gdal.org/

On Ubuntu these prerequisites can be installed with:

  `sudo apt-get install libboost-dev libboost-program-options-dev libgdal-dev cmake build-essential`


### Building
```
cmake . # first time only to initialise build files
cmake --build .
```

