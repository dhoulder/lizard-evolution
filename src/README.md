# DREaD_ds C++ implementation

This directory contains C++ source code that implements the DREaD_ds
evolution model and provides an API in the form of a runtime library.

The API is defined in ./model.h

### Development notes

* Remember to to keep the version string in `constants.h` up to date


### Prerequisites

* C++11 compiler. For gcc, version 4.8 or later

* Boost C++ library. See https://www.boost.org/
  See `./CMakeLists.txt` for version information.

* `cmake`
  See https://cmake.org/

* `GDAL`
  See  https://www.gdal.org/

On Ubuntu these prerequisites can be installed with:

  ```
  sudo apt-get install \
       libboost-dev libboost-program-options-dev libboost-filesystem-dev \
       libgdal-dev cmake build-essential
  ```


### Building
```
cmake . # first time only to initialise build files
cmake --build .
make test # Run test suite
```
