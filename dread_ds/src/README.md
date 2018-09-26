# DREaD_ds C++ implementation

This directory contains C++ source code that implements the DREaD_ds
evolution model and provides an API in the form of a runtime library.

See ./dread_ds.h for API usage.

### Prerequisites

* Boost C++ library. See https://www.boost.org/
  See `./CMakeLists.txt` for version information.

* `cmake`
  See https://cmake.org/

* A C++ compiler.
  Needs to support C++11

On Ubuntu these prerequisites can be installed with:

  `sudo apt-get install libboost-dev cmake build-essential`


### Building
```
cmake --build .
```

