# DREaD_ds C++ implementation

This directory contains C++ source code that implements the DREaD_ds
evolution model.

This builds the dread_ds static library (dread_ds.a on *nix).
The R package builds `dreadds.so`

### Prerequisites

* Eigen template library. See http://eigen.tuxfamily.org
  See `./CMakeLists.txt` for version information.

* `cmake`
  See https://cmake.org/

* A C++ compiler.
  Needs to support C++11

On Ubuntu these prerequisites can be installed with:

  `sudo apt-get install libeigen3-dev cmake build-essential`


### Building
```
cmake --build .
```

