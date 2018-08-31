# DREaD_ds C++ implementation

This directory contains C++ source code that implements the DREaD_ds
evolution model.

This builds the dread_ds static library (dread_ds.a on *nix).
The R package builds `dreadds.so`

### Prerequisites
* `cmake`

  On Ubuntu:
  `sudo apt-get install cmake`

* C++ compiler

  Needs to support C++11
  On Ubuntu:
  `sudo apt-get install build-essential`

### Building
```
cmake --build .
```

