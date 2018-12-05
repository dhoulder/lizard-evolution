# Dynamic Range Evolution and Diversification simulation with dynamic speciation (DREaD_ds)

This directory contains source code and documentation for the DREaD_ds
evolution model.

`./src` contains source code for the C++ implementation of the model.
The model implemented in C++ is intended to closely match the model
implemented by `./R-prototype/DREaD_ds.R`, but in a language that allows
for greater computational performance.

`./examples/` contains sample configuration and input files.

`./r-package/` contains files that provide an R API that allows the
model to be called from the R interpreter.

`./R-prototype/` contains background information and support files for
the original model prototype written in R.


## Building and running he model

`./src/README.md` contains instructions for building both the library
that implements the model, and a standalone executable that provides a
command-line interface to it.

To see the usage information for the standalone executable `dreadds`,
run it without arguments. `dreadds` uses a configuration file as part
of its input. See `examples/example.conf` for the configuration file
format and allowable options.

`./r-package/README.md` contains instructions for building the R API.
