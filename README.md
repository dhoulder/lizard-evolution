# Dynamic Range Evolution and Diversification simulation with dynamic speciation (DREaD_ds)

This directory contains source code and documentation for the DREaD_ds
evolution model. See `README-*.md` for site- and platform-specific
information.

`./src` contains source code for the C++ implementation of the model
and several utility programs. The model implemented in C++ is
intended to closely match the model implemented by
`./R-prototype/DREaD_ds.R`, but in a language that allows for greater
computational performance.

`./examples/` contains sample configuration and input files.

`./scripts/` contains shell scripts to assist in running the model and
pre-processing input data. Run these without arguments for usage
information.

`./r-package/` contains files that provide an R API that allows the
model to be called from the R interpreter.

`./R-prototype/` contains background information and support files for
the original model prototype written in R.


## Building and running he model

`./src/README.md` contains instructions for building both the library
that implements the model, and a standalone executable (`dreadds`)
that provides a command-line interface to it.

To see the usage information and a description of all configuration
options for the `dreadds` executable, run
```
src/dreadds -h
```

Configuration options are usually supplied to `dreadds` by means of a
a configuration file. Command-line options override options in the
configuration file. See `examples/example.conf` for the configuration
file format and comments that explain how the configuration options
are used.

The build procedure will also generate `src/grids2tsv` which can be
used to convert rasters of climate values to vector form.


`./r-package/README.md` contains instructions for building the R API.
