# Dynamic Range Evolution and Diversification simulation with dynamic speciation (DREaD_ds)

This directory contains source code and documentation for the DREaD_ds
evolution model. It implements a spatially explicit dynamic
intraspecific macroevolutionary model to simulate expected spatial
patterns of species and phylogenetic diversity. The primary purpose is
to select appropriate models which generate expected spatial patterns
best matching observed diversity. It is used to determine the processes
responsible for the spatial distribution of biological diversity and
phylogenetic endemism.

* See `./README-*.md` for site- and platform-specific information.

* `./src` contains source code for the C++ implementation of the model
  and several utility programs. The model implemented in C++ is
  intended to closely match the model implemented by
  `./R-prototype/DREaD_ds.R`, but in a language that allows for
  greater computational performance.

* `./examples/` contains sample configuration and input files.

* `./r-package/` contains files that provide an API for the R
  statistical computing language that allows the model to be run within
  the R interpreter.

* `./scripts/` contains shell scripts to assist in running the model
  and pre-processing input data. Run these without arguments for usage
  information.

* `./R-prototype/` contains background information and support files
  for the original model prototype written in R.


## Building and running he model

The model can be run using either the `dreadds` standalone executable,
or the R API.

Both the command-line interface and the R-language interface require
input files of gridded environment values to specify the simulation
environment.

Model parameters are usually specified by means of a configuration
file, but they can also be specified directly, either through
command-line arguments, or function arguments when using the R
API. See `examples/example.conf` for the configuration file format and
comments that explain how the configuration options are used. In R,
see the package help pages for details.

* `./src/README.md` contains instructions for building both the library
  that implements the model, and the command-line interface.

  To see the usage information and a description of all configuration
  options for the `dreadds` executable, run
  ```
  src/dreadds -h
  ```

  The build procedure will also generate `src/grids2tsv` which can be
  used to convert rasters of climate values to vector form.

* `./r-package/README.md` contains instructions for building and
  using the R API.
