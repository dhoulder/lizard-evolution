# DREaD_ds C++ implementation and R API

This directory contains C++ source code that implements the DREaD_ds
evolution model and an R package that allows it to be called from
R. See `../README.md`

The model implemented here is intended to closely match the model
implemented by `../DREaD_ds.R`, but in a language that allows for
greater computational performance.

Source code for the C++ implementation of the model can be found in
`./src/`

Sample configuration and input files can be found in `./examples/`.

Files in `./r-package/` provide an R API that allows the model to be
called from the R interpreter.


## The model

`./src/README.md` contains instructions for building both the library
that implements the model, and a standalone executable that provides a
command line interface to it.

To see the usage information for the standalone executable `dreadds`,
run it without arguments. `dreadds` uses a configuration file as part
of its input. See `examples/example.conf` for the configuration file
format and allowable options.


## The R package (R API for the model)

### Prerequisites

* The DREaD_ds C++ library (`*.a`, `*.so` or `*.dll`) and
  header(s). See `./src/README.md` for build instructions.

* R and the Rcpp package. See
  http://dirk.eddelbuettel.com/code/rcpp.html

  The package build process also needs pdflatex.
  On Ubuntu 18.04 use::
  ```
  apt-get install r-base r-cran-rcpp texlive-latex-base
  ```

  On Windows, install Rtools from https://cran.r-project.org/bin/windows/Rtools/

### Building the R package
```
cd r-package
```

If you have changed the API in rcpp_*.cpp, use
```
# may need to clobber dreadds/R/RcppExports.R first
mv -i dreadds/R/RcppExports.R /tmp/

R
...
> require(Rcpp)
> compileAttributes("dreadds")
```
to update the package skeleton first.

```
# Need to pass in the location of the library and header via environment. See Makevars
DREAD_DS_SRC=$(cd ../src/ && pwd) R CMD build dreadds
```

### Installing
```
# cd to directory containing this README.md
DREAD_DS_SRC=$(cd src && pwd) R
> install.packages("r-package/dreadds_1.0.tar.gz", repos = NULL)
```

### Running
```
R
...
> require(dreadds)
# See…
> help(dreadds)
# for usage information
```
