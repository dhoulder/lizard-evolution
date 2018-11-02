# DREaD_ds C++ implementation and R API

This directory contains C++ source code that implements the DREaD_ds
evolution model and an R package that allows it to be called from
R. See `../README.md`

The model implemented here is intended to closely match the model
implemented by `../DREaD_ds.R`, but in a language that allows for
greater computational performance and efficiency.

Source code for the C++ implementation of the model can be found in
`./src/`

Sample configuration and input files can be found in `./examples/`.

Files in `./r-package/` provide an R API that allows input data, model
configuration, output and visualisation to be handled in R.

This separation means that the model is only loosely coupled to R,
which opens up the possibility of using the model by, say, calling it
from a standalone executable that handles all the model inputs and
outputs.


## The model

See `./src/README.md` for details on building the library that
implements the model. This can be called from a standalone executable
if required, or called from R using the R package below.


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

### Building
```
cd r-package
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
> require(dreadds)
# run stub
> dreadds()
```
