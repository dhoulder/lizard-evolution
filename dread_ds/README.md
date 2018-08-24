# DREaD_ds C++ implementation

This directory contains C++ source code that implements the DREaD_ds
evolution model. See `../README.md`

The model implemented here is intended to closely match the model
implemented by `../DREaD_ds.R`, but in a language that allows for
greater computational performance and efficiency.

Source code for the model can be found in `./src/`

Files in `./r-package/` provide an R API that allows input data, model
configuration, output and visualisation to be handled in R.

This separation means that the model is only loosely coupled to R,
which opens up the possibility of using the model by, say, calling it
from a standalone executable that handles all the model inputs and
outputs.

## The R package

### Prerequisites

Install R and the Rcpp package. See
http://dirk.eddelbuettel.com/code/rcpp.html

The package build process also needs pdflatex.
On Ubuntu 18.04 use::
```
apt-get install r-base r-cran-rcpp texlive-latex-base
```

### Building
```
cd r-package
R CMD check dreadds # may report errors generating PDF
R CMD build dreadds
```

### Installing
```
# cd to directory containing this README.md
R
> install.packages("r-package/dreadds_1.0.tar.gz", repos = NULL)
```

### Running
```
R
> require(dreadds)
# run stub
> dreadds()
```
