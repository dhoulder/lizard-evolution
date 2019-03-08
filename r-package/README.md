## The R package (R API for the model)

### Prerequisites

* The DREaD_ds C++ library (`*.a`, `*.so` or `*.dll`) and
  header(s). See `../src/README.md` for build instructions.

* R and the Rcpp package. See
  http://dirk.eddelbuettel.com/code/rcpp.html

  Rcpp can be installed for an individual user with
  ```
  # You may need to set your environment first with "module load ..."
  # or similar so that the correct version of R and the C++ compiler is
  # available. See ../README-*.md and your local system documentation.
  R
  install.packages('Rcpp')
  ```

  It can also be installed system-wide

  * On Ubuntu 18.04 use::
    ```
    apt-get install r-base r-cran-rcpp texlive-latex-base
    ```
    On RHEL, Fedora and CentOS, install the R-Rcpp package. You may need
    to enable the EPEL package repository first. See
    https://fedoraproject.org/wiki/EPEL

  * On Windows, install Rtools from
    https://cran.r-project.org/bin/windows/Rtools/

### Building the R package

* Make sure you've set up your compilation environment first.
  See `../README-*.md` and `../README.md`.
  The build process needs to be able to find the libraries that the DREaD_ds
  library depends on.

* Change your working directory to the one containing this `README.md`
  if you're not already here. e.g.
  ```
  cd r-package/
  ```

* The build process needs to know where to find the DREaD_ds C++
  library and header files (see `dreadds/src/Makevars`). This is done
  by setting the `DREAD_DS_SRC` environment variable. Assuming your
  working directory is `r-package` and you're using the `/bin/bash`
  shell, use…
  ```
  export DREAD_DS_SRC=$(cd ../src/ && pwd)
  export PKG_CXXFLAGS=-std=c++11
  ```
  On Windows, use the `SET` command to set `DREAD_DS_SRC` to the absolute path
  of `../src`

* Use
  ```
  R CMD build dreadds
  ```
  to create the package file (`dreadds*.tar.gz`). Note that this doesn't compile
  any C++ code at this stage.

* Install the package with
  ```
  R
  > install.packages("dreadds_1.0.tar.gz", repos = NULL)
  ```
  Note that `install.packages("dreadds", repos = NULL)` can also be used to compile and
  install the package without having to wrap everything up in an archive first.

### Running
```
R
...
> require(dreadds)
# See…
> help(dreadds)
# for usage information
```


### Modifying the R package

* If you have changed the API in `rcpp_*.cpp`, do…
  ```
  mv -i dreadds/src/RcppExports.cpp ~/

  R
  ...
  > require(Rcpp)
  > compileAttributes("dreadds")
  > q()
  ```
  …to update the package skeleton first. Note that those files are under
  version control
