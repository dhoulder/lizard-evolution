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
  # available. Check your local system documentation.
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

* ```
  # cd to the directory containing this README.md
  ```

* (Note: if you have changed the API in `rcpp_*.cpp`, do…
  ```
  # may need to clobber dreadds/R/RcppExports.R first
  mv -i dreadds/R/RcppExports.R /tmp/

  R
  ...
  > require(Rcpp)
  > compileAttributes("dreadds")
  ```
  …to update the package skeleton first. Note that those files are under
  version control)


* ```
  # Need to pass in the location of the library and header via environment. See Makevars
  DREAD_DS_SRC=$(cd ../src/ && pwd) R CMD build dreadds
  ```

### Installing
* ```
  # cd to the directory containing this README.md

  DREAD_DS_SRC=$(cd ../src && pwd) R
  > install.packages("dreadds_1.0.tar.gz", repos = NULL)
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
