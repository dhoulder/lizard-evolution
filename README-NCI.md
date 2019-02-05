# Instructions for building and running on nci.org.au systems

## On raijin.nci.org.au

### To build the dreadds library and executable

First, get the source code:
```
git clone git@github.com:DanRosauer/DREaD_ds.git
cd DREaD_ds/
git checkout master
```

To build using the Intel compiler. Use this recipe if you intend to
build the R package as well.

```
module purge
module load pbs # In case you want to submit jobs etc.
module load gcc/4.9.0
module load intel-cc/2018.3.222
module load gdal/2.2.2
module load boost/1.64.0
module load cmake/3.8.2
export BOOST_LIBRARYDIR=$BOOST_LIBRARYDIR/Intel
```

To build using gcc:
```
module purge # in case Intel compilers are loaded
module load pbs # In case you want to submit jobs etc.

module load gcc/4.9.0 # could probably use intel-cc instead
module load gdal/2.2.2
module load boost/1.64.0
module load cmake/3.8.2
export BOOST_LIBRARYDIR=$BOOST_LIBRARYDIR/GNU
```

To compile and link:
```
cd src

# If you've previously built with a different compiler (e.g. gcc vs
# intel-cc), first do:
#   rm CMakeCache.txt
# first
cmake .
cmake --build .
```

### To run the dreadds executable

```
./src/dreadds --help

# See ./README.md` for more information.
```

### To build and install R packages
```
# R was built with the Intel compilers, so we need to use those too

# If not already loaded, do...
module load gcc/4.9.0
module load intel-cc/2018.3.222
module load gdal/2.2.2
module load boost/1.64.0
export BOOST_LIBRARYDIR=$BOOST_LIBRARYDIR/Intel
# and then...
module load intel-fc/2018.3.222
module load R/3.5.1
module load proj/4.9.3
export PKG_CXXFLAGS=-std=c++11

R
install.packages('yaml')
# rgdal allows reading/writing common GIS file formats
install.packages('rgdal')
```

To use:

```
R
library(yaml)
data <- yaml.load("[some, yaml, as, a string]")
data <- yaml.load_file('test/1/out/out1-stats.yml')

library(rgdal)
# show supported file formats
gdalDrivers()
data <- readOGR('some-vector-data.ext')
```
