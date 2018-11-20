# Instructions for building and running on nci.org.au systems

## On raijin.nci.org.au (as of Nov. 2018)

### To build


This recipe uses gcc. You could probably use the Intel C++ compiler
instead and it might produce faster code.

```
git clone git@github.com:DanRosauer/DREaD_ds.git
cd DREaD_ds/
git checkout devel-cpp

module purge # in case Intel compilers are loaded (conflict with gcc)
module load pbs # In case you want to submit jobs etc.

module load gcc/4.9.0 # could probably use intel-cc instead
module load gdal/2.2.2
module load boost/1.64.0
module load cmake/3.8.2
export BOOST_LIBRARYDIR=$BOOST_LIBRARYDIR/GNU # or $BOOST_LIBRARYDIR/Intel for intel-cc

cd dread_ds/src
cmake . # only required once to setup build files
cmake --build .
```

### To run

```
cd  ~/DREaD_ds/dread_ds/src
./dreadds --help

# See ../examples/ for config file format
```



