# Instructions for building and running on nci.org.au systems

## On raijin.nci.org.au

### To build


This recipe uses gcc. You could probably use the Intel C++ compiler
instead and it might produce faster code.

```
git clone git@github.com:DanRosauer/DREaD_ds.git
cd DREaD_ds/
git checkout master

module purge # in case Intel compilers are loaded (conflict with gcc)
module load pbs # In case you want to submit jobs etc.

module load gcc/4.9.0 # could probably use intel-cc instead
module load gdal/2.2.2
module load boost/1.64.0
module load cmake/3.8.2
export BOOST_LIBRARYDIR=$BOOST_LIBRARYDIR/GNU # or $BOOST_LIBRARYDIR/Intel for intel-cc

cd src
cmake .
cmake --build .
```

### To run

```
./src/dreadds --help

# See ./README.md` for more information.
```
