# Instructions for building and running on nci.org.au systems

## On raijin.nci.org.au (as of Nov. 2018)

### To build

```
git clone git@github.com:DanRosauer/DREaD_ds.git
cd DREaD_ds/
git checkout devel-cpp

module load gcc/4.9.0
module load gdal/2.2.2
module load boost/1.64.0
module load cmake/3.8.2
export BOOST_LIBRARYDIR=$BOOST_LIBRARYDIR/GNU  # FIXME? required?

cd dread_ds/src
cmake . # only required once to setup build files
cmake --build .
```

### To run

```
cd  ~/DREaD_ds/dread_ds/src
./test_dd $config_filename $n_iterations $output_path

# See ../examples/ for config format
```



