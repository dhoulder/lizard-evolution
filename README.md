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
environment. See https://www.gdal.org/formats_list.html for a list of
acceptable file formats.

Model parameters are usually specified by means of a configuration
file, but they can also be specified directly, either through
command-line arguments, or function arguments when using the R
API.

See below for details about the configuration options and file format.


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

### Configuration options and file format

The definitive list of configuration options and their meanings can be
obtained by running

```
src/dreadds -h
```

See `examples/example.conf` for an example configuration file and
comments that explain how the configuration options are used. In R,
see the package help pages for details. The R API provides a more
R-like syntax for many of the configuration options by means of
function arguments.

Options can be specified in any order. This includes those options
that are used to specify repeatable sets of values in a section, such
as `[env]` or `[species]`. For these cases, ordinal correspondence of
the repeatable options is used to form each set.

In other words, the first set is constructed from the first
occurrences of its constituent options, the second set from the second
occurrences and so on. In order to maintain readability, it's
probably best to group each set together.

For example, the following two snippets are equivalent, although the
first one is easier to read and maintain

```
[species]
# Recommended layout
# First species
niche-centre = 12, 34
niche-breadth = 1.2, 3.4
# other options for first species…

# Second species
niche-centre = 13, 35
niche-breadth = 1.3, 3.5
# other options for second species…

```

```
[species]
# This works but is not recommended
niche-centre = 12, 34 # for first species
niche-centre = 13, 35 # for second species
niche-breadth = 1.2, 3.4 # for first species
niche-breadth = 1.3, 3.5 # for second species
# remaining options for both species…
```

Options on the command line override those in the config file. For
options that can be repeated (for example `--species.niche-centre`),
the entire list of those options is overridden, effectively overriding
all of those options that were specified in the config file.

This config file is parsed by the Boost program_options library. See
https://www.boost.org/doc/libs/1_68_0/doc/html/program_options/overview.html#id-1.3.31.5.10.2

Comments are prefixed with a `#` character and can be inserted anywhere.

`key=value` specifies a value for an option, equivalent to the
command-line argument `--key=value`

`[something]` starts a section.

Below a `[section]`, `key=value` is equivalent to the command line
syntax `--section.key=value`.

Leading and trailing whitespace is stripped from both keys and values.
**Note that there is no quoting syntax, thus no way to specify a
literal `#` character, nor explicit leading or trailing spaces.** Also
note that shell metacharacters such as `~` and `$` are also
interpreted literally.
