# This makefile just provides some convenience targets for development.
# The C++ API is built by cmake. See ./src/README.md
# For details on building and using the R package, see ./README.md

SHELL=/bin/bash

R_PKG_TGZ=dreadds_1.0.tar.gz # Must match version in r-package/dreadds/DESCRIPTION

# Emacs TAGS file
CPP_SOURCES=src/*.cpp src/*.h r-package/dreadds/src/*.cpp
TAGS: $(CPP_SOURCES)
	etags $(CPP_SOURCES)

# R package
rp: r-package/$(R_PKG_TGZ)

r-package/$(R_PKG_TGZ): src/*.h src/libdreadds.a r-package/dreadds/src/*.cpp
	cd r-package/ && DREAD_DS_SRC=$$(cd ../src/ && pwd) R CMD build dreadds
	DREAD_DS_SRC=$$(cd src && pwd) R --no-save <<< 'install.packages("r-package/$(R_PKG_TGZ)", repos = NULL)'

