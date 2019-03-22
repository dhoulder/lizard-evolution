\name{dreadds-package}
\alias{dreadds-package}
\alias{dreadds}
\docType{package}
\title{
Dynamic Range Evolution and Diversification simulation with dynamic speciation
}
\description{
  Implements a spatially explicit dynamic intraspecific
  macroevolutionary model to simulate expected spatial
  patterns of species and phylogenetic diversity. The primary
  purpose is to select appropriate models which
  generate expected spatial patterns best matching observed
  diversity. Used to determine the processes responsible for the
  spatial distribution of biological diversity and phylogenetic
  endemism.
}
\details{
  The package provides the createDreadDS() constructor. See the
  examples below and the documentation for createDreadDS().
}
\author{
  Model design and original work by Dan Rosauer <dan.rosauer@anu.edu.au>

  C++ and R implementation by David Houlder <david@davidhoulder.com>
}

\seealso{
  \code{\link{createDreadDS}}
  }

\references{
  \url{https://github.com/DanRosauer/DREaD_ds}
}
\keyword{ package }

\examples{
  \dontrun{
  > require(dreadds)
  Loading required package: dreadds
  > m <- createDreadDS('test.conf', c('-o', 'out'))
  # See help(createDreadDS) for details
}
}
