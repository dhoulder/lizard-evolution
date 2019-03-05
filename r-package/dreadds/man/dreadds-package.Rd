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
  The package provides the dreadds_model() constructor. See the
  examples below and the documentation for dreadds_model().
}
\author{
David Houlder, david@davidhoulder.com, david.houlder@anu.edu.au
}
\references{
  https://github.com/DanRosauer/DREaD_ds
}
\keyword{ package }

\examples{
  \dontrun{
  > require(dreadds)
  Loading required package: dreadds
  > m <- dreadds_model('test.conf', c('-o', 'out'))
  # See help(dreadds_model) for details
}
}
