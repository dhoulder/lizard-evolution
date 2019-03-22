\name{createDreadDS}
\alias{createDreadDS}
\docType{package}
\title{
  Dynamic Range Evolution and Diversification simulation with dynamic speciation model
}

\usage{
 m <- createDreadDS(config.file = NULL,
                    genetic.dims = NULL,
                    gene.flow.max.distance = NULL,
                    niche.evolution.rate = NULL,
                    dispersal.min = NULL,
                    output.dir = NULL,
                    output.file.prefix = NULL,
                    iterations = NULL,
                    check.speciation = NULL
                    args = character(0))
}

\arguments{
  \item{config.file}{
    Pathname of DREaD_ds model configuration file. See examples/*.conf in the source
    repository.
  }
  \item{genetic.dims}{Overrides value from configuration file}
  \item{gene.flow.max.distance}{Overrides value from configuration file}
  \item{niche.evolution.rate}{Overrides value from configuration file}
  \item{dispersal.min}{Overrides value from configuration file}
  \item{output.dir}{Overrides value from configuration file}
  \item{output.file.prefix}{Overrides value from configuration file}
  \item{iterations}{Overrides value from configuration file}
  \item{check.speciation}{Overrides value from configuration file}
  \item{args}{
    Character vector of additional configuration options. e.g. c(' --csv-precision', 4)
    Run the standalone model executable as "dreadds --help" for a complete list
  }
}

\value{A model instance. See below for methods}

\description{
  Creates an instance of a DREaD_ds model.
}

\section{Methods}{
  \describe{

    \item{$runAll()}{
      Runs the model until its configured iteration limit is reached or
      the total population reaches zero.

      Returns -1 if no steps could be performed due to the configured
      iteration limit being reached.  In all other circumstances, the
      total cell count after the final step is returned.}
    
    \item{$runSteps(iterations)}{
      Perform one or more steps of the simulation.

      Performs up to the specified number of steps, stopping prematurely
      only if the total population reaches zero. The configured
      iteration limit is ignored.

      Returns -1 if iterations is less than 1. In all other
      circumstances, the total cell count after the last step is
      returned.}

    \item{$save()}{
      Save CSV and YAML output files. Not strictly needed when using the
      R interface as all of the model's internal state can be accessed
      using model methods}

    \item{$getEnv()}{
      Returns a list of matrices containing the current values of the
      simulation environment at each grid cell. The first row of each
      matrix corresponds to the northern edge of the simulation area.}

    \item{$getExtent()}{
      Return the geographical bounding box of the simulation area's cell
      edges in x-min, x-max, y-min, y-max order.}

    \item{$getDemes()}{
      Returns a list of data frames describing the population at each
      grid cell. Each data frame corresponds to an extant species.}

    \item{$getPhylogeny()}{
      Returns the current species phylogeny as a character vector using
      Newick format.}

    \item{$getSpecies()}{
      Returns a data frame describing the characteristics of all species
      in the model.}
  }
}

\section{Properties}{
  \describe{
    \item{$step}{
      The current step number of the model. Will be 0 before the first
      call to $runSteps()}
  }
}

\examples{
\dontrun{
  m <- createDreadDS('test.conf',  output.dir = 'out', niche.evolution.rate =  0.02)
  > m$runSteps(1)
  [1] 14
  > m$runSteps(1)
  [1] 9
  > m$runSteps(1)
  [1] 6
  > m$save() # writes files
  > environment_matrices <- m$getEnv()
  > demes <- m$getDemes()
  > m$getPhylogeny()
  [1] "species_1;"
  > species_stats <- m$getSpecies()
  }
}
