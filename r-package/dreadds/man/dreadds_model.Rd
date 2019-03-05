\name{dreadds_model}
\alias{dreadds_model}
\docType{package}
\title{
Dynamic Range Evolution and Diversification simulation with dynamic speciation model
}

\usage{
 m <- dreadds_model(path_to_config_file, configuration_options = character(0))
}

\arguments{
\item{path_to_config_file}{
  Pathname of DREaD_ds model configuration file. See examples/*.conf in the source
  repository}

\item{configuration_options}{
  Character vector of additional configuration options. e.g. c('-o', 'out-dir')
  Run "dreadds --help" for a complete list
}
}

\value{A dreaddsModel instance. See below for methods}

\description{
Creates an instance of a DREaD_ds model.
}

\section{Methods}{
  \describe{
    \item{$do_step()}{

    Perform one step of the simulation. Returns the total number of occupied cells}

    \item{$save()}{

    Save output files. Not strictly needed when using the R interface as all of
    the model's internal state can be accessed using model methods}

    \item{$get_env()}{

    Returns a list of matrices containing the current values of the simulation environment
    at each grid cell.
    }

    \item{$get_demes()}{

    Returns a list of data frames describing the population at each grid cell. Each data frame
    corresponds to an extant species.
    }

    \item{$get_phylogeny()}{

    Returns the current species phylogeny as a character vector using Newick format.
    }

    \item{$get_species()}{

    Returns a data frame describing the characteristics of all species in the model.
    }
}
}

\examples{
\dontrun{
  m <- dreadds_model('test.conf', c('-o', 'out', '--niche-evolution-rate', 0.02))
  > m$do_step()
  [1] 14
  > m$do_step()
  [1] 9
  > m$do_step()
  [1] 6
  > m$save() # writes files
  > environment_matrices <- m$get_env()
  > demes <- m$get_demes()
  > m$get_phylogeny()
  [1] "species_1;"
  > species_stats <- m$get_species()
  }
}
