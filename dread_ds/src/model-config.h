// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_MODEL_CONFIG_H
#define DREADDS_MODEL_CONFIG_H

#include <vector>
#include <string>

#include "model-limits.h"
#include "exceptions.h"

namespace DreadDs {

  struct EnvParams {
    std::string grid_filename;
    double ramp = 0.0; // linear environment change per time step
    float sine_period = 4.0f; // wavelength of sinusoidal environment
    // change in time steps
    float sine_offset = 0.0f; // shift sinusoidal change by this
    // fraction of the wavelength. Use 0.25
    // for cos()
    float sine_amplitude = 0.0f; // maximum swing of sinusoidal environment change
  };

  typedef std::vector<EnvParams> EnvParamsVec;

  struct NicheSpec {
    float centre;
    float tolerance;
  };

  struct SpeciesParameters {
    float max_dispersal_radius;
    std::vector<NicheSpec> niche;
    float genetics[max_genetic_dims] = {0.0f};

    float north, south, east, west; // Initial bounding box
  };


  class Config {
  public:
    // Parameters for a simulation run.
    int verbosity = 1;
    int env_dims = 0; // must be <= max_env_dims
    int genetic_dims = max_genetic_dims; // <= max_genetic_dims
    float gene_flow_clip = 0.001f;
    float gene_flow_zero_distance = 5.0f;
    float gene_drift_sd = 1.0f; // FIXME use timestep size
    float niche_evolution_rate = 0.1;
    float dispersal_min = 0.02;
    std::string output_dir = "./";
    std::string output_file_prefix = "";
    int n_iterations = 0;
    EnvParamsVec env_params;
    std::vector<SpeciesParameters> initial_species;

    Config() {}

    Config(int argc, const char *argv[]);
  };

}

#endif
