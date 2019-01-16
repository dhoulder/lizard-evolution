// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_MODEL_CONFIG_H
#define DREADDS_MODEL_CONFIG_H

#include <vector>
#include <memory>
#include <string>
#include <cmath>

#include "constants.h"
#include "env-params.h"

namespace DreadDs {

  class Environment; // opaque reference below

  struct NicheSpec {
    NicheSpec(float c, float t):
      centre(c),
      tolerance(t) {}

    float centre;
    float tolerance;
  };

  struct SpeciesParameters {
    float max_dispersal_radius;
    float genetics[max_genetic_dims] = {0.0f};

    enum NicheMode {NICHE_CENTRE, RANDOM_CELL};
    NicheMode niche_mode;
    std::vector<NicheSpec> niche;
    enum BoundsMode {NSEW, RANDOM_SIZE};
    BoundsMode bounds_mode;
    float random_rect_min = NAN, random_rect_max = NAN;
    float north, south, east, west; // Initial bounding box
  };

  class Config {
  private:
    std::vector<SpeciesParameters> initial_species;
    void set_params_from_env(SpeciesParameters &sp,
			     const Environment &env) const;
  public:
    // Parameters for a simulation run.
    int verbosity = 1;
    int env_dims = 0; // must be <= max_env_dims
    int genetic_dims = max_genetic_dims; // <= max_genetic_dims
    float gene_flow_clip = 0.001f;
    float gene_flow_zero_distance = 5.0f;
    float niche_evolution_rate = 0.1;
    float dispersal_min = 0.02;
    std::string output_dir = "./";
    std::string output_file_prefix = "";
    int n_iterations = 0;
    std::vector<std::shared_ptr<BaseEnvParams>> env_params;
    int check_speciation = 1;

    Config() {}

    Config(int argc, const char *argv[]);
    std::vector<SpeciesParameters> get_initial_species(
      const Environment &env) const;
  };
}

#endif
