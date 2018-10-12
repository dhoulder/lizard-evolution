// -*-c++-*-

#ifndef DREADDS_MODEL_CONFIG_H
#define DREADDS_MODEL_CONFIG_H

namespace DreadDs {

  static const int max_env_dims = 2; // Max environment layers (e.g. 2
  // for temperature, precipitation)
  static const int max_genetic_dims= 3; // Max number of abstract genetic axes


  struct EnvChange {
    double ramp = 0.0; // linear environment change per time step
    float sine_period = 4.0f; // wavelength of sinusoidal environment
    // change in time steps
    float sine_offset = 0.0f; // shift sinusoidal change by this
    // fraction of the wavelength. Use 0.25
    // for cos()
    float sine_amplitude = 0.0f; // maximum swing of sinusoidal environment change
  };


  struct Config {
    // Parameters for a simulation run.
    int debug = 4; // FIXME
    int env_dims = 1; // must be <= max_env_dims
    int genetic_dims = max_genetic_dims; // <= max_genetic_dims
    EnvChange env_change[max_env_dims];

    float gene_flow_threshold = 0.001f;
    float gene_flow_zero_distance = 5.0f;
    float niche_evolution_rate = 0.1;
    float gene_drift_sd = 1.0f; // FIXME use timestep size
  };

}

#endif
