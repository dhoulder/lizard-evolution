// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_MODEL_LIMITS_H
#define DREADDS_MODEL_LIMITS_H

namespace DreadDs {
  static const char *version = "species-comp 0.1"; // NOTE: keep this up to date
  static const float gene_drift_sd = 1.0f; //Standard deviation of
  // Brownian motion gene drift per time-step
  static const float default_gene_flow_clip = 0.001;
  static const int max_env_dims = 2; // Max environment layers (e.g. 2
  // for temperature, precipitation)
  static const int max_genetic_dims= 5; // Max number of abstract genetic axes
}

#endif
