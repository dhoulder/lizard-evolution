// -*- coding: utf-8; Mode: c++ -*-
#ifndef DREADDS_ENV_PARAMS_H
#define DREADDS_ENV_PARAMS_H

#include <string>
#include <map>

#include <boost/filesystem.hpp>

namespace DreadDs {

  class Environment; // for opaque pointer below

  class BaseEnvParams {
  public:
    std::string grid_filename;

    BaseEnvParams(const std::string &fn):
      grid_filename(fn) {}
    virtual std::string get_filename(int step_offset) const= 0;
    virtual float update_environment(Environment *env,
                                     int step_offset,
                                     int layer) const = 0;
    virtual ~BaseEnvParams() = default;
  };

  class TsEnvParams: public BaseEnvParams {
    // for environment supplied as time series of grids
  public:
    std::string ts_dir;
    int ts_start = 0;
    int ts_step = 1;
    std::map<int, boost::filesystem::path> ts_dir_table;

    using BaseEnvParams::BaseEnvParams;

    virtual std::string get_filename(int step_offset) const;
    virtual float update_environment(Environment *env,
                                     int step_offset,
                                     int layer) const;
    void scan_ts_dir();
  };

  class ExEnvParams: public BaseEnvParams {
    // for environment specified as base + offset
  public:
    double ramp = 0.0; // linear environment change per time step
    float sine_period = 4.0f; // wavelength of sinusoidal environment
    // change in time steps
    float sine_offset = 0.0f; // shift sinusoidal change by this
    // fraction of the wavelength. Use 0.25 for cos()
    float sine_amplitude = 0.0f; // maximum swing of sinusoidal
                                 // environment change

    using BaseEnvParams::BaseEnvParams;
    virtual std::string get_filename(int step_offset) const;
    virtual float update_environment(Environment *env,
                                     int step_offset,
                                     int layer) const;
  };
}

#endif
