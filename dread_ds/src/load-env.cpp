// -*- coding: utf-8 -*-

#include "model.h"

using namespace DreadDs;
namespace DreadDs {

  std::unique_ptr <EnvMatrix> load_env(const filename_vec &env_inputs) {

    for (const char *filename: env_inputs) {
      // FIXME STUB
      // load cells from file
    }

    int rows = 5, cols = 6;
    return std::unique_ptr<EnvMatrix> (new EnvMatrix(boost::extents[rows][cols]));
  }
}
