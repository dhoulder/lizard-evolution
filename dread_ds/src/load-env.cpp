// -*- coding: utf-8 -*-

#include "model.h"

using namespace DreadDs;
namespace DreadDs {

  std::unique_ptr <EnvMatrix> load_env(const filename_vec &env_inputs) {

    bool first = true;
    std::unique_ptr<EnvMatrix> em;
    for (const char *filename: env_inputs) {
      // FIXME STUB load cells from file

      // first file defines bb. subsequent files must match
      // first file creates em and v[0], subsequent files fill in v[1], v[2], ...

      if (first) {
	first = false;
	int rows = 5, cols = 6; // FIXME
	em = std::unique_ptr<EnvMatrix> (new EnvMatrix(boost::extents[rows][cols]));
      } else {
	// FIXME
      }
    }

    // FIXME fail if empty (first still true)

    return em;

  }
}
