// -*- coding: utf-8 -*-

#include "model.h"

using namespace DreadDs;
namespace DreadDs {

  std::unique_ptr <Environment> load_env(const filename_vec &env_inputs) {

    bool first = true;
    std::unique_ptr<Environment> env;
    for (const char *filename: env_inputs) {
      // FIXME WIP load cells from file

      // first file defines bb. subsequent files must match
      // first file creates env and v[0], subsequent files fill in v[1], v[2], ...

      if (first) {
	first = false;
	int rows = 5, cols = 6; // FIXME STUB

	env = std::unique_ptr<Environment> (new Environment(rows, cols));


      } else {
	// FIXME WIP
      }
    }

    // FIXME WIP fail if empty (first still true)

    return env;

  }
}
