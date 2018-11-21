// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_OUTPUT_FILE_H
#define DREADDS_OUTPUT_FILE_H

#include <stdio.h>
#include <string>

#include "model-config.h"

namespace DreadDs {
  FILE *open_output_file(const Config &conf,
			 std::string filename_suffix);
}

#endif
