// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_EXCEPTIONS_H
#define DREADDS_EXCEPTIONS_H

#include <stdexcept>

namespace DreadDs {

  class ApplicationError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
  };

  class ConfigError : public ApplicationError {
  public:
    using ApplicationError::ApplicationError;
  };

}

#endif
