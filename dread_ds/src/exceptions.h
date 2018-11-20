// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_EXCEPTIONS_H
#define DREADDS_EXCEPTIONS_H

#include <stdexcept>

namespace DreadDs {

  class ApplicationException : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
  };

  class ConfigError : public ApplicationException {
  public:
    using ApplicationException::ApplicationException;
  };

  class UsageException: public ApplicationException {
  public:
    using ApplicationException::ApplicationException;
  };

}

#endif
