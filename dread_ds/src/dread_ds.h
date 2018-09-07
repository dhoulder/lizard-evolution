/**
 */

#ifndef DREADDS_H
#define DREADDS_H

#include <memory>

class DreadDs  {
 public:
  const char *const version = "0.01";
  DreadDs();
  ~DreadDs();

 private:
  class Impl;
  std::unique_ptr<Impl> impl;
};

#endif
