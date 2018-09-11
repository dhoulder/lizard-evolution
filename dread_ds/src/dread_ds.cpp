/**

 */
#include <iostream>

#include <Eigen/Dense>

#include "dread_ds.h"


using Eigen::MatrixXd;


struct DreadDs::Impl {
  MatrixXd env;
};

DreadDs::DreadDs(): impl(new(Impl)) {
  // STUB
  impl->env = MatrixXd::Random(3,3);
  std::cout << "constructor" << std::endl;
  std::cout << impl->env << std::endl;
}


DreadDs::~DreadDs() = default;
