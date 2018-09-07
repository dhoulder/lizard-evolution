/**

 */

#include <Eigen/Dense>

#include "dread_ds.h"


using Eigen::MatrixXd;


struct DreadDs::Impl {
  MatrixXd env;
};

DreadDs::DreadDs(): impl(new(Impl)) {
   // STUB
}


DreadDs::~DreadDs() = default;
