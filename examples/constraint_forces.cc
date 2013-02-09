#include <iostream>
#include <Eigen/Dense>
#include "bicycle.h"
#include "whipple.h"

int main() {
  Bicycle b;
  Whipple w;
  b.set_parameters_from_whipple(w);
  //Eigen::Matrix<double, 6, 1, Eigen::RowMajor> forces;
  std::cout << b.compute_contact_forces();
}
