#include <iostream>
#include <Eigen/Dense>
#include "bicycle.h"
#include "whipple.h"

using namespace bicycle;

int main()
{
  Bicycle b;
  Whipple w;
  b.set_parameters(w);
  b.solve_configuration_constraint_and_set_state();
  b.set_speed(4, -1.0/w.rR);
  b.solve_velocity_constraints_and_set_state();
  
  Vector r = Vector::Zero(22);
  Vector cf = b.steady_constraint_forces();
  r[4] = cf[0];
  r[5] = cf[1];
  r[6] = cf[2];
  r[14] = cf[3];
  r[15] = cf[4];
  r[16] = cf[5];
  r[20] = cf[6];
  r[21] = 9.81;
  b.set_inputs(r);

  std::cout << "State derivatives:" << std::endl;
  std::cout << b.state_derivatives() << std::endl;
} // main()


