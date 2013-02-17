#include <iostream>
#include "bicycle.h"
#include "whipple.h"

using namespace bicycle;

int main() {
  Bicycle b;      // Initial state is zero, all parameters zero, 9.81 gravity
  Whipple w;      // Benchmark Whipple parameter values
  b.set_parameters(w);  // set Whipple parameters
  
  Bicycle::coordinates q = Bicycle::coordinates::Zero();
  q[1] = 0.02;
  q[3] = 0.03;
  b.set_coordinates(q);
  // default dependent coordinate is pitch angle
  b.solve_configuration_constraint_and_set_state();
  //std::cout << b;
  Bicycle::speeds u = Bicycle::speeds::Zero();
  u[4] = -1.0/0.3;
  b.set_speeds(u);
  std::cout << b;
  b.solve_velocity_constraints_and_set_state();
  std::cout << b;
  std::cout << b.mass_matrix_full().block<12, 8>(8, 0);
  std::cout << b.independent_state_matrix();
}
