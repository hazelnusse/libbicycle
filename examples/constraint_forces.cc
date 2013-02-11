#include <iostream>
#include "bicycle.h"
#include "whipple.h"

int main() {
  Bicycle b;      // Initial state is zero, all parameters zero, 9.81 gravity
  Whipple w;      // Benchmark Whipple parameter values
  b.set_parameters_from_whipple(w);  // set Whipple parameters
  // default dependent coordinate is pitch angle
  b.solve_configuration_constraint_and_set_state();
  std::cout << b << std::endl;
  std::cout << b.steady_contact_forces() << std::endl;
  // this should fail with an error to cerr but keep the state unchanged.
  b.set_dependent_coordinate(1);
  b.solve_configuration_constraint_and_set_state();
  b.set_dependent_coordinate(3);
  b.solve_configuration_constraint_and_set_state();
}
