#include <iostream>
#include "bicycle.h"
#include "whipple.h"

int main() {
  Bicycle b;      // Initial state is zero, all parameters zero, 9.81 gravity
  Whipple w;      // Benchmark Whipple parameter values
  b.set_parameters_from_whipple(w);  // set Whipple parameters
  // default dependent coordinate is pitch angle
  b.solve_configuration_constraint_and_set_state();
  Bicycle::state x(Bicycle::state::Zero());

  {
    int dep_speeds[3] = {0, 2, 5};  // yaw, pitch, front wheel rates
    b.set_dependent_speeds(dep_speeds);
    x[12] = -1.0;
    b.set_state(x);
    b.solve_velocity_constraints_and_set_state();
    std::cout << b;
    std::cout << "Front wheel rate should be " << x[12]*w.rR/w.rF << std::endl;
  }
  {
    int dep_speeds[3] = {0, 2, 4};  // yaw, pitch, rear wheel rates
    b.set_dependent_speeds(dep_speeds);
    x[13] = -1.0;
    b.set_state(x);
    b.solve_velocity_constraints_and_set_state();
    std::cout << b;
    std::cout << "Rear wheel rate should be " << x[13]*w.rF/w.rR << std::endl;
  }
}
