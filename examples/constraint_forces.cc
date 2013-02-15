#include <iostream>
#include "bicycle.h"
#include "whipple.h"

using namespace bicycle;

int main() {
  Bicycle b;      // Initial state is zero, all parameters zero, 9.81 gravity
  Whipple w;      // Benchmark Whipple parameter values
  b.set_parameters_from_whipple(w);  // set Whipple parameters
  b.solve_configuration_constraint_and_set_state();
  std::cout << b << std::endl;
  std::cout << b.steady_constraint_forces() << std::endl;
  
  Bicycle::coordinates q = Bicycle::coordinates::Zero();
  q[1] = 0.02;
  q[3] = 0.03;
  b.set_coordinates(q);
  Bicycle::speeds u = Bicycle::speeds::Zero();
  u[4] = -1.0/0.3;
  b.set_speeds(u);
  b.solve_velocity_constraints_and_set_state();
  std::cout << b << std::endl;
  auto cf = b.steady_constraint_forces();
  std::cout << cf << std::endl;
  std::cout << (cf(2, 0) + cf(5, 0)) / -9.81 << std::endl;
}
