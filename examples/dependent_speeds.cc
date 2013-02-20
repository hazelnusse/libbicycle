#include <iostream>
#include "bicycle.h"
#include "whipple.h"

int main()
{
  bicycle::Bicycle b;
  bicycle::Whipple w;
  b.set_parameters(w);
  b.solve_configuration_constraint_and_set_state();

  b.set_dependent_speeds({0, 2, 5});  // yaw, pitch, front wheel rates
  b.set_speed(4, -1.0);               // rear wheel rate
  std::cout << "Residual of velocity constraints: \n"
      << b.solve_velocity_constraints_and_set_state() << std::endl;
  std::cout << "Calculated front wheel rate: " << b.speed(5) << std::endl;
  std::cout << "Front wheel rate should be " << b.speed(4)*w.rR/w.rF << std::endl;

  b.set_dependent_speeds({0, 2, 4});  // yaw, pitch, rear wheel rates
  b.set_speed(5, -1.0);               // front wheel rate
  std::cout << "Residual of velocity constraints: \n"
      << b.solve_velocity_constraints_and_set_state() << std::endl;
  std::cout << "Calculated rear wheel rate: " << b.speed(4) << std::endl;
  std::cout << "Rear wheel rate should be " << b.speed(5)*w.rF/w.rR << std::endl;
}
