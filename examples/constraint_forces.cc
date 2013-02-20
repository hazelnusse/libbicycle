#include <iostream>
#include "bicycle.h"
#include "whipple.h"

int main()
{
  bicycle::Bicycle b;
  bicycle::Whipple w;

  b.set_parameters(w);
  b.solve_configuration_constraint_and_set_state();

  bicycle::Vector constraint_forces(7);
  constraint_forces = b.steady_constraint_forces();
  std::cout << constraint_forces << std::endl;
  std::cout << "Sum of normal forces divided by g: " << std::endl;
  std::cout << -(constraint_forces[2] + constraint_forces[5]) / w.g << std::endl;
  std::cout << "Total mass: " << std::endl;
  std::cout << w.mR + w.mB + w.mH + w.mF << std::endl;
}
