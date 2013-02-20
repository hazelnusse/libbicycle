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
  b.set_dependent_speeds({0, 2, 5});
  
  Matrix A_min(4, 4);
  A_min << 0, 0, 1, 0,
           0, 0, 0, 1,
           0, 0, 0, 0,
           0, 0, 0, 0;

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
  std::cout.precision(16);
  for (int i = 0; i < 11; ++i) { // iterate over all speeds in Meijaard et al.
    b.set_speed(4, -i/w.rR);
    b.solve_velocity_constraints_and_set_state();

    Matrix A_full = b.mass_matrix_full()
                     .block<14, 14>(0, 0)
                     .fullPivHouseholderQr()
                     .solve(b.independent_state_matrix().block<14, 10>(0, 0));

    A_min(2, 0) = A_full(9, 1);
    A_min(2, 1) = A_full(9, 2);
    A_min(2, 2) = A_full(9, 7);
    A_min(2, 3) = A_full(9, 8);
    A_min(3, 0) = A_full(11, 1);
    A_min(3, 1) = A_full(11, 2);
    A_min(3, 2) = A_full(11, 7);
    A_min(3, 3) = A_full(11, 8);

    std::cout << "v = " << i << ": "
      << A_min.eigenvalues().transpose() << std::endl;
  } // for i
} // main()

