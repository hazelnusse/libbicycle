#include "bicycle.h"
//#include <algorithm>

namespace bicycle {

RowMajorMatrix Bicycle::steady_constraint_forces() const
{
 // indices of ground contact forces and steer torque
  const std::set<int> cf_indices = {4, 5, 6, 14, 15, 16, 20};
  RowMajorMatrix gif_steady(12, 1),    // Generalized inertia forces
                 gaf_dr_full(12, 22);  // Input coefficient matrix
  gif_ud_zero(gif_steady.data());      // populate gif_steady
  gaf_dr(gaf_dr_full.data());          // populate gaf_dr_full
  // put rows associated with dependent speeds at the bottom
  gif_steady = P_u_.transpose() * gif_steady;
  gaf_dr_full = P_u_.transpose() * gaf_dr_full;
  
  // constraint coefficient matrix associated with auxilliary speeds
  RowMajorMatrix C_aux(7, 3);
  {
    RowMajorMatrix tmp = Bd_inverse_Bi().transpose(); // 9 by 3
    C_aux.block<1, 3>(0, 0) = tmp.block<1, 3>(1, 0); //  steer rate
    C_aux.block<6, 3>(1, 0) = tmp.block<6, 3>(3, 0); //  contact point rates
  }

  Eigen::Matrix<int, 22, 1> indices;
  std::iota(indices.data(), indices.data() + 22, 0); // 0,1,...,21
  std::stable_partition(indices.data(), indices.data() + 22,
      [&cf_indices](int elem) { return !cf_indices.count(elem); });
  // Permute columns so constraint forces are at end
  gaf_dr_full = gaf_dr_full * PermutationMatrix<22>(indices);

  RowMajorMatrix A(7, 7);
  A.block<6, 7>(0, 0) = gaf_dr_full.block<6, 7>(3, 15)
                + C_aux.block<6, 3>(1, 0) * gaf_dr_full.block<3, 7>(9, 15);
  A.block<1, 7>(6, 0) = gaf_dr_full.block<1, 7>(1, 15)
                + C_aux.block<1, 3>(0, 0) * gaf_dr_full.block<3, 7>(9, 15);

  RowMajorMatrix b(7, 1);
  // Inertia terms
  b.block<6, 1>(0, 0) = - (gif_steady.block<6, 1>(3, 0)
                + C_aux.block<6, 3>(1, 0) * gif_steady.block<3, 1>(9, 0));
  b.block<1, 1>(6, 0) = - (gif_steady.block<1, 1>(1, 0)
                + C_aux.block<1, 3>(0, 0) * gif_steady.block<3, 1>(9, 0));
  // b = - (gif_steady.block<7, 1>(3, 0) + C_aux * gif_steady.block<3, 1>(9, 0));

  b.block<6, 1>(0, 0) -= (gaf_dr_full.block<6, 15>(3, 0)
                + C_aux.block<6, 3>(1, 0) * gaf_dr_full.block<3, 15>(9, 0)) * all_inputs_except_contact_forces();
  b.block<1, 1>(6, 0) -= (gaf_dr_full.block<1, 15>(1, 0)
                + C_aux.block<1, 3>(0, 0) * gaf_dr_full.block<3, 15>(9, 0)) * all_inputs_except_contact_forces();

  // Active force terms
  //b -= (gaf_dr_full.block<6, 16>(3, 0)
  //      + C_aux * gaf_dr_full.block<3, 16>(9, 0))*all_inputs_except_contact_forces();
  // TODO: Add term involving mass matrix and du/dt for non-steady case

  return A.fullPivHouseholderQr().solve(b);
}

void Bicycle::solve_configuration_constraint_and_set_state(double ftol, int iter) {
  const double df_min = 1e-14;
  Matrix<double, 8, 1> df;
  double f, q_d_prev = state_[dependent_coordinate_]; // initial state

  int i = 0;
  do {
    f_c(&f);                        // evaluate f
    f_c_dq(df.data());              // evaluate df
    
    if (fabs(df[dependent_coordinate_]) < 1e-14) {
      std::cerr << "Derivative w.r.t. dependent coordinate q["
                << dependent_coordinate_ << "] is less than " << df_min << "."
                << std::endl << "This indicates this coordinate cannot "
                "effectively raise the front wheel contact point and a different"
                << std::endl << "coordinate should be selected as the "
                "dependent coordinate.  The coordinate has not been changed."
                << std::endl;
      state_[dependent_coordinate_] = q_d_prev; // Restore the coordinate
      break;
    }
    state_[dependent_coordinate_] -= f/df[dependent_coordinate_]; // Newton step
  } while ( (++i < iter) & (fabs(f) > ftol) );
}

void Bicycle::solve_velocity_constraints_and_set_state()
{
  Matrix<double, 9, 1> u_independent;
  u_independent = (state_.block<12, 1>(8, 0).transpose() * P_u_).block<1, 9>(0, 0).transpose();
  Matrix<double, 3, 1> u_dependent(Bd_inverse_Bi() * u_independent);
  int i = 0;
  for (auto it = dependent_speeds_.begin();
       it != dependent_speeds_.end(); ++it) {
    state_[*it + kNumberOfCoordinates] = u_dependent(i++, 0);
  }
}

RowMajorMatrix Bicycle::Bd_inverse_Bi() const
{
  RowMajorMatrix B(3, 12);
  f_v_coefficient(B.data());  // compute velocity constraint coefficient matrix
  B = (B*P_u_).eval();        // move dependent columns to the end
  FullPivHouseholderQR<RowMajorMatrix> decomposition;
  decomposition.compute(B.block<kNumberOfVelocityConstraints,
    kNumberOfVelocityConstraints>(0,
              kNumberOfSpeeds - kNumberOfVelocityConstraints));  // compute decomposition of B_d
  RowMajorMatrix C(kNumberOfVelocityConstraints, kNumberOfSpeeds - kNumberOfVelocityConstraints);
  for (int i = 0; i < kNumberOfSpeeds - kNumberOfVelocityConstraints; ++i) {
    C.block<kNumberOfVelocityConstraints, 1>(0, i) = decomposition.solve(-B.block<kNumberOfVelocityConstraints, 1>(0, i));
  }
  return C;
}

bool Bicycle::is_dependent_index(int i) const
{
  return dependent_speeds_.count(i);
}

} // namespace bicycle
