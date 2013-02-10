#include "bicycle.h"
using namespace Eigen;

Matrix<double, 6, 1> Bicycle::compute_contact_forces() const
{
  Matrix<double, 12, 1> gif;
  Matrix<double, 12, 22, RowMajor> gaf_dr_full;
  
  gif_ud_zero(gif.data());    // compute coriolis/centripetal inertia forces
  gaf_dr(gaf_dr_full.data()); // compute generalized active force input matrix

  // TODO: modify to allow for external forces, external torques, and internal
  // torques to be accounted for.  The current implementation assumes
  // everything besides contact forces is equal to zero.
  Matrix<double, 12, 7, RowMajor> gaf_dr_reduced;
  gaf_dr_reduced << gaf_dr_full.block<12, 3>(0, 4),  // Rear contact forces
                    gaf_dr_full.block<12, 3>(0, 14), // Front contact forces
                    gaf_dr_full.block<12, 1>(0, 21); // Gravitational term

  Matrix<double, 3, 7, RowMajor> gaf_dr_d; // Dependent active forces
  Matrix<double, 9, 7, RowMajor> gaf_dr_i; // Independent active forces
  Matrix<double, 3, 1> gif_ud_zero_d;      // Dependent coriolis/centripetal
  Matrix<double, 9, 1> gif_ud_zero_i;      // Independent coriolis/centripetal

  for (int i = 0, id = 0, ii = 0; i < 12; ++i) {
    if (i == dependent_speeds_[0] ||
        i == dependent_speeds_[1] ||
        i == dependent_speeds_[2]) {
      gaf_dr_d.block<1, 7>(id, 0) = gaf_dr_reduced.block<1, 7>(i, 0);
      gif_ud_zero_d[id] = gif[i];
      ++id;
    } else {
      gaf_dr_i.block<1, 7>(ii, 0) = gaf_dr_reduced.block<1, 7>(i, 0);
      gif_ud_zero_i[ii] = gif[i];
      ++ii;
    }
  }

  Matrix<double, 9, 3, RowMajor> C(compute_Bd_inverse_Bi().transpose());
  Matrix<double, 9, 7, RowMajor> gaf_dr_constrained;
  Matrix<double, 9, 1> rhs;
  gaf_dr_constrained = (gaf_dr_i + (C*gaf_dr_d));
  rhs = - (gif_ud_zero_i
         + C*gif_ud_zero_d
         + gaf_dr_constrained.block<9, 1>(0, 6)*g_);

  return gaf_dr_constrained.block<6, 6>(3, 0)
          .fullPivHouseholderQr().solve(rhs.block<6, 1>(3, 0));
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

void Bicycle::solve_velocity_constraints_and_set_state() {
  Matrix<double, 3, 9, RowMajor> C(compute_Bd_inverse_Bi());
  Matrix<double, 9, 1> u_independent;
  for (int i = 0, ii = 0; i < 12; ++i) {
    if (i == dependent_speeds_[0] ||
        i == dependent_speeds_[1] ||
        i == dependent_speeds_[2]) {
      continue;
    } else {
      u_independent[ii++] = state_[i + 8];
    }
  }
  Matrix<double, 3, 1> u_dependent(C * u_independent);

  state_[dependent_speeds_[0] + 8] = u_dependent[0];
  state_[dependent_speeds_[1] + 8] = u_dependent[1];
  state_[dependent_speeds_[2] + 8] = u_dependent[2];
}

Matrix <double, 3, 9> Bicycle::compute_Bd_inverse_Bi() const {
  Matrix<double, 3, 12, RowMajor> B;
  f_v_coefficient(B.data());  // compute velocity constraint coefficient matrix
  Matrix<double, 3, 3, RowMajor> B_d;      // Dependent columns
  Matrix<double, 3, 9, RowMajor> B_i;      // Independent columns
  for (int i = 0, ii = 0, id = 0; i < 12; ++i) {
    if (i == dependent_speeds_[0] ||
        i == dependent_speeds_[1] ||
        i == dependent_speeds_[2]) {
      B_d.block<3, 1>(0, id++) = B.block<3, 1>(0, i);
    } else {
      B_i.block<3, 1>(0, ii++) = B.block<3, 1>(0, i);
    }
  }
  FullPivHouseholderQR <Matrix<double, 3, 3> > decomposition;
  decomposition.compute(B_d);  // compute decomposition of B_d
  Matrix<double, 3, 9, RowMajor> C;
  for (int i = 0; i < 9; ++i) {
    C.block<3, 1>(0, i) = -decomposition.solve(B_i.block<3, 1>(0, i));
  }
  return C;
}
