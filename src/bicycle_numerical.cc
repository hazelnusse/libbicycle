#include "bicycle.h"

namespace bicycle {

using namespace Eigen;

Matrix<double, 6, 1> Bicycle::steady_contact_forces() const
{
  Matrix<double, 12, 1> gif;
  Matrix<double, 12, 22, RowMajor> gaf_dr_full;
  
  gif_ud_zero(gif.data());    // compute coriolis/centripetal inertia forces
  gaf_dr(gaf_dr_full.data()); // compute generalized active force input matrix

  // TODO: modify to allow for external forces, external torques, and internal
  // torques to be accounted for.  The current implementation assumes
  // everything besides contact forces is equal to zero.
  // TODO: modify to allow forces to be formed without assuming du/dt = 0
  Matrix<double, 12, 7, RowMajor> gaf_dr_reduced;
  gaf_dr_reduced << gaf_dr_full.block<12, 3>(0, 4),  // Rear contact forces
                    gaf_dr_full.block<12, 3>(0, 14), // Front contact forces
                    gaf_dr_full.block<12, 1>(0, 21); // Gravitational term

  Matrix<double, 3, 7, RowMajor> gaf_dr_d; // Dependent active forces
  Matrix<double, 9, 7, RowMajor> gaf_dr_i; // Independent active forces
  Matrix<double, 3, 1> gif_ud_zero_d;      // Dependent coriolis/centripetal
  Matrix<double, 9, 1> gif_ud_zero_i;      // Independent coriolis/centripetal

  for (int i = 0, id = 0, ii = 0; i < 12; ++i) {
    if (is_dependent_index(i)) {
      gaf_dr_d.block<1, 7>(id, 0) = gaf_dr_reduced.block<1, 7>(i, 0);
      gif_ud_zero_d[id] = gif[i];
      ++id;
    } else {
      gaf_dr_i.block<1, 7>(ii, 0) = gaf_dr_reduced.block<1, 7>(i, 0);
      gif_ud_zero_i[ii] = gif[i];
      ++ii;
    }
  }

  Matrix<double, 9, 3, RowMajor> C(Bd_inverse_Bi().transpose());
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

void Bicycle::solve_velocity_constraints_and_set_state()
{
  Matrix<double, 3, 9, RowMajor> C(Bd_inverse_Bi());
  Matrix<double, 9, 1> u_independent;
  for (int i = 0, ii = 0; i < 12; ++i) {
    if (is_dependent_index(i)) {
      continue;
    } else {
      u_independent[ii++] = state_[i + kNumberOfCoordinates];
    }
  }
  Matrix<double, 3, 1> u_dependent(C * u_independent);

  state_[dependent_speeds_[0] + kNumberOfCoordinates] = u_dependent[0];
  state_[dependent_speeds_[1] + kNumberOfCoordinates] = u_dependent[1];
  state_[dependent_speeds_[2] + kNumberOfCoordinates] = u_dependent[2];
}

Matrix<double, 3, 9, RowMajor> Bicycle::Bd_inverse_Bi() const
{
  Matrix<double, 3, 12, RowMajor> B;
  f_v_coefficient(B.data());  // compute velocity constraint coefficient matrix
  Matrix<double, 3, 3, RowMajor> B_d;      // Dependent columns
  Matrix<double, 3, 9, RowMajor> B_i;      // Independent columns
  for (int i = 0, ii = 0, id = 0; i < 12; ++i) {
    if (is_dependent_index(i)) {
      B_d.block<3, 1>(0, id++) = B.block<3, 1>(0, i);
    } else {
      B_i.block<3, 1>(0, ii++) = B.block<3, 1>(0, i);
    }
  }
  FullPivHouseholderQR <Matrix<double, 3, 3, RowMajor> > decomposition;
  decomposition.compute(B_d);  // compute decomposition of B_d
  Matrix<double, 3, 9, RowMajor> C;
  for (int i = 0; i < 9; ++i) {
    C.block<3, 1>(0, i) = decomposition.solve(-B_i.block<3, 1>(0, i));
  }
  return C;
}

Matrix<double, 8, 8, RowMajor> Bicycle::M_qq() const
{
  // \nabla_q f_0
  return Matrix<double, 8, 8>::Identity();
}

Matrix<double, 3, 8, RowMajor> Bicycle::M_uqc() const
{
  // \nabla_{dq/dt} f_a
  VectorXd B_dq(108, 1);           // Dynamically allocated
  f_v_coefficient_dq(B_dq.data()); // Populate the raw data
  Matrix<double, 3, 8, RowMajor> B_j_dq(Matrix<double, 3, 8, RowMajor>::Zero());

  for (int j = 0; j < kNumberOfSpeeds; ++j) { // column index
     B_j_dq.block<3, 3>(0, 1) += state_[j + kNumberOfCoordinates] *
        Map<Matrix<double, 3, 3, RowMajor>,
            Unaligned,
            Stride<3*kNumberOfSpeeds, 1> >(&B_dq[3*j]);
  }
  return B_j_dq;
}

Matrix<double, 3, 12, RowMajor> Bicycle::M_uuc() const
{
  // \nabla_{du/dt} f_a
  Matrix<double, 3, 12, RowMajor> B;
  f_v_coefficient(B.data());
  return B;
}

Matrix<double, 9, 8, RowMajor> Bicycle::M_uqd() const
{
  // \nabla_{dq/dt} f_3
  return Matrix<double, 9, 8, RowMajor>::Zero();
}

Matrix<double, 9, 12, RowMajor> Bicycle::M_uud() const
{
  // \nabla_{du/dt} f_2
  // This is the constrained du/dt coefficient matrix from Kane's equations
  Matrix<double, kNumberOfSpeeds, kNumberOfSpeeds, RowMajor> m;
  gif_dud(m.data());
  Matrix<double, 3, 12, RowMajor> md;
  Matrix<double, 9, 12, RowMajor> mi;
  for (int i = 0, ii = 0, id = 0; i < kNumberOfSpeeds; ++i) {
    if (is_dependent_index(i)) {
      md.block<1, 12>(id++, 0) = m.block<1, 12>(i, 0);
    } else {
      mi.block<1, 12>(ii++, 0) = m.block<1, 12>(i, 0);
    }

  }
  return mi + Bd_inverse_Bi().transpose() * md;
}

Matrix<double, 8, 8, RowMajor> Bicycle::A_qq() const
{
  // - \nabla_{q} (f_0 + f_1)
  // Since f_0 = I_{8x8} * \dot{q}, only the f_1 term appears in the gradient
  Matrix<double, 8, 8, RowMajor> m;
  f_1_dq(m.data());
  return -m;
}

Matrix<double, 8, 12, RowMajor> Bicycle::A_qu() const
{
  // - \nabla_u f_1
  Matrix<double, 8, 12, RowMajor> m;
  f_1_du(m.data());
  return -m;
}
  
Matrix<double, 3, 8, RowMajor> Bicycle::A_uqc() const
{
  // - \nabla_q f_a
  // TODO:  Rederive without assuming lean rate, pitch rate and steer rate are
  // zero.
  // NOTE: This is zero when lean rate, pitch rate, and steer rate are zero.
  // For non-steady equilibriums, this is not the case.
  return Matrix<double, 3, 8, RowMajor>::Zero();
}

Matrix<double, 3, 12, RowMajor> Bicycle::A_uuc() const
{
  // TODO: Verify that this returns a matrix of zeros when lean rate, pitch
  // rate, and steer rate are all zero (steady turn or steady forward cruise)
  VectorXd B_dq(108, 1);           // Dynamically allocated
  f_v_coefficient_dq(B_dq.data()); // Populate the raw data
  Matrix<double, 3, 12, RowMajor> B_j_dq(Matrix<double, 3, 12, RowMajor>::Zero());

  for (int j = 0; j < kNumberOfSpeeds; ++j) { // column index
     B_j_dq.block<3, 1>(0, j) = Map<Matrix<double, 3, 3, RowMajor>,
                                    Unaligned,
                                    Stride<3*kNumberOfSpeeds, 1> >(&B_dq[3*j])
                                      * state_.block<3, 1>(kNumberOfCoordinates + 1, 0);
  }
  return B_j_dq;
   
  // Below is the result when it is assumed that lean rate, pitch rate, and
  // steer rate are all zero.
  // return Matrix<double, 3, 12, RowMajor>::Zero();
}

Matrix<double, 9, 8, RowMajor> Bicycle::A_uqd() const
{
  // TODO: Rederive without assuming du/dt = 0
  // - \nabla_q (f_2 + f_3)
  Matrix<double, 9, 8, RowMajor> m = Matrix<double, 9, 8, RowMajor>::Zero();

  Matrix<double, Dynamic, Dynamic, RowMajor> m_gaf_dq(12, 8);
  gaf_dq(m_gaf_dq.data());
  Matrix<double, Dynamic, Dynamic, RowMajor> m_gif_ud_zero_dq(12, 8);
  gif_ud_zero_dq(m_gif_ud_zero_dq.data());

  for (int i = 0, id = 0, ii = 0; i < kNumberOfSpeeds; ++i) {
    if (is_dependent_index(i)) {
      ++id;
      continue;
    } else {
      m.block<1, 8>(ii, 0) -= (m_gaf_dq.block<1, 8>(i, 0)
                             + m_gif_ud_zero_dq.block<1, 8>(i, 0));
      ++ii;
    }
  }

  //for (int i = 0; i < 3; ++i) {
  //  Matrix<double, 3
  //}
  return m;
}

Matrix<double, 9, 12, RowMajor> Bicycle::A_uud() const
{
  // - \nabla_u (f_3)
  Matrix<double, 9, 12, RowMajor> m = Matrix<double, 9, 12, RowMajor>::Zero();

  return m;

}

Matrix<double, 20, 20, RowMajor> Bicycle::mass_matrix_full() const
{
  Matrix<double, 20, 20, RowMajor> m;
  m.block<kNumberOfCoordinates, kNumberOfCoordinates>(0, 0) = M_qq();

  m.block<kNumberOfCoordinates, kNumberOfSpeeds>(0, kNumberOfCoordinates)
    = Matrix<double, kNumberOfCoordinates, kNumberOfSpeeds>::Zero();

  m.block<kNumberOfVelocityConstraints, kNumberOfCoordinates>
      (kNumberOfCoordinates, 0) = M_uqc();
  
  m.block<kNumberOfVelocityConstraints, kNumberOfSpeeds>(kNumberOfCoordinates,
      kNumberOfCoordinates) = M_uuc();
  
  m.block<kNumberOfSpeeds - kNumberOfVelocityConstraints, kNumberOfCoordinates>
      (kNumberOfCoordinates + kNumberOfVelocityConstraints, 0) = M_uqd();
  
  m.block<kNumberOfSpeeds - kNumberOfVelocityConstraints, kNumberOfSpeeds>
    (kNumberOfCoordinates + kNumberOfVelocityConstraints, kNumberOfCoordinates) = M_uud();
  return m;
}

Matrix<double, 20, 16, RowMajor> Bicycle::independent_state_matrix() const
{
  Matrix<double, 20, 16, RowMajor> m;
  
  m.block<8, 8>(0, 0) = A_qq();

  return m;
}

bool Bicycle::is_dependent_index(int i) const
{
  return (i == dependent_speeds_[0]) ||
         (i == dependent_speeds_[1]) ||
         (i == dependent_speeds_[2]);
}

} // namespace bicycle
