#include "bicycle.h"

namespace bicycle {

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
  // - \nabla_q (f_2 + f_3)
  // If we assume du/dt = 0, then -\nabla_q(f_2) = 0
  // TODO: Rederive without assuming du/dt = 0
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
  
RowMajorMatrix Bicycle::C_0() const
{
  RowMajorMatrix fcdq(1, 8);
  f_c_dq(fcdq.data());

  return RowMajorMatrix::Identity(8, 8)
    - ((P_qd() / (fcdq * P_qd())(0, 0)) * fcdq);
}

RowMajorMatrix Bicycle::P_qd() const {
  return P_q_.block<kNumberOfCoordinates, kNumberOfConfigurationConstraints>
                (0, kNumberOfCoordinates - kNumberOfConfigurationConstraints);
}

RowMajorMatrix Bicycle::P_qi() const {
  return P_q_.block<kNumberOfCoordinates,
         kNumberOfCoordinates - kNumberOfConfigurationConstraints>(0, 0);
}

RowMajorMatrix Bicycle::P_ud() const {
  return P_u_.block<kNumberOfSpeeds, kNumberOfVelocityConstraints>
                (0, kNumberOfSpeeds - kNumberOfVelocityConstraints);
}

RowMajorMatrix Bicycle::P_ui() const {
  return P_u_.block<kNumberOfSpeeds,
         kNumberOfSpeeds - kNumberOfVelocityConstraints>(0, 0);
}

} // namespace bicycle
