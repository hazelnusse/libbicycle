#include "bicycle.h"

namespace bicycle {

RowMajorMatrix Bicycle::M_qq() const
{
  // \nabla_q f_0
  // This is true for all conditions of bicycle motion
  return RowMajorMatrix::Identity(n, n);
}

RowMajorMatrix Bicycle::M_uqc() const
{
  // \nabla_{dq/dt} f_a
  // This result is true even for du/dt != 0
  return f_v_dq();
}

RowMajorMatrix Bicycle::M_uuc() const
{
  // \nabla_{du/dt} f_a
  // This is true assuming f_a = d/dt (f_v) = d/dt (B(q) * u)
  // This is true even when du/dt != 0
  RowMajorMatrix mat(m, o);
  f_v_coefficient(mat.data());
  return mat;
}

RowMajorMatrix Bicycle::M_uqd() const
{
  // \nabla_{dq/dt} f_3
  // This is true even when du/dt != 0
  return RowMajorMatrix::Zero(o - m, n);
}

RowMajorMatrix Bicycle::M_uud() const
{
  // \nabla_{du/dt} f_2
  // This is the constrained du/dt coefficient matrix from Kane's equations
  // This is true even when du/dt != 0
  RowMajorMatrix mat(o, o);
  gif_dud(mat.data());
  mat = P_u_.transpose() * mat; // reorder the rows
  return mat.block<o - m, o>(0, 0)
       + Bd_inverse_Bi().transpose() * mat.block<m, o>(o - m, 0);
}

RowMajorMatrix Bicycle::A_qq() const
{
  // - \nabla_{q} (f_0 + f_1)
  // Since f_0 = I_{8x8} * \dot{q}, only the f_1 term appears in the gradient
  // This is true even when du/dt != 0
  RowMajorMatrix mat(n, n);
  f_1_dq(mat.data());
  return -mat;
}

RowMajorMatrix Bicycle::A_qu() const
{
  // - \nabla_u f_1
  // This is true even when du/dt != 0
  RowMajorMatrix mat(n, o);
  f_1_du(mat.data());
  return -mat;
}
  
RowMajorMatrix Bicycle::A_uqc() const
{
  // - \nabla_q f_a
  // TODO:  Rederive without assuming lean rate, pitch rate, steer rate and
  // du/dt = 0
  // Assumption:  lean rate, pitch rate, steer rate are zero and du/dt = 0
  return RowMajorMatrix::Zero(m, n);
}

RowMajorMatrix Bicycle::A_uuc() const
{
  // - \nabla_u f_a
  // This is true even when du/dt != 0
  VectorXd B_dq(m*o*n_min, 1);
  f_v_coefficient_dq(B_dq.data()); // Populate the raw data
  RowMajorMatrix B_j_dq = RowMajorMatrix::Zero(m, o);
  // lean rate, pitch rate, steer rate
  // note negative sign!!!
  Vector3d u_tilde = -state_.block<n_min, 1>(n + 1, 0);

  // Iterate over each column of B(q)
  for (int j = 0; j < o; ++j) { // column index
     B_j_dq.block<m, 1>(0, j) = Map<RowMajorMatrix, Unaligned,
                  Stride<n_min*o, 1>>(&B_dq[n_min*j], n_min, n_min) * u_tilde;
  }
  return B_j_dq;
   
  // Result when lean rate, pitch rate, and steer rate are all zero:
  // return RowMajorMatrix::Zero(m, o);
}

RowMajorMatrix Bicycle::A_uqd() const
{
  // - \nabla_q (f_2 + f_3)
  // Assume du/dt = 0, then \nabla_q(f_2) = 0
  // TODO: Rederive without assuming du/dt = 0
  RowMajorMatrix mat(o - m, n);

  // Generalized active forces partial derivatives
  RowMajorMatrix m_gaf_dq(o, n);
  gaf_dq(m_gaf_dq.data());
  m_gaf_dq = P_u_.transpose() * m_gaf_dq; // reorder rows

  // Generalized inertia forces, coriolis/centripetal
  RowMajorMatrix m_gif_ud_zero_dq(o, n);
  gif_ud_zero_dq(m_gif_ud_zero_dq.data());
  m_gif_ud_zero_dq = P_u_.transpose() * m_gif_ud_zero_dq; // reorder rows

  // Constraint coefficient matrix
  RowMajorMatrix B(m, o);
  f_v_coefficient(B.data());
  B = B * P_u_; // reorder columns so dependents speeds are at end.
  RowMajorMatrix B_i(B.block<m, o - m>(0, 0));
  // compute decomposition of B_d
  RowMajorMatrix B_d(B.block<m, m>(0, o - m));
  FullPivHouseholderQR<RowMajorMatrix> decomposition;
  decomposition.compute(B_d);

  // Independent portion (Term 1)
  mat = -(m_gaf_dq.block<o - m, n>(0, 0) + m_gif_ud_zero_dq.block<o - m, n>(0, 0));

  // Dependent portion (Term 2)
  //
  // Term 2 a:  \nabla_q[C] * (F_d + F_{dqu}^*)
  for (int k = 1; k < 4; ++k) {
    mat.block<o - m, 1>(0, k);

  }
  // Term 2 b:   C * \nabla_q[F_d + F_{dqu}^*]
  RowMajorMatrix C = decomposition.solve(B_i).transpose();
  mat.block<o - m, n>(0, 0) -= C*(m_gaf_dq.block<m, n>(o - m, 0) +
                                  m_gif_ud_zero_dq.block<m, n>(o - m, 0));

  return mat;
}

RowMajorMatrix Bicycle::A_uud() const
{
  // - \nabla_u (f_3)
  RowMajorMatrix mat = RowMajorMatrix::Zero(o - m, o);

  return mat;
}

RowMajorMatrix Bicycle::mass_matrix_full() const
{
  RowMajorMatrix mat(n + o, n + o);
  mat.block<n, n>(0, 0) = M_qq();

  mat.block<n, o>(0, n) = RowMajorMatrix::Zero(n, o);

  mat.block<m, n>(n, 0) = M_uqc();
  
  mat.block<m, o>(n, n) = M_uuc();
  
  mat.block<o - m, n>(n + m, 0) = M_uqd();
  
  mat.block<o - m, o>(n + m, n) = M_uud();
  return mat;
}

RowMajorMatrix Bicycle::independent_state_matrix() const
{
  RowMajorMatrix mat(n + o, n + o - l - m);
  
  mat.block<8, 8>(0, 0) = A_qq() * C_0();

  return mat;
}
  
RowMajorMatrix Bicycle::C_0() const
{
  RowMajorMatrix fcdq(l, n);
  f_c_dq(fcdq.data());

  return RowMajorMatrix::Identity(8, 8) - ((P_qd() / (fcdq * P_qd())(0, 0)) * fcdq);
}

RowMajorMatrix Bicycle::C_1() const
{
  RowMajorMatrix B(m, o);
  f_v_coefficient(B.data());
  B = B * P_u_;
  FullPivHouseholderQR<RowMajorMatrix> decomposition;
  decomposition.compute(B.block<m, m>(0, o - m));

  // Constraint coefficient matrix partials w.r.t. q
  //RowMajorMatrix 

}

RowMajorMatrix Bicycle::P_qd() const {
  return P_q_.block<n, l>(0, n - l);
}

RowMajorMatrix Bicycle::P_qi() const {
  return P_q_.block<n, n - l>(0, 0);
}

RowMajorMatrix Bicycle::P_ud() const {
  return P_u_.block<o, m>(0, o - m);
}

RowMajorMatrix Bicycle::P_ui() const {
  return P_u_.block<o, o - m>(0, 0);
}

} // namespace bicycle
