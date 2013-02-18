#include "bicycle.h"

namespace bicycle {

Matrix Bicycle::M_qq() const
{
  // \nabla_q f_0
  // This is true for all conditions of bicycle motion
  return Matrix::Identity(n, n);
}

Matrix Bicycle::M_uqc() const
{
  // \nabla_{dq/dt} f_a
  // This result is true even for du/dt != 0
  return f_v_dq();
}

Matrix Bicycle::M_uuc() const
{
  // \nabla_{du/dt} f_a
  // This is true assuming f_a = d/dt (f_v) = d/dt (B(q) * u)
  // This is true even when du/dt != 0
  Matrix mat(m, o);
  f_v_du(mat.data());
  return mat;
}

Matrix Bicycle::M_uqd() const
{
  // \nabla_{dq/dt} f_3
  // This is true even when du/dt != 0
  return Matrix::Zero(o - m, n);
}

Matrix Bicycle::M_uud() const
{
  // \nabla_{du/dt} f_2
  // This is the constrained du/dt coefficient matrix from Kane's equations
  // This is true even when du/dt != 0
  Matrix mat(o, o);
  gif_dud(mat.data());
  mat = P_u_.transpose() * mat; // reorder the rows
  return mat.block<o - m, o>(0, 0)
       + Bd_inverse_Bi().transpose() * mat.block<m, o>(o - m, 0);
}

Matrix Bicycle::A_qq() const
{
  // - \nabla_{q} (f_0 + f_1)
  // Since f_0 = I_{8x8} * \dot{q}, only the f_1 term appears in the gradient
  // This is true even when du/dt != 0
  Matrix mat(n, n);
  f_1_dq(mat.data());
  return -mat;
}

Matrix Bicycle::A_qu() const
{
  // - \nabla_u f_1
  // This is true even when du/dt != 0
  Matrix mat(n, o);
  f_1_du(mat.data());
  return -mat;
}
  
Matrix Bicycle::A_uqc() const
{
  // - \nabla_q f_a
  // TODO:  Rederive without assuming lean rate, pitch rate, steer rate and
  // du/dt = 0
  // Assumption:  lean rate, pitch rate, steer rate are zero and du/dt = 0
  return Matrix::Zero(m, n);
}

Matrix Bicycle::A_uuc() const
{
  // - \nabla_u f_a
  // This is true even when du/dt != 0
  Vector B_dq(m * o * n_min, 1);
  f_v_dudq(B_dq.data()); // Populate the raw data
  Matrix B_j_dq = Matrix::Zero(m, o);
  // lean rate, pitch rate, steer rate
  // note negative sign!!!
  ::Eigen::Vector3d u_tilde = -state_.block<n_min, 1>(n + 1, 0);

  // Iterate over each column of B(q)
  for (int j = 0; j < o; ++j) { // column index
     B_j_dq.block<m, 1>(0, j) = Map<Matrix, Unaligned,
                  Stride<n_min*o, 1>>(&B_dq[n_min*j], n_min, n_min) * u_tilde;
  }
  return B_j_dq;
}

Matrix Bicycle::A_uqd() const
{
  // - \nabla_q (f_2 + f_3)
  // Assume du/dt = 0, then \nabla_q(f_2) = 0
  // TODO: Rederive without assuming du/dt = 0
  Matrix mat = Matrix::Zero(o - m, n);

  // Generalized active forces
  Vector gaf_vec(o);
  gaf(gaf_vec.data());
  gaf_vec = P_u_.transpose() * gaf_vec;

  // Generalized active forces partial derivatives
  Matrix gaf_dq_mat(o, n_min);
  gaf_dq(gaf_dq_mat.data());
  gaf_dq_mat = P_u_.transpose() * gaf_dq_mat; // reorder rows
  
  // Generalized inertia forces, coriolis/centripetal
  Vector gif_ud_zero_vec(o);
  gif_ud_zero(gif_ud_zero_vec.data());
  gif_ud_zero_vec = P_u_.transpose() * gif_ud_zero_vec; // reorder rows

  // Generalized inertia forces, coriolis/centripetal, partial derivatives
  Matrix gif_ud_zero_dq_mat(o, n_min);
  gif_ud_zero_dq(gif_ud_zero_dq_mat.data());
  gif_ud_zero_dq_mat = P_u_.transpose() * gif_ud_zero_dq_mat; // reorder rows

  // Constraint coefficient matrix
  Matrix B(m, o);               // Constraint coefficient matrix
  f_v_du(B.data());             // Populate matrix
  B = B * P_u_;                 // reorder columns
  Matrix B_i(B.block<m, o - m>(0, 0)); // independent columns of B
  Matrix B_d(B.block<m, m>(0, o - m)); // dependent columns of B
  FullPivHouseholderQR<Matrix> decomposition;  // QR decomposition so we can
  decomposition.compute(B_d);                  // solve B_d * x = rhs for x

  // Independent portion (Term 1)
  mat.block<o - m, n_min>(0, 1) = -(gaf_dq_mat.block<o - m, n_min>(0, 0)
        + gif_ud_zero_dq_mat.block<o - m, n_min>(0, 0));

  // Dependent portion (Term 2)
  //
  // Term 2 a:  \nabla_q[C] * (F_d + F_{dqu}^*)
  //for (int k = 1; k < 4; ++k) {
 //   mat.block<o - m, 1>(0, k);
//
 // }
  // Term 2 b:   C * \nabla_q[F_d + F_{dqu}^*]
  //Matrix C = decomposition.solve(B_i).transpose();
  //mat.block<o - m, n>(0, 0) -= C*(m_gaf_dq.block<m, n>(o - m, 0) +
  //                                m_gif_ud_zero_dq.block<m, n>(o - m, 0));
//
  return mat;
}

Matrix Bicycle::A_uud() const
{
  // - \nabla_u (f_3)
  Matrix mat(o, o);
  gif_ud_zero_du(mat.data());
  mat = P_u_.transpose() * mat; // reorder rows

  return -(mat.block<o - m, o>(0, 0)
           + Bd_inverse_Bi().transpose() * mat.block<m, o>(o - m, 0));
}


Matrix Bicycle::B_u() const
{
  // - \nabla_r (f_3)
  // True when du/dt != 0
  Matrix mat(o, s);
  gaf_dr(mat.data());
  mat = P_u_.transpose() * mat; // reorder rows

  return -(mat.block<o - m, s>(0, 0)
      + Bd_inverse_Bi().transpose() * mat.block<m, s>(o - m, 0));
}

Matrix Bicycle::mass_matrix_full() const
{
  Matrix mat(n + o, n + o);
  mat.block<n, n>(0, 0) = M_qq();

  mat.block<n, o>(0, n) = Matrix::Zero(n, o);

  mat.block<m, n>(n, 0) = M_uqc();
  
  mat.block<m, o>(n, n) = M_uuc();
  
  mat.block<o - m, n>(n + m, 0) = M_uqd();
  
  mat.block<o - m, o>(n + m, n) = M_uud();
  return mat;
}

Matrix Bicycle::independent_state_matrix() const
{
  Matrix mat(n + o, n + o - l - m);
  
  mat.block<n, n - l>(0, 0) = (A_qq() + A_qu() * C_1()) * C_0();

  mat.block<n, o - m>(0, n - l) = A_qu() * C_2();

  mat.block<m, n - l>(n - l, 0) = (A_uqc() + A_uuc()*C_1()) * C_0();
  
  mat.block<m, o - m>(n - l, n - l) = A_uuc() * C_2();

  mat.block<o - m, n - l>(n + m, 0) = (A_uqd() + A_uud() * C_1()) * C_0();

  mat.block<o - m, o - m>(n + m, n - l) = A_uud() * C_2();

  return mat;
}
  
Matrix Bicycle::C_0() const
{
  // [I_{n x n} - P_{qd}*(\nabla_q f_c * P_{qd})^{-1} \nabla_q f_c] * P_{qi}
  // shape: n x (n - l)
  Matrix fcdq(l, n);
  f_c_dq(fcdq.data());

  return (Matrix::Identity(8, 8)
         - ((P_qd() / (fcdq * P_qd())(0, 0)) * fcdq)) * P_qi();
}

Matrix Bicycle::C_1() const
{
  // - P_{ud} * (\nabla_u f_v * P_{ud})^{-1} * \nabla_q f_v
  // shape o x n
  Matrix B(m, o);
  f_v_du(B.data());
  FullPivHouseholderQR<Matrix> decomposition;
  decomposition.compute((B * P_u_).block<m, m>(0, o - m));
  return - (P_ud() * decomposition.solve(f_v_dq()));
}

Matrix Bicycle::C_2() const
{
  // [I_{o x o} - P_{ud} * (\nabla_u f_v)^{-1} * \nabla_u f_v] * P_{ui}
  // shape: o x (o - m)
  Matrix B(m, o);
  f_v_du(B.data());
  FullPivHouseholderQR<Matrix> decomposition;
  decomposition.compute((B * P_u_).block<m, m>(0, o - m));
  return (Matrix::Identity(o, o) - (P_ud() * decomposition.solve(B)))
    * P_ui();
}

Matrix Bicycle::P_qd() const {
  return P_q_.block<n, l>(0, n - l);
}

Matrix Bicycle::P_qi() const {
  return P_q_.block<n, n - l>(0, 0);
}

Matrix Bicycle::P_ud() const {
  return P_u_.block<o, m>(0, o - m);
}

Matrix Bicycle::P_ui() const {
  return P_u_.block<o, o - m>(0, 0);
}

} // namespace bicycle
