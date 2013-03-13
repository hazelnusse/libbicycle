#include <algorithm>
#include <queue>
#include "bicycle.h"

namespace bicycle {

Vector Bicycle::steady_no_slip_constraint_forces() const
{
  const Matrix C = Bd_inverse_Bi().transpose(); // 9 by 3 constraint matrix

  Vector gif_cc(o);                   // Generalized inertia (coriolis/centripetal)
  gif_ud_zero_steady(gif_cc.data());         // populate gif_cc
  gif_cc = P_u_.transpose() * gif_cc; // reorder rows
  Vector gif_cc_constrained(o - m);
  gif_cc_constrained = gif_cc.block<o - m, 1>(0, 0)    // indepedent rows
                     + C * gif_cc.block<m, 1>(o - m, 0);   //  dependent rows

  Matrix gaf_dr_full(o, s);           // Input coefficient matrix
  gaf_dr(gaf_dr_full.data());         // populate gaf_dr_full
  gaf_dr_full = P_u_.transpose() * gaf_dr_full; // reorder rows
  Matrix gaf_dr_full_constrained(o - m, s); 
  gaf_dr_full_constrained = gaf_dr_full.block<o - m, s>(0, 0)  // independent rows
                      + C * gaf_dr_full.block<m, s>(o - m, 0); //   dependent rows


  Matrix gaf_dr_c_constrained(o - m, 7);     // Columns associated with constraint forces
  gaf_dr_c_constrained << gaf_dr_full_constrained.col(4),    // rear longitudinal force
                          gaf_dr_full_constrained.col(5),    // rear lateral force
                          gaf_dr_full_constrained.col(6),    // rear normal force
                          gaf_dr_full_constrained.col(14),   // front longitudinal force
                          gaf_dr_full_constrained.col(15),   // front lateral force
                          gaf_dr_full_constrained.col(16),   // front normal force
                          gaf_dr_full_constrained.col(20);   // steer torque

  Matrix gaf_dr_a_constrained(o - m, s - 7); // Columns associated with active forces
  gaf_dr_a_constrained << gaf_dr_full_constrained.col(0),    // rear wheel torque
                          gaf_dr_full_constrained.col(1),    // rear x torque
                          gaf_dr_full_constrained.col(2),    // rear y torque
                          gaf_dr_full_constrained.col(3),    // rear z torque
                          gaf_dr_full_constrained.col(7),    // rear x force
                          gaf_dr_full_constrained.col(8),    // rear y force
                          gaf_dr_full_constrained.col(9),    // rear z force
                          gaf_dr_full_constrained.col(10),   // front wheel torque
                          gaf_dr_full_constrained.col(11),   // front x torque
                          gaf_dr_full_constrained.col(12),   // front y torque
                          gaf_dr_full_constrained.col(13),   // front z torque
                          gaf_dr_full_constrained.col(17),   // front x force
                          gaf_dr_full_constrained.col(18),   // front y force
                          gaf_dr_full_constrained.col(19),   // front z force
                          gaf_dr_full_constrained.col(21);   // gravity

  // At this point we have a system of 9 equations with seven unknowns.  The
  // last 6 equations come from the generalized speeds associated with the
  // contact point velocities and need to be used when solving for constraint
  // forces.  However, the first three equations come from the independent
  // speeds which may change depending on parameters or configuration.
  // Therefore, it isn't safe to assume which of the first three equations we
  // need.  Instead of picking one, we can form a Jacobi SVD and do a least
  // squares solution.  This approach should be very numerically robust.
  Vector rhs = -(gif_cc_constrained
      + gaf_dr_a_constrained * all_inputs_except_constraint_forces());

  JacobiSVD<Matrix, ::Eigen::FullPivHouseholderQRPreconditioner>
    svd(gaf_dr_c_constrained, ::Eigen::ComputeFullU | ::Eigen::ComputeFullV);
  return svd.solve(rhs);
}

Vector Bicycle::no_slip_contact_forces() const
{
  // TODO: Implement me.
  Vector cf(6);
  return cf;
}

std::pair<int, double> Bicycle::solve_configuration_constraint_and_set_state(double ftol, int iter) {
  const double df_min = 1e-14;
  Vector df(n);
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
  return std::pair<int, double>(i, f);
}

Vector Bicycle::solve_velocity_constraints_and_set_state()
{
  Vector u_i = (P_u_.transpose() * state_.block<o, 1>(n, 0)).block<o - m, 1>(0, 0);

  Matrix B(m, o);
  f_v_du(B.data());  // compute velocity constraint coefficient matrix
  B = B * P_u_;      // move dependent columns to the end
  FullPivHouseholderQR<Matrix> dec(B.block<m, m>(0, o - m));
  Vector u_d = -dec.solve(B.block<m, o - m>(0, 0) * u_i);

  int i = 0;
  for (auto it = dependent_speeds_.begin();
           it != dependent_speeds_.end(); ++it)
    state_[*it + n] = u_d[i++];
  
  // return the residual (should be zero or nearly zero)
  return B.block<m, m>(0, o - m) * u_d + B.block<m, o - m>(0, 0) * u_i;
}

Matrix Bicycle::Bd_inverse_Bi() const
{
  Matrix B(m, o);
  f_v_du(B.data());  // compute velocity constraint coefficient matrix
  B = B * P_u_;      // move dependent columns to the end

  // Decompose B_d
  FullPivHouseholderQR<Matrix> decomposition(B.block<m, m>(0, o - m));

  // Solve B_d * u_d = - B_i * u_i for (-B_d^-1 * B_i)
  return decomposition.solve(-B.block<m, o - m>(0, 0));
}

Matrix Bicycle::f_v_dq() const
{
  Matrix mat = Matrix::Zero(m, n);
  Vector B_dq_raw(m * o * n_min);
  f_v_dudq(B_dq_raw.data()); // Populate the raw data

  // Perform 3 matrix multiplies of 3 x 12 * 12 x 1
  // Each matrix multiply results in the gradient of f_v with respect to lean,
  // pitch, or steer
  //for (int k = 1; k < n_min + 1; ++k) { // lean, pitch, steer
  //  mat.block<m, 1>(0, k) = Map<Matrix, Unaligned, Stride<m*o, n_min>>
  //    (B_dq_raw.data() + k - 1, m, o) * state_.block<o, 1>(n, 0);
  //}

  // This can also be done by iterating over the 12 speeds and accumulating
  // the product of a 3x3 matrix multiplied by each speed.
  for (int j = 0; j < o; ++j) {
    Matrix tmp = Map<Matrix, Unaligned, Stride<n_min*o, 1>>(B_dq_raw.data() + n_min*j, n_min, n_min);
    tmp *= state_[n + j];
    mat.block<3, 3>(0, 1) += tmp;
  }
  // both approaches give the same result.

  return mat;
}

std::set<int> Bicycle::best_dependent_speeds() const
{
  Matrix B(m, o);
  f_v_du(B.data());
  B = B.block<m, 6>(0, 0);
  JacobiSVD<Matrix> svd(B, ComputeThinV);

  int r = svd.nonzeroSingularValues();
  if (r < m) {
    std::cerr << "Not all constraints are active. Row rank of the constraint "
                 "matrix is " << r << std::endl;
  }
  Matrix R = svd.matrixV();
  Matrix d = R.rowwise().squaredNorm();
  
  std::priority_queue<std::pair<double, int>> q;
  for (int i = 0; i < 6; ++i)
    q.push(std::pair<double, int>(d(i, 0), i));
 
  std::set<int> indices;
  for (int i = 0; i < m; ++i) {
    indices.insert(q.top().second);
    q.pop();
  }

  return indices;
}

int Bicycle::best_dependent_coordinate() const
{
  Matrix df(n, 1);
  f_c_dq(df.data());
  df = df.cwiseAbs();
  return std::distance(df.data(), std::max_element(df.data(), df.data() + n));
}

bool Bicycle::is_dependent_index(int i) const
{
  return dependent_speeds_.count(i);
}

Vector Bicycle::coordinate_derivatives() const
{
  Vector rhs(n);
  f_1(rhs.data());
  return -rhs;
}

Vector Bicycle::speed_derivatives() const
{
  const Matrix C = Bd_inverse_Bi().transpose();
  Matrix mm(o, o);
  Vector forcing(o);

  // Portion of mass matrix associated with acceleration constraints
  Matrix B(m, o);
  f_v_du(B.data());
  mm.block<m, o>(0, 0) = B;

  // Portion of forcing vector associated with acceleration constraints
  Matrix B_dot(m, o);
  f_v_dudt(B_dot.data());
  forcing.block<m, 1>(0, 0) = -B_dot * state_.block<o, 1>(n, 0);
  
  // Portion of mass matrix associated with dynamic equations
  Matrix mm_d(o, o);
  gif_dud(mm_d.data());
  mm_d = P_u_.transpose() * mm_d;  // reorder rows: independent then dependent
  mm.block<o - m, o>(m, 0) = mm_d.block<o - m, o>(0, 0) +
                         C * mm_d.block<m, o>(o - m, 0); 

  // Portion of forcing vector associated with coriolis/centripetal terms from  dynamic equations
  Vector gif_unconstrained(o);
  gif_ud_zero(gif_unconstrained.data());
  gif_unconstrained = P_u_.transpose() * gif_unconstrained; // reorder rows
  forcing.block<o - m, 1>(m, 0) = -(gif_unconstrained.block<o - m, 1>(0, 0) +
                                C * gif_unconstrained.block<m, 1>(o - m, 0));

  // Portion of forcing vector associated with generalized active forces
  Vector gaf_unconstrained(o);
  gaf(gaf_unconstrained.data());
  gaf_unconstrained = P_u_.transpose() * gaf_unconstrained; // reorder rows
  forcing.block<o - m, 1>(m, 0) -= (gaf_unconstrained.block<o - m, 1>(0, 0) +
                                C * gaf_unconstrained.block<m, 1>(o - m, 0));

  return mm.fullPivHouseholderQr().solve(forcing);
}

Vector Bicycle::state_derivatives() const
{
  Vector dxdt(o);
  dxdt << coordinate_derivatives(), speed_derivatives();
  return dxdt;
}

} // namespace bicycle
