#include <algorithm>
#include <queue>
#include "bicycle.h"

namespace bicycle {

Matrix Bicycle::steady_constraint_forces() const
{
 // indices of ground contact forces and steer torque
  const std::set<int> cf_indices = {4, 5, 6, 14, 15, 16, 20};
  Matrix gif_steady(o, 1),    // Generalized inertia forces
                 gaf_dr_full(o, s);  // Input coefficient matrix
  gif_ud_zero(gif_steady.data());      // populate gif_steady
  gaf_dr(gaf_dr_full.data());          // populate gaf_dr_full
  // put rows associated with dependent speeds at the bottom
  gif_steady = P_u_.transpose() * gif_steady;
  gaf_dr_full = P_u_.transpose() * gaf_dr_full;
  
  // constraint coefficient matrix associated with auxilliary speeds
  Matrix C_aux(7, 3);
  {
    Matrix tmp = Bd_inverse_Bi().transpose(); // 9 by 3
    C_aux.block<1, m>(0, 0) = tmp.block<1, m>(1, 0); //  steer rate
    C_aux.block<6, m>(1, 0) = tmp.block<6, m>(3, 0); //  contact point rates
  }

  ::Eigen::Matrix<int, s, 1> indices;
  std::iota(indices.data(), indices.data() + s, 0); // 0,1,...,21
  std::stable_partition(indices.data(), indices.data() + s,
      [&cf_indices](int elem) { return !cf_indices.count(elem); });
  // Permute columns so constraint forces are at end
  gaf_dr_full = gaf_dr_full * PermutationMatrix<s>(indices);

  Matrix A(7, 7);
  A.block<6, 7>(0, 0) = gaf_dr_full.block<6, 7>(3, 15)
                + C_aux.block<6, 3>(1, 0) * gaf_dr_full.block<3, 7>(9, 15);
  A.block<1, 7>(6, 0) = gaf_dr_full.block<1, 7>(1, 15)
                + C_aux.block<1, 3>(0, 0) * gaf_dr_full.block<3, 7>(9, 15);

  Matrix b(7, 1);
  // Inertia terms
  b.block<6, 1>(0, 0) = - (gif_steady.block<6, 1>(3, 0)
                + C_aux.block<6, 3>(1, 0) * gif_steady.block<3, 1>(9, 0));
  b.block<1, 1>(6, 0) = - (gif_steady.block<1, 1>(1, 0)
                + C_aux.block<1, 3>(0, 0) * gif_steady.block<3, 1>(9, 0));

  b.block<6, 1>(0, 0) -= (gaf_dr_full.block<6, 15>(3, 0)
                + C_aux.block<6, 3>(1, 0) * gaf_dr_full.block<3, 15>(9, 0)) * all_inputs_except_constraint_forces();
  b.block<1, 1>(6, 0) -= (gaf_dr_full.block<1, 15>(1, 0)
                + C_aux.block<1, 3>(0, 0) * gaf_dr_full.block<3, 15>(9, 0)) * all_inputs_except_constraint_forces();

  // TODO: Add term involving mass matrix and du/dt for non-steady case

  return A.fullPivHouseholderQr().solve(b);
}

void Bicycle::solve_configuration_constraint_and_set_state(double ftol, int iter) {
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
}

void Bicycle::solve_velocity_constraints_and_set_state()
{
  Vector u_independent(o - m);
  u_independent = (state_.block<o, 1>(n, 0).transpose() * P_u_).block<1, o - m>(0, 0).transpose();
  Vector u_dependent = Bd_inverse_Bi() * u_independent;
  int i = 0;
  for (auto it = dependent_speeds_.begin();
       it != dependent_speeds_.end(); ++it) {
    state_[*it + n] = u_dependent(i++, 0);
  }
}

Matrix Bicycle::Bd_inverse_Bi() const
{
  Matrix B(m, o);
  f_v_du(B.data());  // compute velocity constraint coefficient matrix
  B = B * P_u_;      // move dependent columns to the end

  // Decompose B_d
  FullPivHouseholderQR<Matrix> decomposition;
  decomposition.compute(B.block<m, m>(0, o - m));

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
  // This could also be done by iterting over the 12 speeds and accumulating
  // the product of a 3x3 matrix multiplied by each speed.
  for (int k = 1; k < n_min + 1; ++k) { // lean, pitch, steer
    mat.block<m, 1>(0, k) = Map<Matrix, Unaligned, Stride<m*o, n_min>>
      (B_dq_raw.data() + k - 1, m, o) * state_.block<o, 1>(n, 0);
  }

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

} // namespace bicycle
