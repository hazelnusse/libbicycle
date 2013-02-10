#include <cmath>
#include <iostream>
#include <iomanip>
#include "bicycle.h"
using namespace Eigen;

Bicycle::Bicycle()
  : state_(state::Zero()), ls_(0.0), g_(9.81), steer_torque_(0.0),
  dependent_coordinate_(2), dependent_speeds_{0, 2, 5}
{
}

void Bicycle::set_state(state & x)
{
  state_ = x;
}
  
Matrix<double, 6, 1> Bicycle::compute_contact_forces() const
{
  Matrix<double, 3, 12, RowMajor> B;
  Matrix<double, 12, 1> gif;
  Matrix<double, 12, 22, RowMajor> gaf_dr_full;
  
  f_v_coefficient(B.data());  // compute velocity constraint coefficient matrix
  gif_ud_zero(gif.data());    // compute coriolis/centripetal inertia forces
  gaf_dr(gaf_dr_full.data()); // compute generalized active force input matrix

  // TODO: modify to allow for external forces, external torques, and internal
  // torques to be accounted for.  The current implementation assumes
  // everything besides contact forces is equal to zero.
  Matrix<double, 12, 7, RowMajor> gaf_dr_reduced;
  gaf_dr_reduced << gaf_dr_full.block<12, 3>(0, 4),  // Rear contact forces
                    gaf_dr_full.block<12, 3>(0, 14), // Front contact forces
                    gaf_dr_full.block<12, 1>(0, 21); // Gravitational term

  Matrix<double, 3, 3, RowMajor> B_d;      // Dependent columns
  Matrix<double, 3, 9, RowMajor> B_i;      // Independent columns
  Matrix<double, 3, 7, RowMajor> gaf_dr_d; // Dependent active forces
  Matrix<double, 9, 7, RowMajor> gaf_dr_i; // Independent active forces
  Matrix<double, 3, 1> gif_ud_zero_d;      // Independent coriolis/centripetal
  Matrix<double, 9, 1> gif_ud_zero_i;      // Dependent coriolis/centripetal

  for (int i = 0, ii = 0, id = 0; i < 12; ++i) {
    if (i == dependent_speeds_[0] ||
        i == dependent_speeds_[1] ||
        i == dependent_speeds_[2]) {
      B_d.block<3, 1>(0, id) = B.block<3, 1>(0, i);
      gaf_dr_d.block<1, 7>(id, 0) = gaf_dr_reduced.block<1, 7>(i, 0);
      gif_ud_zero_d[id] = gif[i];
      ++id;
    } else {
      B_i.block<3, 1>(0, ii) = B.block<3, 1>(0, i);
      gaf_dr_i.block<1, 7>(ii, 0) = gaf_dr_reduced.block<1, 7>(i, 0);
      gif_ud_zero_i[ii] = gif[i];
      ++ii;
    }
  }

  Matrix<double, 9, 3, RowMajor> C;
  C = - (B_d.inverse() * B_i).transpose();

  Matrix<double, 9, 7, RowMajor> gaf_dr_constrained;
  Matrix<double, 9, 1> rhs;
  gaf_dr_constrained = (gaf_dr_i + (C*gaf_dr_d));
  rhs = - (gif_ud_zero_i
         + C*gif_ud_zero_d
         + gaf_dr_constrained.block<9, 1>(0, 6)*g_);

  return gaf_dr_constrained.block<6, 6>(3, 0)
          .fullPivHouseholderQr().solve(rhs.block<6, 1>(3, 0));
}

void Bicycle::set_dependent_coordinate(int i) {
  if (!((i == 1) | (i == 2) | (i == 3))) {
    std::cerr << "Invalid dependent coordinate.  Ignoring request." << std::endl
              << "dependent coordinate must be 1, 2, or 3 (lean, pitch, steer)." << std::endl;
    return;
  }
  dependent_coordinate_ = i;
}

std::ostream & operator<<(std::ostream & os,
                          const Bicycle & b)
{
  os << "Rear assembly:" << std::endl
     << b.rear_
     << "Front assembly:" << std::endl
     << b.front_;
  os.precision(16);
  os << "ls  = " << std::setw(25) << b.ls_ << std::endl
     << "g   = " << std::setw(25) << b.g_  << std::endl
     << "State:" << std::endl 
     << " x[0] = " << std::setw(25) <<  b.state_[0] << "    (yaw)" << std::endl
     << " x[1] = " << std::setw(25) <<  b.state_[1] << "    (lean)" << std::endl
     << " x[2] = " << std::setw(25) <<  b.state_[2] << "    (pitch)" << std::endl
     << " x[3] = " << std::setw(25) <<  b.state_[3] << "    (steer)" << std::endl
     << " x[4] = " << std::setw(25) <<  b.state_[4] << "    (rear wheel angle)" << std::endl
     << " x[5] = " << std::setw(25) <<  b.state_[5] << "    (front wheel angle)" << std::endl
     << " x[6] = " << std::setw(25) <<  b.state_[6] << "    (x of rear wheel contact)" << std::endl
     << " x[7] = " << std::setw(25) <<  b.state_[7] << "    (y of rear wheel contact)" << std::endl
     << " x[8] = " << std::setw(25) <<  b.state_[8] << "    (yaw rate)" << std::endl
     << " x[9] = " << std::setw(25) <<  b.state_[9] << "    (lean rate)" << std::endl
     << "x[10] = " << std::setw(25) << b.state_[10] << "    (pitch rate)" << std::endl
     << "x[11] = " << std::setw(25) << b.state_[11] << "    (steer rate)" << std::endl
     << "x[12] = " << std::setw(25) << b.state_[12] << "    (rear wheel rate)" << std::endl
     << "x[13] = " << std::setw(25) << b.state_[13] << "    (front wheel rate)" << std::endl
     << "x[14] = " << std::setw(25) << b.state_[14] << "    (rear wheel contact longitudinal rate)" << std::endl
     << "x[15] = " << std::setw(25) << b.state_[15] << "    (rear wheel contact lateral rate)" << std::endl
     << "x[16] = " << std::setw(25) << b.state_[16] << "    (rear wheel contact vertical rate)" << std::endl
     << "x[17] = " << std::setw(25) << b.state_[17] << "    (front wheel contact longitudinal rate)" << std::endl
     << "x[18] = " << std::setw(25) << b.state_[18] << "    (front wheel contact lateral rate)" << std::endl
     << "x[19] = " << std::setw(25) << b.state_[19] << "    (front wheel contact vertical rate)" << std::endl;
  return os;
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
