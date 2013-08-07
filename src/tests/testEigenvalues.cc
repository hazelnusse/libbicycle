#include <cmath>
#include <complex>
#include <set>
#include <Eigen/Dense>
#include "gtest/gtest.h"
#include "benchmark_eigenvalues.h"
#include "benchmark_steady_turns.h"
#include "bicycle.h"
#include "wheelassemblygyrostat.h"
#include "whipple.h"

using namespace Eigen;

namespace std {
bool operator<(std::complex<double> l, std::complex<double> r)
{
  double l_norm = ::std::norm(l);
  double r_norm= ::std::norm(r);
  if (l_norm != r_norm) { // they have different magnitudes
    return l_norm < r_norm;
  } else {                  // magnitudes are equal
    if (l == r) {    // they are exactly equal
      return false;
    } else {         // they have equal magnitude, but are not equal, sort by phase angle
      // atan2 returns a number [-pi, +pi]
      double l_theta = ::std::atan2(l.imag(), l.real()),
             r_theta = ::std::atan2(r.imag(), r.imag());
      return l_theta < r_theta;
    }
  }
}
}

TEST(UprightSteadyForwardCruise, BenchmarkEigenvalues)
{
  bicycle::Bicycle b;
  bicycle::Whipple w;
  b.set_parameters(w);
  b.solve_configuration_constraint_and_set_state();
  bicycle::Vector cf = b.steady_no_slip_constraint_forces();
  VectorXd r = VectorXd::Zero(22);
  r[4] = cf[0];
  r[5] = cf[1];
  r[6] = cf[2];
  r[14] = cf[3];
  r[15] = cf[4];
  r[16] = cf[5];
  r[20] = cf[6];
  r[21] = 9.81;
  b.set_inputs(r);
  
  bicycle::Matrix A_min(4, 4);
  A_min << 0, 0, 1, 0,
           0, 0, 0, 1,
           0, 0, 0, 0,
           0, 0, 0, 0;

  Map<Eigen::Matrix<std::complex<double>,11, 4, ColMajor>>
        e(reinterpret_cast<std::complex<double> *>(evals));

  for (int i = 0; i < 11; ++i) {
    b.set_speed(4, -i/w.rR);
    b.solve_velocity_constraints_and_set_state();
    bicycle::Matrix A_full = b.mass_matrix_full()
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

    MatrixXcd ei_calculated = A_min.eigenvalues().transpose();
    Matrix<std::complex<double>, 1, 4> ei = e.row(i);
    std::set<std::complex<double>> unmatched_eigenvalues;
    for (int j = 0; j < 4; ++j)
      unmatched_eigenvalues.insert(ei(0, j));

    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        if (std::abs(ei(0, j) - ei_calculated(0, k)) < 1e-12)
          unmatched_eigenvalues.erase(ei(0, j));
      }
    }
    EXPECT_TRUE(unmatched_eigenvalues.empty());
  }
}

TEST(HandsFreeTurns, BenchmarkEigenvalues)
{
  // TODO: Finish implementing checks of all Basu-Mandal steady turning
  // results.
  bicycle::Bicycle b;
  bicycle::Whipple w;
  b.set_parameters(w);
  b.set_coordinates_basu_mandal(Map<Matrix<double, 9, 1>>(q));
  b.solve_configuration_constraint_and_set_state();
  b.set_dependent_speeds({0, 2, 4});  // yaw, pitch, rear wheel rates dependent
  b.set_speeds_basu_mandal(Map<Matrix<double, 9, 1>>(q_dot));
  b.solve_velocity_constraints_and_set_state();
  bicycle::Vector cf = b.steady_no_slip_constraint_forces();
  VectorXd r = VectorXd::Zero(22);
  r[4] = cf[0];
  r[5] = cf[1];
  r[6] = cf[2];
  r[14] = cf[3];
  r[15] = cf[4];
  r[16] = cf[5];
  r[20] = cf[6];
  r[21] = 9.81;
  b.set_inputs(r);

}
