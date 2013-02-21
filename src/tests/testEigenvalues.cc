#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include "gtest/gtest.h"
#include "benchmark_eigenvalues.h"
#include "bicycle.h"
#include "wheelassemblygyrostat.h"
#include "whipple.h"

using namespace Eigen;

bool complex_less_than(std::complex<double> l, std::complex<double> r)
{
  double l_mag_2 = l.real()*l.real() + l.imag()*l.imag();
  double r_mag_2 = r.real()*r.real() + r.imag()*r.imag();
  if (l_mag_2 != r_mag_2) { // they have different magnitudes
    return l_mag_2 < r_mag_2;
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

TEST(UprightSteadyForwardCruise, BenchmarkEigenvalues)
{
  bicycle::Bicycle b;
  bicycle::Whipple w;
  b.set_parameters(w);
  b.solve_configuration_constraint_and_set_state();
  bicycle::Vector cf = b.steady_constraint_forces();
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
    std::sort(ei_calculated.data(), ei_calculated.data() + 4, complex_less_than);
    Matrix<std::complex<double>, 1, 4> ei = e.row(i);
    std::sort(ei.data(), ei.data() + 4, complex_less_than);
    for (int j = 0; j < 4; ++j) {
      EXPECT_NEAR(ei(0, j).real(), ei_calculated(0, j).real(), 1e-12);
      EXPECT_NEAR(ei(0, j).imag(), ei_calculated(0, j).imag(), 1e-12);
    }
  }
}
