#ifndef BICYCLE_H
#define BICYCLE_H
#include <Eigen/Dense>
#include "wheelassemblygyrostat.h"

class Bicycle {
 public:
  typedef Eigen::Matrix<double, 20, 1> state;

  Bicycle();



 private:
  // Generated code, don't expose as this is an implementation detail
  void f_c(double m[1]);
  void f_c_dq(double m[8]);
  void f_v_coefficient(double m[36]);
  void f_v_coefficient_dq(double m[108]);
  void f_v_coefficient_dqdq(double m[324]);
  void gif_dud(double m[]);
  void gif_ud_zero(double m[]);
  void gif_ud_zero_dqdu(double m[]);
  void gaf(double m[]);
  void gaf_dqdr(double m[]);

  static const int kNumberOfCoordinates = 8;
  static const int kNumberOfSpeeds = 12;
  static const int kNumberOfVelocityConstraints = 3;
  WheelAssemblyGyrostat rear_, front_;
  double x_[20];
  double ls_, g_, T_s;
  int dependent_coordinate_index, dependent_speeds_indices[3];
};
#endif
