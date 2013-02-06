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
  double ls_, g_;
  double T_rw, T_s, T_fw, G_Rx, G_Ry, G_Rz, F_Rx, F_Ry, F_Rz,
         T_Rx, T_Ry, T_Rz, G_Fx, G_Fy, G_Fz, F_Fx, F_Fy, F_Fz,
         T_Fx, T_Fy, T_Fz;

  int dependent_coordinate_index, dependent_speeds_indices[3];
};
#endif
