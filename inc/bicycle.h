#ifndef BICYCLE_H
#define BICYCLE_H
#include <iostream>
#include <Eigen/Dense>
#include "wheelassemblygyrostat.h"

class Whipple;

class Bicycle {
 public:
  typedef Eigen::Matrix<double, 20, 1> state;

  Bicycle();

  void set_state(state & x);
  void set_parameters_from_whipple(const Whipple & w);

  friend std::ostream & operator<<(std::ostream & os,
                                   const Bicycle & b);

  void set_rear_contact_forces() {}
  void set_rear_frame_forces() {}
  void set_rear_frame_torques() {}
  void set_front_contact_forces() {}
  void set_front_frame_forces() {}
  void set_front_frame_torques() {}
  void set_steer_torque() {}

  Eigen::Matrix<double, 6, 1> compute_contact_forces() const;


  // This is to ensure state has 128-bit alignment and hence vectorizable
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  // Generated code, don't expose as this is an implementation detail
  void f_c(double m[1]) const;
  void f_c_dq(double m[8]) const;
  void f_v_coefficient(double m[36]) const;
  void f_v_coefficient_dq(double m[108]) const;
  void f_v_coefficient_dqdq(double m[324]) const;
  void kinematic_odes_rhs(double m[8]) const;
  void gif_dud(double m[144]) const;
  void gif_ud_zero(double m[12]) const;
  void gif_ud_zero_dqdu(double m[240]);
  void gaf(double m[12]) const;
  void gaf_dq(double m[96]) const;
  void gaf_dr(double m[264]) const;
  void gaf_dqdr(double m[360]) const;

  static const int kNumberOfCoordinates = 8;
  static const int kNumberOfSpeeds = 12;
  static const int kNumberOfVelocityConstraints = 3;
  state state_;
  WheelAssemblyGyrostat rear_, front_;
  double ls_, g_, steer_torque_;
  int dependent_coordinate_, dependent_speeds_[3];
};
#endif
