#ifndef BICYCLE_H
#define BICYCLE_H
#include <iostream>
#include <Eigen/Dense>
#include "wheelassemblygyrostat.h"

class Whipple;

class Bicycle {
 public:
  typedef Eigen::Matrix<double, 8, 1> coordinates;
  typedef Eigen::Matrix<double, 12, 1> speeds;
  typedef Eigen::Matrix<double, 20, 1> state;

  Bicycle();

  void set_state(const state & x);
  void set_state(int i, double xi);

  void set_coordinates(const coordinates & q);
  void set_coordinate(int i, double qi);

  void set_speeds(const speeds & u);
  void set_speed(int i, double ui);

  void set_parameters(const WheelAssemblyGyrostat & rear,
                      const WheelAssemblyGyrostat & front,
                      double ls, double g);
  void set_parameters_from_whipple(const Whipple & w);

  void set_dependent_coordinate(int i);
  void set_dependent_speeds(int dependent_speed_indices[3]);

  friend std::ostream & operator<<(std::ostream & os,
                                   const Bicycle & b);

  Eigen::Matrix<double, 6, 1> compute_contact_forces() const;
  Eigen::Matrix <double, 3, 9> compute_Bd_inverse_Bi() const;
  void solve_configuration_constraint_and_set_state(double ftol=1e-14,
                                                    int iter=20);
  void solve_velocity_constraints_and_set_state();



  // This is to ensure state has 128-bit alignment and hence vectorizable
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  // Generated code, don't expose as this is an implementation detail
  void gc_r_ogl(double m[16]) const;
  void wc_r_ogl(double m[16]) const;
  void mc_r_ogl(double m[16]) const;
  void gc_f_ogl(double m[16]) const;
  void wc_f_ogl(double m[16]) const;
  void mc_f_ogl(double m[16]) const;
  void N_ogl(double m[16]) const;
  void f_c(double m[1]) const;
  void f_c_dq(double m[8]) const;
  void f_v_coefficient(double m[36]) const;
  void f_v_coefficient_dq(double m[108]) const;
  void f_v_coefficient_dqdq(double m[324]) const;
  void kinematic_odes_rhs(double m[8]) const;
  void gif_dud(double m[144]) const;
  void gif_ud_zero(double m[12]) const;
  void gif_ud_zero_dqdu(double m[240]) const;
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
  double azimuth, elevation, twist, cam_x, cam_y, cam_z;
  int dependent_coordinate_, dependent_speeds_[3];
};

#include "bicycle_priv.h"
#endif
