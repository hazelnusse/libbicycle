#ifndef BICYCLE_H
#define BICYCLE_H
#include "wheelassemblygyrostat.h"

class Bicycle {
 public:
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
  void gaf_dq(double m[]);
  void gaf_dr(double m[]);

  WheelAssemblyGyrostat rear_, front_;
  double x_[8];
  double ls_, g_;
  int q_dep_index, u_dep_index[3];
};
#endif
