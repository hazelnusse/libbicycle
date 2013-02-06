#include "bicycle.h"

Bicycle::Bicycle()
  : rear_{}, front_{}, x_{}, ls_(0.0), g_(9.81),
  T_rw(0.0),  T_s(0.0), T_fw(0.0), G_Rx(0.0), G_Ry(0.0), G_Rz(0.0), F_Rx(0.0),
  F_Ry(0.0), F_Rz(0.0), T_Rx(0.0), T_Ry(0.0), T_Rz(0.0), G_Fx(0.0), G_Fy(0.0),
  G_Fz(0.0), F_Fx(0.0), F_Fy(0.0), F_Fz(0.0), T_Fx(0.0), T_Fy(0.0), T_Fz(0.0),
  dependent_coordinate_index(2), dependent_speeds_indices{0, 2, 5}
{
}
