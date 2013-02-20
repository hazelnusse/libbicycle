#include <cmath>
#include "bicycle.h"

namespace bicycle {

Matrix Bicycle::points_of_interest() const
{
  Matrix mat(7, 3);
  rear_wheel_center_point(mat.data());
  rear_mass_center_point(mat.data() + 3);
  rear_steer_axis_point(mat.data() + 6);
  front_wheel_center_point(mat.data() + 9);
  front_mass_center_point(mat.data() + 12);
  front_steer_axis_point(mat.data() + 15);
  front_ground_contact_point(mat.data() + 18);
  return mat;
}

} // namespace bicycle
