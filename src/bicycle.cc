#include "bicycle.h"

Bicycle::Bicycle()
  : rear_{}, front_{}, x_{}, ls_(0.0), g_(9.81), T_s(0.0),
  dependent_coordinate_index(2), dependent_speeds_indices{0, 2, 5}
{
}
