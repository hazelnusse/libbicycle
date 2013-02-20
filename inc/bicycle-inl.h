#ifndef BICYCLE_PRIV_H
#define BICYCLE_PRIV_H

namespace bicycle {

inline void Bicycle::set_state(int i, double xi)
{
  state_[i] = xi;
}

inline Vector Bicycle::state() const
{
  return state_;
}

inline void Bicycle::set_coordinate(int i, double qi)
{
  state_[i] = qi;
}

inline void Bicycle::set_speed(int i, double ui)
{
  state_[n + i] = ui;
}

inline WheelAssemblyGyrostat Bicycle::rear_parameters() const
{
  return rear_;
}

inline WheelAssemblyGyrostat Bicycle::front_parameters() const
{
  return front_;
}

inline double Bicycle::steer_axis_offset() const
{
  return ls_;
}

} // namespace bicycle

#endif

