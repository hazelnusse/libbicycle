#ifndef BICYCLE_PRIV_H
#define BICYCLE_PRIV_H

namespace bicycle {

inline void Bicycle::set_state(int i, double xi)
{
  state_[i] = xi;
}

inline void Bicycle::set_coordinate(int i, double qi)
{
  state_[i] = qi;
}

inline void Bicycle::set_speed(int i, double ui)
{
  state_[n + i] = ui;
}

} // namespace bicycle

#endif
