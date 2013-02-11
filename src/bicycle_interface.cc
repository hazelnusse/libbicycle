#include <iostream>
#include <iomanip>
#include "bicycle.h"

Bicycle::Bicycle()
  : state_(state::Zero()), ls_(0.0), g_(9.81), steer_torque_(0.0),
  dependent_coordinate_(2), dependent_speeds_{0, 2, 5}
{
}

void Bicycle::set_state(const state & x)
{
  state_ = x;
}

void Bicycle::set_coordinates(const coordinates & q)
{
  state_.block<8, 1>(0, 0) = q;
}

void Bicycle::set_speeds(const speeds & u)
{
  state_.block<12, 1>(8, 0) = u;
}
  
void Bicycle::set_parameters(const WheelAssemblyGyrostat & rear,
                             const WheelAssemblyGyrostat & front,
                             double ls, double g)
{
  rear_ = rear;
  front_ = front;
  ls_ = ls;
  g_ = g;
}

void Bicycle::set_dependent_coordinate(int i)
{
  if (!((i == 1) | (i == 2) | (i == 3))) {
    std::cerr << "Invalid dependent coordinate.  Ignoring request." << std::endl
              << "dependent coordinate must be 1, 2, or 3 (lean, pitch, steer)." << std::endl;
    return;
  }
  dependent_coordinate_ = i;
}

void Bicycle::set_dependent_speeds(int indices[3])
{
  if (indices[0] == indices[1] ||
      indices[0] == indices[2] ||
      indices[1] == indices[2]) {
    std::cerr << "Invalid dependent speeds.  Speed indices must be unique."
      << std::endl;
  } else {
    for (int i = 0; i < 3; ++i) {
      dependent_speeds_[i] = indices[i];
    }
  }
}

std::ostream & operator<<(std::ostream & os,
                          const Bicycle & b)
{
  os << "Rear assembly:" << std::endl
     << b.rear_
     << "Front assembly:" << std::endl
     << b.front_;
  os.precision(16);
  os << "ls  = " << std::setw(25) << b.ls_ << std::endl
     << "g   = " << std::setw(25) << b.g_  << std::endl
     << "State:" << std::endl 
     << " x[0] = " << std::setw(25) <<  b.state_[0] << "    (yaw)" << std::endl
     << " x[1] = " << std::setw(25) <<  b.state_[1] << "    (lean)" << std::endl
     << " x[2] = " << std::setw(25) <<  b.state_[2] << "    (pitch)" << std::endl
     << " x[3] = " << std::setw(25) <<  b.state_[3] << "    (steer)" << std::endl
     << " x[4] = " << std::setw(25) <<  b.state_[4] << "    (rear wheel angle)" << std::endl
     << " x[5] = " << std::setw(25) <<  b.state_[5] << "    (front wheel angle)" << std::endl
     << " x[6] = " << std::setw(25) <<  b.state_[6] << "    (x of rear wheel contact)" << std::endl
     << " x[7] = " << std::setw(25) <<  b.state_[7] << "    (y of rear wheel contact)" << std::endl
     << " x[8] = " << std::setw(25) <<  b.state_[8] << "    (yaw rate)" << std::endl
     << " x[9] = " << std::setw(25) <<  b.state_[9] << "    (lean rate)" << std::endl
     << "x[10] = " << std::setw(25) << b.state_[10] << "    (pitch rate)" << std::endl
     << "x[11] = " << std::setw(25) << b.state_[11] << "    (steer rate)" << std::endl
     << "x[12] = " << std::setw(25) << b.state_[12] << "    (rear wheel rate)" << std::endl
     << "x[13] = " << std::setw(25) << b.state_[13] << "    (front wheel rate)" << std::endl
     << "x[14] = " << std::setw(25) << b.state_[14] << "    (rear wheel contact longitudinal rate)" << std::endl
     << "x[15] = " << std::setw(25) << b.state_[15] << "    (rear wheel contact lateral rate)" << std::endl
     << "x[16] = " << std::setw(25) << b.state_[16] << "    (rear wheel contact vertical rate)" << std::endl
     << "x[17] = " << std::setw(25) << b.state_[17] << "    (front wheel contact longitudinal rate)" << std::endl
     << "x[18] = " << std::setw(25) << b.state_[18] << "    (front wheel contact lateral rate)" << std::endl
     << "x[19] = " << std::setw(25) << b.state_[19] << "    (front wheel contact vertical rate)" << std::endl;
  return os;
}
