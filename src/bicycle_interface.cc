#include <algorithm>
#include <set>
#include <iostream>
#include <iomanip>
#include "bicycle.h"

namespace bicycle {

// We must have definitions of static const members even if we define them in
// the header.  Weird linker errors can (and have occured if they are declared
// but not defined.
const int Bicycle::kNumberOfCoordinates;
const int Bicycle::kNumberOfConfigurationConstraints;
const int Bicycle::kNumberOfSpeeds;
const int Bicycle::kNumberOfVelocityConstraints;
const int Bicycle::kNumberOfInputs;

Bicycle::Bicycle()
  : state_(state::Zero()), ls_(0.0), g_(9.81), steer_torque_(0.0),
    dependent_coordinate_(2), dependent_speeds_{0, 2, 5},
    azimuth(0.0), elevation(0.0), twist(0.0),
    cam_x(0.0), cam_y(0.0), cam_z(0.0)
{
  update_permutations();
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
  update_coordinate_permutation();
}

void Bicycle::set_dependent_speeds(const std::set<int> & speed_indices)
{
  if (speed_indices.size() != kNumberOfVelocityConstraints) {
    std::cerr << "Three dependent speeds must be specified." << std::endl;
  } else if (std::all_of(speed_indices.begin(),
               speed_indices.end(),
               [](int i){
               return (i < 0) || (i > 5);})) {
    std::cerr << "Dependent speeds must be in [0, 5]." << std::endl;
  } else if(speed_indices.count(3)) {
    std::cerr << "Steer rate must be an independent speed." << std::endl;
  } else if(speed_indices.count(4) && speed_indices.count(5)) {
    std::cerr << "Only one wheel rate may be dependent." << std::endl;
  } else if(!speed_indices.count(4) && !speed_indices.count(5)) {
    std::cerr << "At least one wheel rate must be dependent." << std::endl;
  } else {
    dependent_speeds_ = speed_indices;
    update_speed_permutation();
  }
}

void Bicycle::update_coordinate_permutation()
{
  Eigen::Matrix<int, kNumberOfCoordinates, 1> s;
  std::iota(s.data(), s.data() + kNumberOfCoordinates, 0); // fill with 0...11
  // Move independent coordinate to front, dependent coordinate to back,
  // preserving relative order
  std::stable_partition(s.data(),
      s.data() + kNumberOfCoordinates,
      [&](int elem) { return elem != dependent_coordinate_;});

  P_q_ = PermutationMatrix<kNumberOfCoordinates>(s).toDenseMatrix().cast<double>();
}

void Bicycle::update_speed_permutation()
{
  Eigen::Matrix<int, kNumberOfSpeeds, 1> s;
  std::iota(s.data(), s.data() + kNumberOfSpeeds, 0); // fill with 0...11
  // Move independent speeds to front, dependent speeds to back, preserving
  // relative order
  std::stable_partition(s.data(),
      s.data() + kNumberOfSpeeds,
      [&](int elem) { return !dependent_speeds_.count(elem); });

  P_u_ = PermutationMatrix<kNumberOfSpeeds>(s).toDenseMatrix().cast<double>();
}

void Bicycle::update_permutations()
{
  update_coordinate_permutation();
  update_speed_permutation();
}

RowMajorMatrix Bicycle::all_inputs_except_constraint_forces() const
{
  RowMajorMatrix r(15, 1);
  r(0, 0) = rear_.Tw;
  r(1, 0) = rear_.Tx;
  r(2, 0) = rear_.Ty;
  r(3, 0) = rear_.Tz;
  r(4, 0) = rear_.Fx;
  r(5, 0) = rear_.Fy;
  r(6, 0) = rear_.Fz;
  r(7, 0) = front_.Tw;
  r(8, 0) = front_.Tx;
  r(9, 0) = front_.Ty;
  r(10, 0) = front_.Tz;
  r(11, 0) = front_.Fx;
  r(12, 0) = front_.Fy;
  r(13, 0) = front_.Fz;
  r(14, 0) = g_;
  return r;
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

} // namespace bicycle
