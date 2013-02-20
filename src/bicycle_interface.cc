#include <algorithm>
#include <set>
#include <iostream>
#include <iomanip>
#include "bicycle.h"

namespace bicycle {

// We must have definitions of static const members even if we define them in
// the header.  Weird linker errors can (and have occured if they are declared
// but not defined.
const int Bicycle::n;
const int Bicycle::n_min;
const int Bicycle::l;
const int Bicycle::o;
const int Bicycle::m;
const int Bicycle::s;

Bicycle::Bicycle()
  : state_(Eigen::Matrix<double, n + o, 1>::Zero()), ls_(0.0), g_(9.81),
  steer_torque_(0.0), dependent_coordinate_(2), dependent_speeds_{0, 2, 5},
  azimuth(0.0), elevation(0.0), twist(0.0), cam_x(0.0), cam_y(0.0), cam_z(0.0)
{
  update_permutations();
}

void Bicycle::set_state(const Vector & x)
{
  state_ = x;
}

void Bicycle::set_coordinates(const Vector & q)
{
  state_.block<n, 1>(0, 0) = q;
}

Vector Bicycle::coordinates() const
{
  return state_.block<n, 1>(0, 0);
}

double Bicycle::coordinate(int i) const
{
  return state_[i];
}

void Bicycle::set_speeds(const Vector & u)
{
  state_.block<o, 1>(n, 0) = u;
}

Vector Bicycle::speeds() const
{
  return state_.block<o, 1>(n, 0);
}

double Bicycle::speed(int i) const
{
  return state_[n + i];
}

void Bicycle::set_inputs(const Vector & r)
{
  rear_.Tw = r[0];
  rear_.Tx = r[1];
  rear_.Ty = r[2];
  rear_.Tz = r[3];
  rear_.Gx = r[4];
  rear_.Gy = r[5];
  rear_.Gz = r[6];
  rear_.Fx = r[7];
  rear_.Fy = r[8];
  rear_.Fz = r[9];
  front_.Tw = r[10];
  front_.Tx = r[11];
  front_.Ty = r[12];
  front_.Tz = r[13];
  front_.Gx = r[14];
  front_.Gy = r[15];
  front_.Gz = r[16];
  front_.Fx = r[17];
  front_.Fy = r[18];
  front_.Fz = r[19];
  steer_torque_ = r[20];
  g_ = r[21];
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
              << "dependent coordinate must be 1, 2, or 3 (lean, pitch, steer)."
              << std::endl;
    return;
  }
  dependent_coordinate_ = i;
  update_coordinate_permutation();
}

void Bicycle::set_dependent_speeds(const std::set<int> & speed_indices)
{
  if (speed_indices.size() != m) {
    std::cerr << "Three dependent speeds must be specified." << std::endl;
  } else if (std::all_of(speed_indices.begin(), speed_indices.end(),
               [](int i){ return (i < 0) || (i > 5); })) {
    std::cerr << "Dependent speeds must be in [0, 5]." << std::endl;
  } else if(speed_indices.count(3)) {
    std::cerr << "Steer rate must be an independent speed." << std::endl;
  } else if(!speed_indices.count(4) && !speed_indices.count(5)) {
    std::cerr << "At least one wheel rate must be dependent." << std::endl;
  } else {
    dependent_speeds_ = speed_indices;
    update_speed_permutation();
  }
}

void Bicycle::update_coordinate_permutation()
{
  ::Eigen::Matrix<int, n, 1> indices;
  std::iota(indices.data(), indices.data() + n, 0); // fill with 0...7
  // Move independent coordinate to front, dependent coordinate to back,
  // preserving relative order
  std::stable_partition(indices.data(), indices.data() + n,
                        [&](int si) { return si != dependent_coordinate_;});

  P_q_ = PermutationMatrix<n>(indices).toDenseMatrix().cast<double>();
}

void Bicycle::update_speed_permutation()
{
  ::Eigen::Matrix<int, o, 1> indices;
  std::iota(indices.data(), indices.data() + o, 0); // fill with 0...11
  // Move independent speeds to front, dependent speeds to back, preserving
  // relative order
  std::stable_partition(indices.data(), indices.data() + o,
                      [&](int si) { return !dependent_speeds_.count(si); });

  P_u_ = PermutationMatrix<o>(indices).toDenseMatrix().cast<double>();
}

void Bicycle::update_permutations()
{
  update_coordinate_permutation();
  update_speed_permutation();
}

Vector Bicycle::all_inputs_except_constraint_forces() const
{
  Vector r(15);
  r[0] = rear_.Tw;
  r[1] = rear_.Tx;
  r[2] = rear_.Ty;
  r[3] = rear_.Tz;
  r[4] = rear_.Fx;
  r[5] = rear_.Fy;
  r[6] = rear_.Fz;
  r[7] = front_.Tw;
  r[8] = front_.Tx;
  r[9] = front_.Ty;
  r[10] = front_.Tz;
  r[11] = front_.Fx;
  r[12] = front_.Fy;
  r[13] = front_.Fz;
  r[14] = g_;
  return r;
}

std::ostream & operator<<(std::ostream & os, const Bicycle & b)
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
