#ifndef BICYCLE_H
#define BICYCLE_H
#include <iostream>
#include <set>
#include <Eigen/Dense>
#include "wheelassemblygyrostat.h"
#include "whipple.h"

namespace bicycle {

typedef double Real;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMajorMatrix;
//typedef RowMajorMatrix Matrix;
using namespace Eigen;

/** A class representing a dynamic model of a bicycle.
 *
 * This class encapsulates the state and the physical parameters in addition to
 * useful functionality for simulating bicycle dynamics, solving constraints,
 * constraint forces, determining equilibria, and the stability of those
 * equilibria.
 *
 * Most member functions do not have side effects; i.e., they do not modify the
 * state or the parameters.  Users of the class should set the state and
 * parameters to what they want, then request that calculations be performed.
 */
class Bicycle {
 public:
  /** Number of generalized coordinates.
   *
   * Includes dependent and independent coordinates.  This includes cyclic
   * coordinates.
   */
  static const int n = 8;

  /** Number of configuration constraints.
   */
  static const int l = 1;

  /** Number of non ignorable coordinates
   */
  static const int n_min = 3;

  /** Number of generalized speeds.
   *
   * Includes dependent and independent speeds.
   */
  static const int o = 12;

  /** Number of velocity constraints.
   *
   */
  static const int m = 3;

  /** Number of exogenous inputs.
   *
   */
  static const int s = 22;

  /** 8 x 1 matrix of doubles with ordering:
   *
   * - Rear assembly yaw angle \f$\psi\f$
   * - Rear assembly roll angle \f$\phi\f$
   * - Rear assembly pitch angle \f$\theta\f$
   * - Steer angle \f$\delta\f$
   * - Rear wheel angle \f$\theta_R\f$
   * - Front wheel angle \f$\theta_F\f$
   * - Rear wheel contact point \f$x\f$
   * - Rear wheel contact point \f$y\f$
   *
   * Note that the condition that the front wheel touch the ground plane
   * imposes a nonlinear relationship between \f$\phi\f$, \f$\theta\f$, and
   * \f$\delta\f$.  The default choice for the dependent coordinate is
   * \f$\theta\f$ though for highly leaned configurations this may not be
   * appropriate.
   */
  typedef Eigen::Matrix<double, n, 1> coordinates;

  /** 12 x 1 matrix of doubles with ordering:
   *
   * - Rear assembly yaw angle rate \f$\dot{\psi}\f$
   * - Rear assembly roll angle rate \f$\dot{\phi}\f$
   * - Rear assembly pitch angle rate \f$\dot{\theta}\f$
   * - Steer angle \f$\dot{\delta}\f$
   * - Rear wheel angle \f$\dot{\theta}_R\f$
   * - Front wheel angle \f$\dot{\theta}_F\f$
   * - Rear wheel contact point longitudinal velocity \f$v_{Rx}\f$
   * - Rear wheel contact point lateral velocity \f$v_{Ry}\f$
   * - Rear wheel contact point vertical velocity \f$v_{Rz}\f$
   * - Front wheel contact point longitudinal velocity \f$v_{Fx}\f$
   * - Front wheel contact point lateral velocity \f$v_{Fy}\f$
   * - Front wheel contact point vertical velocity \f$v_{Fz}\f$
   *
   * Note that the condition that the common steer axis point of the front and
   * rear assemblies have the same velocity imposes a linear relationship
   * between all speeds.  Three scalar equations relate these twelve speeds;
   * and hence only nine can be chosen independently.  The default choice for
   * the dependent speeds is \f$\dot{\psi}\f$, \f$\dot{\theta}\f$,
   * \f$\dot{\theta}_F\f$.
   */
  typedef Eigen::Matrix<double, o, 1> speeds;

  /** 20 x 1 matrix of doubles with ordering:
   *
   * - Generalized coordinates (8 x 1 matrix of doubles)
   * - Generalized speeds (12 x 1 matrix of doubles)
   */
  typedef Eigen::Matrix<double, n + o, 1> state;

  /** 22 x 1 matrix of doubles with ordering:
   *
   * - Rear wheel axle torque
   * - Rear assembly torque about x direction
   * - Rear assembly torque about y direction
   * - Rear assembly torque about z direction
   * - Rear tire ground contact forces in x direction
   * - Rear tire ground contact forces in y direction
   * - Rear tire ground contact forces in z direction
   * - Rear assembly force at mass center in x direction
   * - Rear assembly force at mass center in y direction
   * - Rear assembly force at mass center in z direction
   * - Front wheel axle torque
   * - Front assembly torque about x direction
   * - Front assembly torque about y direction
   * - Front assembly torque about z direction
   * - Front tire ground contact forces in x direction
   * - Front tire ground contact forces in y direction
   * - Front tire ground contact forces in z direction
   * - Front assembly force at mass center in x direction
   * - Front assembly force at mass center in y direction
   * - Front assembly force at mass center in z direction
   * - Steer torque
   * - Gravitational acceleration
   *
   * The last entry may seem a bit counterintuitive, but if gravity is viewed
   * as an input that can vary, then it makes sense to consider it as an input.
   *
   */
  typedef Eigen::Matrix<double, s, 1> inputs;

  /** Default constructor.
   *
   * Initializes a bicycle with zero state and zero for all parameters except
   * gravity which is set by default to 9.81 m/s^2.  The dependent coordinate
   * is set to the pitch angle \f$\theta\f$, and the dependent speeds are set
   * to be the yaw rate \f$\dot{\psi}\f$, pitch rate \f$\dot{\theta}\f$, and
   * front wheel rate \f$\dot{\theta}_F\f$.
   */
  Bicycle();

  /** Set bicycle state.
   *
   * \param[in] x The state to set the bicycle to.
   *
   * No checking is performed to determine if the assigned state satisfies the
   * configuration or velocity constraints.
   */
  void set_state(const state & x);

  /** Set i-th state to xi.
   *
   * \param[in] i Index of state to set.
   * \param[in] xi Value assigned to state i.
   */
  void set_state(int i, double xi);

  /** Set coordinates.
   *
   * \param[in] q Coordinates to set.
   */
  void set_coordinates(const coordinates & q);

  /** Set i-th coordinate to qi.
   *
   * \param[in] i Index of coordinate to set.
   * \param[in] qi Value assigned to coordinate i.
   */
  void set_coordinate(int i, double qi);

  /** Set speeds.
   *
   * \param[in] u Speeds to set.
   */
  void set_speeds(const speeds & u);

  /** Set i-th speed to ui.
   *
   * \param[in] i Index of speed to set.
   * \param[in] ui Value assigned to speed i.
   */
  void set_speed(int i, double ui);

  /** Set the physical parameters.
   *
   * \param[in] rear Rear wheel assembly gyrostat.
   * \param[in] front Front wheel assembly gyrostat.
   * \param[in] ls Steer axis offset.
   * \param[in] g Acceleration due to gravity.
   */
  void set_parameters(const WheelAssemblyGyrostat & rear,
                      const WheelAssemblyGyrostat & front,
                      double ls, double g);

  /** Set physical parameters from Whipple model.
   *
   * \param[in] w Whipple parameters.
   */
  void set_parameters(const Whipple & w);

  /** Set the dependent coordinate.
   *
   * \param[in] i Index of dependent coordinate.
   *
   * Must be 1, 2, or 3, which
   * correspond to lean \f$\phi\f$, pitch \f$\theta\f$, and steer \f$\delta\f$,
   * respectively.  This setting affects which coordinate is solved for with 
   * Bicycle::solve_configuration_constraint_and_set_state.
   *
   * \see Bicycle::solve_configuration_constraint_and_set_state
   */
  void set_dependent_coordinate(int i);

  /** Set the dependent speeds.
   *
   * \param[in] s The set of dependent speeds indices.
   */
  void set_dependent_speeds(const std::set<int> & s);

  /** Stream insertion operator.
   *
   * \param[in,out] os Output stream.
   * \param[in] b Bicycle instance.
   *
   * Prints the parameters and state.
   */
  friend std::ostream & operator<<(std::ostream & os,
                                   const Bicycle & b);

  /** Compute steady ground contact forces assuming \f$\dot{u}=0\f$.
   *
   * \pre State and parameters of bicycle are set.
   *
   * \returns The constraint forces as a 7 x 1 matrix with ordering:
   * - \f$G_{Rx}\f$
   * - \f$G_{Ry}\f$
   * - \f$G_{Rz}\f$
   * - \f$G_{Fx}\f$
   * - \f$G_{Fy}\f$
   * - \f$G_{Fz}\f$
   * - \f$\tau_s\f$
   *
   * The first six components of this matrix are with respect to the rear and
   * front wheel yaw frames, respectively. x is parallel to the intersection
   * line of the wheel plane and the ground plane, y is perpendicular to x and
   * in the ground plane, z is perpendicular to the ground plane and points
   * downwards (to the half space opposite that which the bicycle is normally
   * in).
   * */
  RowMajorMatrix steady_constraint_forces() const;

  /** Solve configuration constraint.
   *
   * \param[in] ftol Tolerance used to terminate root finder.
   * \param[in] iter Maximum number of iterations.
   *
   * \pre The dependent coordinate has been properly selected and the value of
   * that coordinate has been set as an initial guess for the root finder.
   *
   * \post If the root finder converges, the dependent coordinate is set to
   * that root.
   *
   * The configuration constraint that the front wheel touch the ground is of
   * the form:
   * \f[
   *   f(\phi, \theta, \delta) = 0
   * \f]
   *
   * A Newton-Raphson root finding method is applied to the configuration
   * constraint, taking the currently set dependent coordinate as the initial
   * condition for the root finder.  If the root finder doesn't converge the
   * state is not changed.  Failure to converge occurs primarily when the wrong
   * dependent coordinate is selected; for some configurations, changing some
   * coordinates has no effect on the front wheel contact height.  Another
   * issue that can occur is that due to the existence of multiple roots, a
   * poor initial guess can lead to convergence to a root other than the
   * desired one.
   */
  void solve_configuration_constraint_and_set_state(double ftol=1e-14,
                                                    int iter=20);

  /** Solve velocity constraints.
   *
   * \pre The coordinates \f$\phi\f$, \f$\theta\f$, and \f$\delta\f$ satisfy
   * the configuration constraint, and a valid choice of dependent speeds has
   * been selected.
   *
   * \post The dependent speeds are solved for and the corresponding state
   * variables are set.
   *
   * The velocity constraints are of the form:
   * \f[
   *   B(q) u = 0
   * \f]
   * Given a set of 3 dependent speeds, we can rearrange this as:
   * \f[
   *   B_{d}(q) u_d + B_{i}(q) u_i = 0
   * \f]
   * which allows the dependent speeds to be computed by solving the linear
   * system:
   * \f[
   *   B_{d}(q) u_d = -B_{i}(q) u_i
   * \f]
   *
   * for \f$u_d\f$. An invalid choice of dependent speeds cause \f$B_d(q)\f$ to
   * be singular or poorly conditioned.
   */
  void solve_velocity_constraints_and_set_state();

  /** Form linearized full mass matrix
   *
   * \returns a 20 x 20 coefficient matrix of dq/dt and du/dt
   */
  RowMajorMatrix mass_matrix_full() const;

  /** Form linearized state matrix
   *
   * \returns a 20 x 16 coefficient matrix
   */
  RowMajorMatrix independent_state_matrix() const;

  // This is to ensure state has 128-bit alignment and hence vectorizable
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  // Generated code, don't expose as this is an implementation detail
  void gc_r_ogl(double ar[16]) const;
  void wc_r_ogl(double ar[16]) const;
  void mc_r_ogl(double ar[16]) const;
  void gc_f_ogl(double ar[16]) const;
  void wc_f_ogl(double ar[16]) const;
  void mc_f_ogl(double ar[16]) const;
  void N_ogl(double ar[16]) const;
  void f_c(double ar[1]) const;
  void f_c_dq(double ar[8]) const;
  void f_v_coefficient(double ar[m * o]) const;
  void f_v_coefficient_dq(double ar[m * o * n_min]) const;
  void f_v_coefficient_dqdq(double ar[m * o * n_min * n_min]) const;
  void kinematic_odes_rhs(double ar[n]) const;
  void f_1(double ar[n]) const;
  void f_1_dq(double ar[n * n]) const;
  void f_1_du(double ar[n * o]) const;
  void gif_dud(double ar[o * o]) const;
  void gif_dud_dq(double ar[o * o * n_min]) const;
  void gif_ud_zero(double ar[o]) const;
  void gif_ud_zero_steady(double ar[o]) const;
  void gif_ud_zero_dq(double ar[o * n]) const;
  void gif_ud_zero_du(double ar[o * o]) const;
  void gaf(double ar[o]) const;
  void gaf_dq(double ar[o * n]) const;
  void gaf_dr(double ar[o * s]) const;

  // Private member functions related to linearization of dynamic equations
  RowMajorMatrix Bd_inverse_Bi() const;
  RowMajorMatrix f_v_dq() const;
  RowMajorMatrix M_qq() const;
  RowMajorMatrix M_uqc() const;
  RowMajorMatrix M_uuc() const;
  RowMajorMatrix M_uqd() const;
  RowMajorMatrix M_uud() const;
  RowMajorMatrix A_qq() const;
  RowMajorMatrix A_qu() const;
  RowMajorMatrix A_uqc() const;
  RowMajorMatrix A_uuc() const;
  RowMajorMatrix A_uqd() const;
  RowMajorMatrix A_uud() const;

  RowMajorMatrix C_0() const;
  RowMajorMatrix C_1() const;
  RowMajorMatrix C_2() const;

  RowMajorMatrix P_qd() const;
  RowMajorMatrix P_qi() const;

  RowMajorMatrix P_ud() const;
  RowMajorMatrix P_ui() const;

  // Private members used for convenience
  bool is_dependent_index(int i) const;
  void update_coordinate_permutation();
  void update_speed_permutation();
  void update_permutations();
  RowMajorMatrix all_inputs_except_constraint_forces() const;
  int best_dependent_coordinate() const;
  std::set<int> best_dependent_speeds() const;

  // Private data
  state state_;                        // Bicycle state, q's then u's
  WheelAssemblyGyrostat rear_, front_; // Wheel gyrostats
  double ls_, g_, steer_torque_;       // Steer axis offset, gravity, torque

  int dependent_coordinate_;
  std::set<int> dependent_speeds_;
  RowMajorMatrix P_q_, P_u_;           // Permutation matrices
  // Camera related variables
  double azimuth, elevation, twist, cam_x, cam_y, cam_z;
};

} // namespace bicycle

#include "bicycle-inl.h"
#endif
