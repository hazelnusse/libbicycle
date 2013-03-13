#ifndef BICYCLE_H
#define BICYCLE_H
#include <iostream>
#include <set>
#include <utility>
#include <Eigen/Dense>
#include "wheelassemblygyrostat.h"
#include "whipple.h"

#ifdef WRAP_PYTHON
#include <Python.h>
#endif

namespace bicycle {

using ::Eigen::ComputeThinV;
using ::Eigen::FullPivHouseholderQR;
using ::Eigen::JacobiSVD;
using ::Eigen::Map;
using ::Eigen::PermutationMatrix;
using ::Eigen::Stride;
using ::Eigen::Unaligned;

/** Dynamically sized column vector of type double
 */
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

/** Dynamically sized matrix of type double
 */
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

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

  /** Default constructor.
   *
   * Initializes a bicycle with zero state and zero for all parameters except
   * gravity which is set by default to 9.81 m/s^2.  The dependent coordinate
   * is set to the pitch angle \f$\theta\f$, and the dependent speeds are set
   * to be the yaw rate \f$\dot{\psi}\f$, pitch rate \f$\dot{\theta}\f$, and
   * front wheel rate \f$\dot{\theta}_F\f$.
   */
  Bicycle();

  /** Set i-th state to xi.
   *
   * \param[in] i Index of state to set.
   * \param[in] xi Value assigned to state i.
   */
  void set_state(int i, double xi);

  /** Set coordinates.
   *
   * \param[in] q Vector of length 8 with the following ordering:
   *
   * - Rear assembly yaw angle \f$q_0\f$
   * - Rear assembly roll angle \f$q_1\f$
   * - Rear assembly pitch angle \f$q_2\f$
   * - Steer angle \f$q_3\f$
   * - Rear wheel angle \f$q_4\f$
   * - Front wheel angle \f$q_5\f$
   * - Rear wheel contact point \f$q_6\f$
   * - Rear wheel contact point \f$q_7\f$
   *
   * Note that the condition that the front wheel touch the ground plane
   * imposes a nonlinear relationship between roll \f$q_1\f$, pitch \f$q_2\f$,
   * and steer \f$q_3\f$.  The default choice for the dependent coordinate is
   * \f$q_2\f$ though for highly leaned configurations this may not be
   * appropriate.
   */
  void set_coordinates(const Vector & q);

  /** Set coordinates following Basu-Mandal et al. 2007 definitions.
   *
   * \param[in] q Vector of length 9 with the following ordering:
   *
   * - Rear wheel center \f$x\f$
   * - Rear wheel center \f$y\f$
   * - Rear wheel center \f$z\f$
   * - Yaw \f$\theta\f$
   * - Lean \f$\psi\f$
   * - Pitch \f$\phi\f$
   * - Steer \f$\psi_f\f$
   * - Rear wheel angle \f$\beta_r\f$
   * - Front wheel angle \f$\beta_f\f$
   */
  void set_coordinates_basu_mandal(const Vector & q_bm);

  /** Set coordinate time derivatives following Basu-Mandal et al. 2007 definitions.
   *
   * \param[in] q_dot Vector of length 9 with the following ordering:
   *
   * - Rear wheel center \f$\dot{x}\f$
   * - Rear wheel center \f$\dot{y}\f$
   * - Rear wheel center \f$\dot{z}\f$
   * - Yaw \f$\dot{\theta}\f$
   * - Lean \f$\dot{\psi}\f$
   * - Pitch \f$\dot{\phi}\f$
   * - Steer \f$\dot{\psi}_f\f$
   * - Rear wheel angle \f$\dot{\beta}_r\f$
   * - Front wheel angle \f$\dot{\beta}_f\f$
   *
   * Note that the Bicycle model uses a smaller number of speeds (6) than
   * Basu-Mandal et al. does (9), so the velocity measure numbers of the rear
   * wheel center in the inertial frame (\f$\dot{x}, \dot{y}, \dot{z}\f$) are
   * not used because this information is redundant.
   */
  void set_speeds_basu_mandal(const Vector & q_dot_bm);

  /** Get coordinates
   *
   * \returns A Vector of length 8 with the following ordering:
   *
   * - Rear assembly yaw angle \f$q_0\f$
   * - Rear assembly roll angle \f$q_1\f$
   * - Rear assembly pitch angle \f$q_2\f$
   * - Steer angle \f$q_3\f$
   * - Rear wheel angle \f$q_4\f$
   * - Front wheel angle \f$q_5\f$
   * - Rear wheel contact point \f$q_6\f$
   * - Rear wheel contact point \f$q_7\f$
   */
  Vector coordinates() const;

  /** Return i-th coordinate
   *
   */
  double coordinate(int i) const;

  /** Set i-th coordinate to qi.
   *
   * \param[in] i Index of coordinate to set.
   * \param[in] qi Value assigned to coordinate i.
   */
  void set_coordinate(int i, double qi);

  /** Set i-th speed to ui.
   *
   * \param[in] i Index of speed to set.
   * \param[in] ui Value assigned to speed i.
   *
   * The 12 speeds are organzied as follows
   *
   * - Rear assembly yaw angle rate \f$\dot{q}_0\f$
   * - Rear assembly roll angle rate \f$\dot{q}_1\f$
   * - Rear assembly pitch angle rate \f$\dot{q}_2\f$
   * - Steer angle rate \f$\dot{q}_3\f$
   * - Rear wheel angle rate \f$\dot{q}_4\f$
   * - Front wheel angle rate \f$\dot{q}_5\f$
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
   * the dependent speeds is yaw rate \f$\dot{q}_0\f$, pitch rate
   * \f$\dot{q}_2\f$, and front wheel rate \f$\dot{q}_5\f$.
   */
  void set_speed(int i, double ui);
  
  /* Set generalized speeds.
   *
   */
  void set_speeds(const Vector & u);

  /** Get the i-th generalized speed
   *
   */
  double speed(int i) const;

  /** Get generalized speeds.
   **/
  Vector speeds() const;

  /** Get system state.
   **/
  Vector state() const;

  /** Set exogenous inputs
   *
   * \param[in] r Vector of length 22 with the following ordering
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
  void set_inputs(const Vector & r);

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

  /** Get rear assembly parameters
   *
   */
  WheelAssemblyGyrostat rear_parameters() const;

  /** Get front assembly parameters
   *
   */
  WheelAssemblyGyrostat front_parameters() const;

  /** Get steer axis offset
   *
   */
  double steer_axis_offset() const;

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

  /** Compute ground contact forces and steer torque under steady turning
   * conditions and assuming no slip.
   *
   * \pre State and parameters of bicycle are set, and the speeds associated
   * with tire contact velocities are zero.
   *
   * \returns The constraint forces as a length 7 Vector with ordering:
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
   */
  Vector steady_no_slip_constraint_forces() const;

  /** Contact forces assuming no slip.
   *
   * \pre State and parameters are set, input forces are applied.
   *
   * \returns The constraint forces acting at each contact point.
   *
   * Note that this doesn't assume that the motion is steady, i.e. this is for
   * arbitrary motions.
   *
   */
  Vector no_slip_contact_forces() const;

  /** Solve configuration constraint.
   *
   * \param[in] ftol Tolerance used to terminate root finder.
   * \param[in] iter Maximum number of iterations.
   * \returns Number of number of Newton iterations completed.
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
  std::pair<int, double> solve_configuration_constraint_and_set_state(double ftol=1e-14,
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
   * \returns Residual of solution:  B_i * u_d + B_d * u_i
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
  Vector solve_velocity_constraints_and_set_state();

  /** Form linearized mass matrix
   *
   * \returns a 20 x 20 matrix of dq/dt and du/dt coefficients
   */
  Matrix mass_matrix_full() const;
  
  /** Form linearized input coefficient matrix
   *
   * \returns a 20 x 22 matrix
   */
  Matrix input_matrix() const { return B_u(); }

  /** Form linearized state matrix
   *
   * \returns a 20 x 16 matrix of dq_i/dt and du_i/dt coefficients
   */
  Matrix independent_state_matrix() const;

  /** Coordinate derivatives
   *
   */
  Vector coordinate_derivatives() const;
  
  /** Speed derivatives
   *
   */
  Vector speed_derivatives() const;

  /** State derivatives
   *
   */
  Vector state_derivatives() const;

  /** Position of points
   *
   * \returns a 7 x 3 matrix with the following seven points, all relative to
   * the rear wheel contact point, expressed in the rear wheel yaw frame
   * coordinates:
   *  - Rear wheel center
   *  - Rear assembly mass center
   *  - Rear assembly steer axis point
   *  - Front wheel center
   *  - Front assembly mass center
   *  - Front assembly steer axis point
   *  - Front wheel contact point
   *
   */
  Matrix points_of_interest() const;

  /** Reference pitch angle.
   *
   * \returns pitch angle when bicycle has zero lean and steer.
   */
  double reference_pitch() const;

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
  void rear_wheel_center_point(double ar[3]) const;
  void rear_mass_center_point(double ar[3]) const;
  void rear_steer_axis_point(double ar[3]) const;
  void front_wheel_center_point(double ar[3]) const;
  void front_mass_center_point(double ar[3]) const;
  void front_steer_axis_point(double ar[3]) const;
  void front_ground_contact_point(double ar[3]) const;
  void q6q7_from_bm(double ar[2], double x_bm, double y_bm) const;
  void f_c(double ar[1]) const;
  void f_c_dq(double ar[n]) const;
  void f_v_du(double ar[m * o]) const;
  void f_v_dudq(double ar[m * o * n_min]) const;
  void f_v_dudt(double ar[36]) const;
  void f_v_dudtdq(double ar[108]) const;
  void kinematic_odes_rhs(double ar[n]) const;
  void f_1(double ar[n]) const;
  void f_1_dq(double ar[n * n]) const;
  void f_1_du(double ar[n * o]) const;
  void gif_dud(double ar[o * o]) const;
  void gif_dud_dq(double ar[o * o * n_min]) const;
  void gif_ud_zero(double ar[o]) const;
  void gif_ud_zero_steady(double ar[o]) const;
  void gif_ud_zero_dq(double ar[o * n_min]) const;
  void gif_ud_zero_du(double ar[o * o]) const;
  void gaf(double ar[o]) const;
  void gaf_dq(double ar[o * n_min]) const;
  void gaf_dr(double ar[o * s]) const;
  void path_radii(double ar[2]) const;
  void ke_pe(double ar[6]) const;
  void xyz_dot_bm(double ar[3]) const;
  void gif_ud_zero_steady_dudu(double ar[84]) const;
  void gif_ud_zero_steady_cross_terms(double ar[7]) const;

  // Private member functions related to linearization of dynamic equations
  Matrix Bd_inverse_Bi() const;
  Matrix f_v_dq() const;
  Matrix M_qq() const;
  Matrix M_uqc() const;
  Matrix M_uuc() const;
  Matrix M_uqd() const;
  Matrix M_uud() const;
  Matrix A_qq() const;
  Matrix A_qu() const;
  Matrix A_uqc() const;
  Matrix A_uuc() const;
  Matrix A_uqd() const;
  Matrix A_uud() const;
  Matrix B_u() const;

  Matrix C_0() const;
  Matrix C_1() const;
  Matrix C_2() const;

  Matrix P_qd() const;
  Matrix P_qi() const;

  Matrix P_ud() const;
  Matrix P_ui() const;

  // Private members used for convenience
  bool is_dependent_index(int i) const;
  void update_coordinate_permutation();
  void update_speed_permutation();
  void update_permutations();
  Vector all_inputs_except_constraint_forces() const;
  int best_dependent_coordinate() const;
  std::set<int> best_dependent_speeds() const;

  // Private data
  Eigen::Matrix<double, n + o, 1> state_;// Bicycle state, q's then u's
  WheelAssemblyGyrostat rear_, front_; // Wheel gyrostats
  double ls_, g_, steer_torque_;       // Steer axis offset, gravity, torque

  int dependent_coordinate_;
  std::set<int> dependent_speeds_;
  Matrix P_q_, P_u_;           // Permutation matrices
  // Camera related variables
  double azimuth, elevation, twist, cam_x, cam_y, cam_z;

#if WRAP_PYTHON
 public:
  void solve_configuration_constraint_and_set_state_python();
  void solve_velocity_constraints_and_set_state_python();
  void mass_matrix_full_python(PyObject * matrix_out) const;
  void independent_state_matrix_python(PyObject * matrix_out) const;
  void input_matrix_python(PyObject * matrix_out) const;
#endif

};

} // namespace bicycle

#include "bicycle-inl.h"
#endif
