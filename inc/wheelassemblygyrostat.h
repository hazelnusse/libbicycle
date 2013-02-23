#ifndef WHEELASSEMBLYGYROSTAT_H
#define WHEELASSEMBLYGYROSTAT_H
#include <iostream>

namespace bicycle {

/** Wheel assembly gyrostat parameters.
 *
 * This class contains the parameters used to describe a cylindrical gyrostat
 * that is symmetric about the x-z plane, and which has a rotor that is
 * represented by a torus.
 *
 * Additionally, a axle torque and torques applied to the carrier (in body
 * fixed coordinates) can be specified, as well as forces applied at the
 * contact point and mass center.  No coordinate information is stored, so this
 * class is simply a container for all the parameters and the applied forces
 * and torques.
 *
 * */
class WheelAssemblyGyrostat {
 public:
  /** Default constructor.
   *
   * All values are zero-initialized by default.
   */
  WheelAssemblyGyrostat();

  /** Stream insertion operator.
   *
   * \param[in,out] os Output stream.
   * \param[in] w WheelAssemblyGyrostat instance.
   *
   * Prints the parameters.
   */
  friend std::ostream & operator<<(std::ostream & os,
                                   const bicycle::WheelAssemblyGyrostat & w);


  double Ixx, /**< Gyrostat total central moment of inertia relative to x direction */
         Iyy, /**< Gyrostat total central moment of inertia relative to y direction */
         Izz, /**< Gyrostat total central moment of inertia relative to z direction */
         Ixz, /**< Gyrostat total central product of inertia relative to x and z directions */
         J,   /**< Gyrostat rotor moment of inertia relative to spin axis */
         m,   /**< Gyrostat total mass */
         R,   /**< Rotor major radius */
         r,   /**< Rotor minor radius */
         a,   /**< Gyrotstat position center of mass from rotor center relative to x direction */
         b,   /**< Gyrotstat position center of mass from rotor center relative to z direction */
         c,   /**< Steer axis location from rotor center relative to x direction */
         Tw,  /**< Torque applied to rotor from carrier about axle axis  */
         Tx,  /**< Torque applied to carrier about x direction */
         Ty,  /**< Torque applied to carrier about y direction */
         Tz,  /**< Torque applied to carrier about z direction */
         Gx,  /**< Contact force applied to rotor in x direction */
         Gy,  /**< Contact force applied to rotor in y direction */
         Gz,  /**< Contact force applied to rotor in z direction */
         Fx,  /**< Force applied to gyrostat mass center in x direction */
         Fy,  /**< Force applied to gyrostat mass center in y direction */
         Fz;  /**< Force applied to gyrostat mass center in z direction */
};

} // namespace bicycle

#endif

