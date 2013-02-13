#ifndef WHIPPLE_H
#define WHIPPLE_H
#include <cmath>

namespace bicycle {
/** Whipple bicycle model parameters.
 *
 * This class contains the parameters described in <a
 * href="http://dx.doi.org/10.1098/rspa.2007.1857">Meijaard et al.</a>. In
 * addition to the parameters described there, the rear and front tire cross
 * sectional radii are added as \f$t_R\f$ and \f$t_F\f$, respectively.
 *
 * */
class Whipple {
 public:
  /** Default constructor.
   *
   * Default values of all parameters are set to the values presented in Table
   * 1 of <a href="http://dx.doi.org/10.1098/rspa.2007.1857">Meijaard et
   * al.</a>
   * */
  Whipple() : w(1.02), c(0.08), lambda(M_PI/10.0), g(9.81), rR(0.3),
    mR(2.0), IRxx(0.0603), IRyy(0.12), xB(0.3), zB(-0.9), mB(85.0), IBxx(9.2),
    IByy(11.0), IBzz(2.8), IBxz(2.4), xH(0.9), zH(-0.7), mH(4.0), IHxx(0.05892),
    IHyy(0.06), IHzz(0.00708), IHxz(-0.00756), rF(0.35), mF(3.0), IFxx(0.1405),
    IFyy(0.28), tR(0.0), tF(0.0) {}

  double w,     /**< wheel base */
         c,     /**< trail */
         lambda,/**< steer axis tilt */
         g,     /**< gravity */
         rR,    /**< rear wheel radius */
         mR,    /**< rear wheel mass */
         IRxx,  /**< rear wheel central moment of inertia relative to any diameter in wheel plane*/
         IRyy,  /**< rear wheel central moment of inertia relative to spin axis */
         xB,    /**< rear frame position center of mass location in x direction from rear contact */
         zB,    /**< rear frame position center of mass location in z direction from rear contact */
         mB,    /**< rear frame mass */
         IBxx,  /**< rear frame central moment of inertia relative to x direction */
         IByy,  /**< rear frame central moment of inertia relative to y direction */
         IBzz,  /**< rear frame central moment of inertia relative to z direction */
         IBxz,  /**< rear frame central product of inertia relative to x and z directions */
         xH,    /**< front frame position center of mass location in x direction from rear contact */
         zH,    /**< front frame position center of mass location in z direction from rear contact */
         mH,    /**< front frame mass */
         IHxx,  /**< front frame central moment of inertia relative to x direction */
         IHyy,  /**< front frame central moment of inertia relative to y direction */
         IHzz,  /**< front frame central moment of inertia relative to z direction */
         IHxz,  /**< front frame central product of inertia relative to x and z directions */
         rF,    /**< front wheel radius */
         mF,    /**< front wheel mass */
         IFxx,  /**< front wheel central moment of inertia relative to any diameter in wheel plane*/
         IFyy,  /**< front wheel central moment of inertia relative to spin axis */
         tR,    /**< rear tire cross section radius */
         tF;    /**< rear tire cross section radius */
};

} // namespace bicycle

#endif
