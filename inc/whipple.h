#ifndef WHIPPLE_H
#define WHIPPLE_H
#include <cmath>

class WhippleParameters {
 public:
  WhippleParameters() : w(1.02), c(0.08), lambda(M_PI/10.0), g(9.81), rR(0.3),
    mR(2.0), IRxx(0.0603), IRyy(0.12), xB(0.3), zB(-0.9), mB(85.0), IBxx(9.2),
    IByy(11.0), IBzz(2.8), IBxz(2.4), xH(0.9), zH(-0.7), mH(4.0), IHxx(0.05892),
    IHyy(0.06), IHzz(0.00708), IHxz(-0.00756), rF(0.35), mF(3.0), IFxx(0.1405),
    IFyy(0.28) {}

  double w, c, lambda, g, rR, mR, IRxx, IRyy, xB, zB, mB, IBxx, IByy, IBzz,
         IBxz, xH, zH, mH, IHxx, IHyy, IHzz, IHxz, rF, mF, IFxx, IFyy, tR, tF;
};

#endif
