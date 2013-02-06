#ifndef WHEELASSEMBLYGYROSTAT_H
#define WHEELASSEMBLYGYROSTAT_H

class WheelAssemblyGyrostat {
 public:
   WheelAssemblyGyrostat() : Ixx(0.0), Iyy(0.0), Izz(0.0), Ixz(0.0), J(0.0),
                               m(0.0), R(0.0), r(0.0), a(0.0), b(0.0), c(0.0) {}
   double Ixx, Iyy, Izz, Ixz,
          J, m, R, r, a, b, c;
};
#endif
