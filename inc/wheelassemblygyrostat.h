#ifndef WHEELASSEMBLYGYROSTAT_H
#define WHEELASSEMBLYGYROSTAT_H

class WheelAssemblyGyrostat {
 public:
   WheelAssemblyGyrostat() : Ixx(0.0), Iyy(0.0), Izz(0.0), Ixz(0.0),   J(0.0),
                               m(0.0),   R(0.0),   r(0.0),   a(0.0),   b(0.0),
                               c(0.0),  Tw(0.0),  Gx(0.0),  Gy(0.0),  Gz(0.0),
                              Fx(0.0),  Fy(0.0),  Fz(0.0),  Tx(0.0),  Ty(0.0),
                              Tz(0.0) {}
   double Ixx, Iyy, Izz, Ixz,
          J, m, R, r, a, b, c,
          Tw, Tx, Ty, Tz, Gx, Gy, Gz, Fx, Fy, Fz;
};
#endif
