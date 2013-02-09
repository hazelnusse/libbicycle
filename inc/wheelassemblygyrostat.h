#ifndef WHEELASSEMBLYGYROSTAT_H
#define WHEELASSEMBLYGYROSTAT_H
#include <iostream>

class WheelAssemblyGyrostat {
 public:
   WheelAssemblyGyrostat();
   friend std::ostream & operator<<(std::ostream & os,
                                    const WheelAssemblyGyrostat & w);

   double Ixx, Iyy, Izz, Ixz,
          J, m, R, r, a, b, c,
          Tw, Tx, Ty, Tz, Gx, Gy, Gz, Fx, Fy, Fz;
};
#endif
