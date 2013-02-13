#include <iomanip>
#include "wheelassemblygyrostat.h"

namespace bicycle {

WheelAssemblyGyrostat::WheelAssemblyGyrostat()
  : Ixx(0.0), Iyy(0.0), Izz(0.0), Ixz(0.0),   J(0.0),
      m(0.0),   R(0.0),   r(0.0),   a(0.0),   b(0.0),
      c(0.0),  Tw(0.0),  Tx(0.0),  Ty(0.0),  Tz(0.0),
     Gx(0.0),  Gy(0.0),  Gz(0.0),
     Fx(0.0),  Fy(0.0),  Fz(0.0)
{
}

std::ostream & operator<<(std::ostream & os,
                          const WheelAssemblyGyrostat & w)
{
  os.precision(16);
  os << "Ixx = " << std::setw(25) << w.Ixx << std::endl
     << "Iyy = " << std::setw(25) << w.Iyy << std::endl
     << "Izz = " << std::setw(25) << w.Izz << std::endl
     << "Ixz = " << std::setw(25) << w.Ixz << std::endl
     << "J   = " << std::setw(25) <<   w.J << std::endl
     << "m   = " << std::setw(25) <<   w.m << std::endl
     << "R   = " << std::setw(25) <<   w.R << std::endl
     << "r   = " << std::setw(25) <<   w.r << std::endl
     << "a   = " << std::setw(25) <<   w.a << std::endl
     << "b   = " << std::setw(25) <<   w.b << std::endl
     << "c   = " << std::setw(25) <<   w.c << std::endl;
  return os;
}

} // namespace bicycle
