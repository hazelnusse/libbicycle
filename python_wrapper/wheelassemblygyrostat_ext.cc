#include <boost/python.hpp>

#include "wheelassemblygyrostat.h"

BOOST_PYTHON_MODULE(wheelassemblygyrostat_ext)
{
    using namespace boost::python;
    using bicycle::WheelAssemblyGyrostat;
    class_<WheelAssemblyGyrostat>("WheelAssemblyGyrostat")
      .def_readwrite("Ixx", &WheelAssemblyGyrostat::Ixx)
      .def_readwrite("Iyy", &WheelAssemblyGyrostat::Iyy)
      .def_readwrite("Izz", &WheelAssemblyGyrostat::Izz)
      .def_readwrite("Ixz", &WheelAssemblyGyrostat::Ixz)
      .def_readwrite("J", &WheelAssemblyGyrostat::J)
      .def_readwrite("m", &WheelAssemblyGyrostat::m)
      .def_readwrite("R", &WheelAssemblyGyrostat::R)
      .def_readwrite("r", &WheelAssemblyGyrostat::r)
      .def_readwrite("a", &WheelAssemblyGyrostat::a)
      .def_readwrite("b", &WheelAssemblyGyrostat::b)
      .def_readwrite("c", &WheelAssemblyGyrostat::c)
      .def_readwrite("Tw", &WheelAssemblyGyrostat::Tw)
      .def_readwrite("Tx", &WheelAssemblyGyrostat::Tx)
      .def_readwrite("Ty", &WheelAssemblyGyrostat::Ty)
      .def_readwrite("Tz", &WheelAssemblyGyrostat::Tz)
      .def_readwrite("Gx", &WheelAssemblyGyrostat::Gx)
      .def_readwrite("Gy", &WheelAssemblyGyrostat::Gy)
      .def_readwrite("Gz", &WheelAssemblyGyrostat::Gz)
      .def_readwrite("Fx", &WheelAssemblyGyrostat::Fx)
      .def_readwrite("Fy", &WheelAssemblyGyrostat::Fy)
      .def_readwrite("Fz", &WheelAssemblyGyrostat::Fz)
      .def(self_ns::str(self));
}
 
