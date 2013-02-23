#include <boost/python.hpp>

#include "bicycle.h"

BOOST_PYTHON_MODULE(bicycle_ext)
{
    using namespace boost::python;
    using namespace bicycle;
    void (Bicycle::*set_parameters_whipple)(const Whipple &) = &Bicycle::set_parameters;
    void (Bicycle::*set_parameters_gyrostats)
      (const WheelAssemblyGyrostat&, const WheelAssemblyGyrostat&,
       double, double) = &Bicycle::set_parameters;

    class_<Bicycle>("Bicycle")
      .def("set_state", &Bicycle::set_state)
      .def("set_parameters", set_parameters_whipple)
      .def("set_parameters", set_parameters_gyrostats)
      .def(self_ns::str(self));
}
 
