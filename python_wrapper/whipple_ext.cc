#include <boost/python.hpp>

#include "whipple.h"

BOOST_PYTHON_MODULE(whipple_ext)
{
    using namespace boost::python;
    class_<bicycle::Whipple>("Whipple")
      .def_readwrite("w", &bicycle::Whipple::w)
      .def_readwrite("c", &bicycle::Whipple::c)
      .def_readwrite("lambda", &bicycle::Whipple::lambda)
      .def_readwrite("g", &bicycle::Whipple::g)
      .def_readwrite("rR", &bicycle::Whipple::rR)
      .def_readwrite("mR", &bicycle::Whipple::mR)
      .def_readwrite("IRxx", &bicycle::Whipple::IRxx)
      .def_readwrite("IRyy", &bicycle::Whipple::IRyy)
      .def_readwrite("xB", &bicycle::Whipple::xB)
      .def_readwrite("zB", &bicycle::Whipple::zB)
      .def_readwrite("mB", &bicycle::Whipple::mB)
      .def_readwrite("IBxx", &bicycle::Whipple::IBxx)
      .def_readwrite("IByy", &bicycle::Whipple::IByy)
      .def_readwrite("IBzz", &bicycle::Whipple::IBzz)
      .def_readwrite("IBxz", &bicycle::Whipple::IBxz)
      .def_readwrite("xH", &bicycle::Whipple::xH)
      .def_readwrite("zH", &bicycle::Whipple::zH)
      .def_readwrite("mH", &bicycle::Whipple::mH)
      .def_readwrite("IHxx", &bicycle::Whipple::IHxx)
      .def_readwrite("IHyy", &bicycle::Whipple::IHyy)
      .def_readwrite("IHzz", &bicycle::Whipple::IHzz)
      .def_readwrite("IHxz", &bicycle::Whipple::IHxz)
      .def_readwrite("rF", &bicycle::Whipple::rF)
      .def_readwrite("mF", &bicycle::Whipple::mF)
      .def_readwrite("IFxx", &bicycle::Whipple::IFxx)
      .def_readwrite("IFyy", &bicycle::Whipple::IFyy)
      .def_readwrite("tR", &bicycle::Whipple::tR)
      .def_readwrite("tF", &bicycle::Whipple::tF);
}
