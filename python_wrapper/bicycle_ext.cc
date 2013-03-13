#include <Python.h>
#include <iostream>
#include <utility>
#include <boost/python.hpp>
#define NPY_NO_DEPRECATED_API 8
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include "bicycle.h"

namespace bp = boost::python;
using bicycle::Bicycle;
using bicycle::Vector;
using bicycle::Whipple;
using bicycle::WheelAssemblyGyrostat;

void Bicycle::solve_configuration_constraint_and_set_state_python()
{
  solve_configuration_constraint_and_set_state();
}

void Bicycle::solve_velocity_constraints_and_set_state_python()
{
  solve_velocity_constraints_and_set_state();
}

void Bicycle::mass_matrix_full_python(PyObject * matrix_out) const
{
  const int N = Bicycle::n + Bicycle::o;
  typedef Eigen::Map<Eigen::Matrix<double, N, N, Eigen::RowMajor>> MapMatrix;
  MapMatrix mm((double *) PyArray_DATA((PyArrayObject *)matrix_out));
  mm = mass_matrix_full();
}

void Bicycle::independent_state_matrix_python(PyObject * matrix_out) const
{
  const int N = Bicycle::n + Bicycle::o;
  const int M = Bicycle::n - Bicycle::l + Bicycle::o - Bicycle::m;
  typedef Eigen::Map<Eigen::Matrix<double, N, M, Eigen::RowMajor>> MapMatrix;
  MapMatrix mm((double *) PyArray_DATA((PyArrayObject *)matrix_out));
  mm = independent_state_matrix();
}

void Bicycle::input_matrix_python(PyObject * matrix_out) const
{
  const int N = Bicycle::o - Bicycle::m;
  const int S = Bicycle::s;
  typedef Eigen::Map<Eigen::Matrix<double, N, S, Eigen::RowMajor>> MapMatrix;
  MapMatrix mm((double *) PyArray_DATA((PyArrayObject *)matrix_out));
  mm = input_matrix();
}

void set_steady_constraint_forces(Bicycle * b)
{
  Vector r = Vector::Zero(22);
  Vector cf = b->steady_no_slip_constraint_forces();
  r[4] = cf[0];
  r[5] = cf[1];
  r[6] = cf[2];
  r[14] = cf[3];
  r[15] = cf[4];
  r[16] = cf[5];
  r[20] = cf[6];
  r[21] = 9.81;
  b->set_inputs(r);
}

BOOST_PYTHON_MODULE(bicycle_ext)
{
    void (Bicycle::*set_parameters_whipple)(const Whipple &) = &Bicycle::set_parameters;
    void (Bicycle::*set_parameters_gyrostats)
      (const WheelAssemblyGyrostat&, const WheelAssemblyGyrostat&,
       double, double) = &Bicycle::set_parameters;

    bp::class_<Bicycle>("Bicycle")
      .def_readonly("n", &Bicycle::n)
      .def_readonly("l", &Bicycle::l)
      .def_readonly("n_min", &Bicycle::n_min)
      .def_readonly("o", &Bicycle::o)
      .def_readonly("m", &Bicycle::m)
      .def_readonly("s", &Bicycle::s)
      .def("set_state", &Bicycle::set_state)
      .def("set_parameters", set_parameters_whipple)
      .def("set_parameters", set_parameters_gyrostats)
      .def("solve_configuration_constraint_and_set_state",
          &Bicycle::solve_configuration_constraint_and_set_state_python)
      .def("solve_velocity_constraints_and_set_state",
          &Bicycle::solve_velocity_constraints_and_set_state_python)
      .def("set_steady_constraint_forces", set_steady_constraint_forces)
      .def("mass_matrix_full", &Bicycle::mass_matrix_full_python)
      .def("independent_state_matrix", &Bicycle::independent_state_matrix_python)
      .def("input_matrix", &Bicycle::input_matrix_python)
      .def(bp::self_ns::str(bp::self));
}
