#include "common.hpp"

namespace bicycle_python {

Map<Vector, Aligned> pyvec(const PyArrayRef &ref)
{
    auto obj = ref.borrow();
    if(!(PyArray_FLAGS(obj) & NPY_ARRAY_CARRAY)) throw std::runtime_error("unmappable ndarray");
    if(PyArray_TYPE(obj) != NPY_DOUBLE) throw std::runtime_error("ndarray is not of type double");
    if(PyArray_NDIM(obj) != 1) throw std::runtime_error("expected a vector");
    double *data = reinterpret_cast<double*>(PyArray_DATA(obj));
    const auto rows = PyArray_DIM(obj, 0);
    return Map<Vector, Aligned>(data, rows);
}

Map<const Vector, Aligned> const_pyvec(const PyArrayRef &ref)
{
    auto obj = ref.borrow();
    if(!(PyArray_FLAGS(obj) & NPY_ARRAY_CARRAY_RO)) throw std::runtime_error("unmappable ndarray");
    if(PyArray_TYPE(obj) != NPY_DOUBLE) throw std::runtime_error("ndarray is not of type double");
    if(PyArray_NDIM(obj) != 1) throw std::runtime_error("expected a vector");
    const double *data = reinterpret_cast<const double*>(PyArray_DATA(obj));
    const auto rows = PyArray_DIM(obj, 0);
    return Map<const Vector, Aligned>(data, rows);
}

Map<Matrix, Aligned> pymat(const PyArrayRef &ref)
{
    auto obj = ref.borrow();
    if(!(PyArray_FLAGS(obj) & NPY_ARRAY_CARRAY)) throw std::runtime_error("unmappable ndarray");
    if(PyArray_TYPE(obj) != NPY_DOUBLE) throw std::runtime_error("ndarray is not of type double");
    if(PyArray_NDIM(obj) != 2) throw std::runtime_error("expected a matrix");
    double *data = reinterpret_cast<double*>(PyArray_DATA(obj));
    const auto rows = PyArray_DIM(obj, 0);
    const auto cols = PyArray_DIM(obj, 1);
    return Map<Matrix, Aligned>(data, rows, cols);
}

Map<const Matrix, Aligned> const_pymat(const PyArrayRef &ref)
{
    auto obj = ref.borrow();
    if(!(PyArray_FLAGS(obj) & NPY_ARRAY_CARRAY_RO)) throw std::runtime_error("unmappable ndarray");
    if(PyArray_TYPE(obj) != NPY_DOUBLE) throw std::runtime_error("ndarray is not of type double");
    if(PyArray_NDIM(obj) != 2) throw std::runtime_error("expected a matrix");
    const double *data = reinterpret_cast<const double*>(PyArray_DATA(obj));
    const auto rows = PyArray_DIM(obj, 0);
    const auto cols = PyArray_DIM(obj, 1);
    return Map<const Matrix, Aligned>(data, rows, cols);
}

int const_matrix_converter(PyObject *obj, void *dest_)
{
    PyObject **dest = reinterpret_cast<PyObject**>(dest_);
    *dest = PyArray_FromAny(obj, PyArray_DescrFromType(NPY_DOUBLE), 1, 2, NPY_ARRAY_IN_ARRAY, NULL);
    return 1;
}

int matrix_converter(PyObject *obj, void *dest_)
{
    PyObject **dest = reinterpret_cast<PyObject**>(dest_);
    *dest = PyArray_FromAny(obj, PyArray_DescrFromType(NPY_DOUBLE), 1, 2, NPY_ARRAY_INOUT_ARRAY, NULL);
    return 1;
}

} // namespace bicycle_python
