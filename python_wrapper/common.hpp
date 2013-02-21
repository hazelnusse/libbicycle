#ifndef LIBTHETA__COMMON_H
#define LIBTHETA__COMMON_H

// Don't ask.
#define PY_ARRAY_UNIQUE_SYMBOL bicycle_python_PyArray_API
#ifndef SKIP_NO_IMPORT
#define NO_IMPORT
#endif


#include <Python.h>
#define NPY_NO_DEPRECATED_API 8
#include <numpy/arrayobject.h>
#include <Eigen/Dense>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <cerrno>
#include <string>
#include <iostream>

#define UNLIKELY(x) __builtin_expect(!!(x), 0)

#ifdef NDEBUG
#define ENSURE(x) do { if(UNLIKELY(!(x))) throw std::runtime_error("FATAL: " #x); } while(false)
#else
#define ENSURE(x) assert(x)
#endif

namespace bicycle_python {

typedef double Real;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMajorMatrix;
typedef RowMajorMatrix Matrix;
using namespace Eigen;

class ErrnoError : public std::runtime_error
{
public:
    ErrnoError(const std::string &prefix = "") :
        std::runtime_error(make_error_message(prefix, errno))
    {
    }

private:
    static std::string make_error_message(const std::string &prefix, int errnum)
    {
        std::string msg;
        if(!prefix.empty()) msg = prefix + ": ";
        msg.append(strerror(errnum));
        return msg;
    }
};

class PythonError : public std::exception
{
public:
    PythonError()
    {
        ENSURE(PyErr_Occurred());
    }

    PythonError(PyObject *type, const std::string &msg)
    {
        PyErr_SetString(type, msg.c_str());
    }
};

// note: return values of const methods are non-const because Python/C API is not const-aware
template <typename PyObjectT>
class PyRef
{
    typedef void (PyRef<PyObjectT>::*dummy_bool_type)() const;
    void dummy_true() const {}

public:
    static PyRef newref(PyObjectT *obj)
    {
        return PyRef(obj);
    }

    static PyRef borrow(PyObjectT *obj)
    {
        Py_XINCREF(obj);
        return PyRef(obj);
    }

    PyRef() :
        obj(NULL)
    {
    }

    PyRef(decltype(nullptr)) :
        obj(NULL)
    {
    }

    PyRef(const PyRef &other) :
        obj(other.obj)
    {
        Py_XINCREF(obj);
    }

    PyRef(PyRef &&other) :
        obj(NULL)
    {
        std::swap(obj, other.obj);
    }

    ~PyRef()
    {
        Py_XDECREF(obj);
    }

    PyRef& operator= (const PyRef &other)
    {
        Py_XINCREF(other.obj);
        Py_XDECREF(obj);
        obj = other.obj;
        return *this;
    }

    PyRef& operator= (PyRef &&other)
    {
        std::swap(obj, other.obj);
    }

    PyObjectT* newref() const
    {
        Py_XINCREF(obj);
        return obj;
    }

    PyObjectT* borrow() const
    {
        return obj;
    }

    PyObjectT& operator* () const
    {
        return *obj;
    }

    PyObjectT* operator-> () const
    {
        return obj;
    }

    template <typename TargetT>
    PyRef<TargetT> cast() const
    {
        return PyRef<TargetT>::borrow(reinterpret_cast<TargetT*>(obj));
    }

    // uses the "safe bool" idiom: http://www.artima.com/cppsource/safebool.html
    operator dummy_bool_type() const
    {
        return obj ? &PyRef::dummy_true : 0;
    }

private:
    PyRef(PyObjectT *obj_) :
        obj(obj_)
    {
    }

    PyObjectT *obj;
};

template <typename PyObjectT>
PyRef<PyObjectT> newref(PyObjectT *obj)
{
    return PyRef<PyObjectT>::newref(obj);
}

template <typename PyObjectT>
PyRef<PyObjectT> borrow(PyObjectT *obj)
{
    return PyRef<PyObjectT>::borrow(obj);
}

typedef PyRef<PyObject> PyObjectRef;
typedef PyRef<PyArrayObject> PyArrayRef;
typedef PyObjectRef PyFunc(const PyObjectRef &self, const PyObjectRef &args);
typedef PyObjectRef PyFuncWithKeywords(const PyObjectRef &self, const PyObjectRef &args, const PyObjectRef &kwargs);

template <PyFunc func>
PyObject* pywrap(PyObject *self, PyObject *args)
{
    try {
        return func(borrow(self), borrow(args)).newref();
    } catch(const PythonError &e) {
        assert(PyErr_Occurred());
        return NULL;
    } catch(const std::exception &e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    } catch(...) {
        PyErr_SetString(PyExc_RuntimeError, "<unknown error>");
        return NULL;
    }
}

template <PyFuncWithKeywords func>
PyObject* pywrap(PyObject *self, PyObject *args, PyObject *kwargs)
{
    try {
        return func(borrow(self), borrow(args), borrow(kwargs)).newref();
    } catch(const PythonError &e) {
        assert(PyErr_Occurred());
        return NULL;
    } catch(const std::exception &e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    } catch(...) {
        PyErr_SetString(PyExc_RuntimeError, "<unknown error>");
        return NULL;
    }
}

template <PyFuncWithKeywords func>
int pyinitproc(PyObject *self, PyObject *args, PyObject *kwargs)
{
    try {
        return func(borrow(self), borrow(args), borrow(kwargs)).newref() == NULL ? -1 : 0;
    } catch(const PythonError &e) {
        assert(PyErr_Occurred());
        return -1;
    } catch(const std::exception &e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return -1;
    } catch(...) {
        PyErr_SetString(PyExc_RuntimeError, "<unknown error>");
        return -1;
    }
}

int matrix_converter(PyObject *obj, void *dest);
int const_matrix_converter(PyObject *obj, void *dest);

Map<Vector, Aligned> pyvec(const PyArrayRef &obj);
Map<Matrix, Aligned> pymat(const PyArrayRef &obj);
Map<const Vector, Aligned> const_pyvec(const PyArrayRef &obj);
Map<const Matrix, Aligned> const_pymat(const PyArrayRef &obj);

inline PyObjectRef pyfloat(Real value)
{
    return newref(PyFloat_FromDouble(value));
}

template <typename Derived>
PyArrayRef pyvec(const DenseBase<Derived> &vec)
{
    ENSURE(vec.cols() == 1);
    npy_intp dims[] = { vec.size() };
    PyArrayRef result = newref(PyArray_EMPTY(1, dims, NPY_DOUBLE, 0)).template cast<PyArrayObject>();
    ENSURE(PyArray_Check(result.borrow()));
    pyvec(result) = vec;
    return result;
}

template <typename Derived>
PyArrayRef pymat(const MatrixBase<Derived> &mat)
{
    npy_intp dims[] = { mat.rows(), mat.cols() };
    PyArrayRef result = newref(PyArray_EMPTY(2, dims, NPY_DOUBLE, 0)).template cast<PyArrayObject>();
    ENSURE(PyArray_Check(result.borrow()));
    pymat(result) = mat;
    return result;
}

} // namespace libtheta

#endif
