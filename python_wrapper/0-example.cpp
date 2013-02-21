#define SKIP_NO_IMPORT
#include "common.hpp"
#include <iostream>

namespace bicycle_python {

PyObjectRef moments_from_samples(const PyObjectRef &self, const PyObjectRef &args)
{
    PyArrayObject *weights_ = NULL;
    PyArrayObject *samples_ = NULL;

    if(!PyArg_ParseTuple(args.borrow(), "O&O&", const_matrix_converter, &weights_,
            const_matrix_converter, &samples_)) {
        return nullptr;
    }

    const auto &weights = const_pyvec(newref(weights_));
    const auto &samples = const_pymat(newref(samples_));
    const unsigned dim = samples.rows();
    const unsigned num_samples = samples.cols();
    ENSURE(weights.size() == num_samples);

    Real total_weight = 0.0;
    Vector mean = Vector::Zero(dim);
    Matrix covar_ = Matrix::Zero(dim, dim);
    auto covar = covar_.selfadjointView<Lower>();
    for(unsigned i = 0; i < num_samples; ++i) {
        const auto weight = weights[i];
        const auto sample = samples.col(i);
        mean += weight * sample;
        covar.rankUpdate(sample, weight);
        total_weight += weight;
    }
    mean /= total_weight;
    covar_ /= total_weight;
    covar.rankUpdate(mean, -1);

    const auto &o_mean = pyvec(mean);
    const auto &o_covar = pymat(Matrix(covar));
    return newref(PyTuple_Pack(2, o_mean.borrow(), o_covar.borrow()));
}

static PyMethodDef bicycle_python_methods[] = {
    { "moments_from_samples", pywrap<moments_from_samples>, METH_VARARGS, NULL },
    { 0 }
};

static struct PyModuleDef bicycle_python_module = {
   PyModuleDef_HEAD_INIT,
   "bicycle",   /* name of module */
   "Example numpy Eigen interface using Python C API.", /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   bicycle_python_methods
};


PyMODINIT_FUNC 
PyInit_bicycle_python()
{
    using namespace bicycle_python;
    import_array();
    PyObject * m = PyModule_Create(&bicycle_python_module);
    return m;
}

} // namespace bicycle_python
