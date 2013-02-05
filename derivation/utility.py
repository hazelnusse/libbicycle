from sympy import ccode, cse, numbered_symbols, symbols, sympify
import re
import numpy as np

class EigenMatrixCodeOutput(object):
    """A class to generate Eigen matrices from SymPy matrices.
    """

    def __init__(self):
        self.s = self._prefix()
        self.subs_dict = {}
        self.regex_list = []
        self.state_prefix = ''

    def generate(self, expressions, functionname):
        """Given a numpy nd-array of sympy expressions, return a function which
        will generate a C-style function that will compute the entries of the
        array.  The array is flattened out into a contiguous array with the
        last index varying fastest (row-major for matrices).
        """

        orig_shape = expressions.shape
        s =  "/*!\n"
        s += "   Computes the n-d array of shape ("
        for d in expressions.shape[:-1]:
            s += "{0}, ".format(d)
        s += "{0})\n\n".format(expressions.shape[-1])
        s += "   @param[out] a C-array of with {0} elements\n".format(expressions.size)
        s += "*/\n"

        expressions = np.array([sympify(ai).subs(self.subs_dict) for ai in expressions.flat])

        s += "void {0}(double m[{1}])".format(functionname, expressions.size)
        s += " {\n"
        repl, redu = cse(expressions.flat, symbols=numbered_symbols("z"))

        if len(repl):
            s += "  double * z = new double[{0}];\n\n".format(len(repl))

        for i, r in enumerate(repl):
            s += "  " + re.sub(r'z(\d+)', r'z[\1]', str(r[0])) + " = "
            tmp = re.sub(r'z(\d+)', r'z[\1]', ccode(r[1])) + ";\n"
            if self.state_prefix:
                tmp = re.sub(self.state_prefix + r'(\d+)',
                             self.state_prefix + r'[\1]', tmp)
            for p, r in self.regex_list:
                test = re.compile(p)
                tmp = test.sub(r, tmp)
            s += tmp

        s += "\n"
        for i, red_i in enumerate(redu):
            s += "  m[{0}] = ".format(i)
            tmp = re.sub(r'z(\d+)', r'z[\1]', ccode(redu[i]))
            if self.state_prefix:
                tmp = re.sub(self.state_prefix + r'(\d+)',
                             self.state_prefix + r'[\1]', tmp)
            for p, r in self.regex_list:
                test = re.compile(p)
                tmp = test.sub(r, tmp)

            s += tmp
            s += ";\n"

        s += "\n  delete [] z;\n}\n\n"

        self.s += s

        return s

    def output(self, filename):
        f = open(filename, "w")
        f.write(self.s)
        f.close()


    def _prefix(self):
        s = """#include <cmath>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

"""
        return s

    def set_states(self, variables, prefix):
        self.state_variables = variables
        self.state_prefix = prefix
        self.state_symbols = symbols(prefix + ":" + str(len(variables)))
        self.state_dict = {vi : si for vi, si in zip(variables, self.state_symbols)}
        self.subs_dict.update(self.state_dict)

    def add_regex(self, pattern, replacement):
        self.regex_list.append([pattern, replacement])

