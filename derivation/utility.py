from sympy import ccode, cse, numbered_symbols, symbols
import re

class EigenMatrixCodeOutput(object):
    """A class to generate Eigen matrices from SymPy matrices.
    """

    def __init__(self):
        self.s = self._prefix()
        self.subs_dict = {}
        self.regex_list = []
        self.state_prefix = ''

    def generate(self, matrix, functionname):
        """Given a SymPy matrix, return a string of C++ code that will compile.

        If filename is given, the string will be written to the filename.
        """

        m, n = matrix.shape

        s = "Matrix<double, {0}, {1}> m {2}()".format(m, n, functionname)
        s += " {\n"
        repl, redu = cse(matrix.subs(self.subs_dict),
                         symbols=numbered_symbols("z"))

        s += "  Matrix<double, {0}, {1}> m;\n".format(m, n)
        if len(repl):
            s += "  Matrix<double, {0}, 1> z;\n".format(len(repl))
        for i, r in enumerate(repl):
            s += "  " + re.sub(r'z(\d+)', r'z(\1)', str(r[0])) + " = "
            tmp = re.sub(r'z(\d+)', r'z(\1)', ccode(r[1])) + ";\n"
            if self.state_prefix:
                tmp = re.sub(self.state_prefix + r'(\d+)',
                             self.state_prefix + r'(\1)', tmp)
            for p, r in self.regex_list:
                test = re.compile(p)
                tmp = test.sub(r, tmp)
            s += tmp

        s += "\n"
        for i in range(m):
            for j in range(n):
                s += "  m({0}, {1}) = ".format(i, j)
                tmp = re.sub(r'z(\d+)', r'z(\1)', ccode(redu[0][n*i + j]))
                if self.state_prefix:
                    tmp = re.sub(self.state_prefix + r'(\d+)',
                                 self.state_prefix + r'(\1)', tmp)
                for p, r in self.regex_list:
                    test = re.compile(p)
                    tmp = test.sub(r, tmp)

                s += tmp
                s += ";\n"

        s += "\n  return m;\n}\n\n"

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

    def add_subsitution_dict(self, d):
        self.subs_dict.update(d)

    def set_states(self, variables, prefix):
        self.state_variables = variables
        self.state_prefix = prefix
        self.state_symbols = symbols(prefix + ":" + str(len(variables)))
        self.state_dict = {vi : si for vi, si in zip(variables, self.state_symbols)}
        self.subs_dict.update(self.state_dict)

    def add_regex(self, pattern, replacement):
        self.regex_list.append([pattern, replacement])

