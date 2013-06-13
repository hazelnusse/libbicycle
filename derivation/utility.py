from sympy import ccode, cse, numbered_symbols, symbols, S
import re
import numpy as np

class NumpyArrayOutput(object):
    def __init__(self,
                 class_name=None,
                 includes=None,
                 enclosing_namespace=None,
                 using_namespace_directives=None,
                 using_declarations=None):
        """A class to generate flat C-arrays from numpy nd arrays.

        If you want certain files to be included or namespaces to be used, pass
        these in as a list of strings.  Include files need to be specified with
        " " or < > characters to indicate whether they are system includes or
        local includes.

        """

        if class_name is not None:
            self._class_name = class_name
        else:
            self._class_name = ""

        self.s = self._prefix(includes, enclosing_namespace,
                              using_namespace_directives, using_declarations)
        self.subs_dict = {}
        self.regex_list = []
        self.state_prefix = ''

    def generate(self, expressions, functionname, const_function=True, callbacks=None):
        """Given a numpy nd-array of sympy expressions, return a function which
        will generate a C-style function that will compute the entries of the
        array.  The array is flattened out into a contiguous array with the
        last index varying fastest (row-major for matrices).

        """

        orig_shape = expressions.shape
        s =  "\n/** Computes the n-d array of shape ("
        for d in expressions.shape[:-1]:
            s += "{0}, ".format(d)
        s += "{0})\n *\n".format(expressions.shape[-1])
        s += " * @param[out] ar a C-array of with {0}".format(expressions.size)
        s += " elements\n */\n"

        expressions_flat = np.zeros((expressions.size,), dtype=object)
        for i, ai in enumerate(expressions.flat):
            if isinstance(ai, int):
                expressions_flat[i] = S(0)
            else:
                expressions_flat[i] = ai.subs(self.subs_dict)

        if self._class_name != "":
            member_function_name = self._class_name + "::" + functionname
            function_signature = "void {0}(double ar[{1}])".format(member_function_name, expressions_flat.size)
        else:
            function_signature = "void {0}(double ar[{1}])".format(functionname, expressions_flat.size)

        if const_function:
            function_signature += " const"
        s += "//  " + function_signature + ";\n"
        s += function_signature + "\n{\n"

        repl, redu = cse(expressions_flat, symbols=numbered_symbols("z"))

        if len(repl):
            s += "  double z[{0}];\n\n".format(len(repl))

        if callbacks is not None:
            for c in callbacks:
                cname, symbolname = c
                s += "  " + cname + "(" + symbolname + ");\n"

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
            s += "  ar[{0}] = ".format(i)
            tmp = re.sub(r'z(\d+)', r'z[\1]', ccode(red_i))
            if self.state_prefix:
                tmp = re.sub(self.state_prefix + r'(\d+)',
                             self.state_prefix + r'[\1]', tmp)
            for p, r in self.regex_list:
                test = re.compile(p)
                tmp = test.sub(r, tmp)

            s += tmp
            s += ";\n"

        self.s += s + "}\n"

        return s

    def output(self, filename):
        """Write all generated code to file.  No checking is done on whether
        the file exists, so this will overwrite existing files.

        """

        if self._enclosing_namespace:
            self.s += "\n}\n"

        f = open(filename, "w")
        f.write(self.s)
        f.close()

    def _prefix(self, includes, enclosing_namespace,
                using_namespace_directives, using_declarations):
        s = ""

        if includes is not None:
            for i in includes:
                s += "#include " + i + "\n"

        if enclosing_namespace is not None:
            s += "\nnamespace " + enclosing_namespace + " {\n"
            self._enclosing_namespace = True
        else:
            self._enclosing_namespace = False

        if using_namespace_directives is not None:
            s += "\n"
            for n in using_namespace_directives:
                s += "using namespace " + n + ";\n"
        
        if using_declarations is not None:
            s += "\n"
            for n in using_declarations:
                s += "using " + n + ";\n"

        return s

    def set_states(self, variables, prefix):
        """Define state variables.

        The generated code is output so that all symbols in the variables
        argument are rewritten to be indexed into a state array with the name
        specified by prefix.

        """
        self.state_variables = variables
        self.state_prefix = prefix
        self.state_symbols = symbols(prefix + ":" + str(len(variables)))
        self.state_dict = {vi : si for vi, si in zip(variables, self.state_symbols)}
        self.subs_dict.update(self.state_dict)

    def add_regex(self, pattern, replacement):
        """Add a pattern and replacement that will be applied to ccode output
        strings.

        This is useful if you want to use different variable names in your
        C/C++ code than what you use in your Python code.
        """
        self.regex_list.append([pattern, replacement])

