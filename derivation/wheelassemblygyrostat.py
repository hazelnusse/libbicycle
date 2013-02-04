from sympy import Symbol

class WheelAssemblyGyrostat(object):
    """A class to hold physical parameters which describe a wheel assembly
    gyrostat."""
    
    def __init__(self, name):
        def nonnegative(s, name):
            return Symbol(s + '_' + name, real=True, nonnegative=True)
        
        def real(s, name):
            return Symbol(s + '_' + name, real=True)

        self.Ixx = nonnegative('Ixx', name)   # x * I * x
        self.Iyy = nonnegative('Iyy', name)   # y * I * y
        self.Izz = nonnegative('Izz', name)   # z * I * z
        self.Ixz = real('Ixz', name)          # x * I * z
        self.J = nonnegative('J', name)       # Wheel spin inertia
        self.m = nonnegative('m', name)       # Gyrostat total mass
        self.R = nonnegative('R', name)       # Wheel radius
        self.r = nonnegative('r', name)       # Tire radius
        self.a = real('a', name)              # r^{mc/wc} * x
        self.b = real('b', name)              # r^{mc/wc} * z
        self.c = real('c', name)              # r^{sa/wc} * x
