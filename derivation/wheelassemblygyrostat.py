from sympy import Symbol

class WheelAssemblyGyrostatParams(object):
    """A class to hold constant physical parameters which describe a wheel
    assembly gyrostat."""

    def __init__(self, name):
        def nonnegative(s, name):
            return Symbol(s + '_' + name + '_param', real=True, nonnegative=True)

        def real(s, name):
            return Symbol(s + '_' + name + '_param', real=True)

        self.Ixx = nonnegative('Ixx', name)   # x * I * x
        self.Iyy = nonnegative('Iyy', name)   # y * I * y
        self.Izz = nonnegative('Izz', name)   # z * I * z
        self.Ixy = real('Ixy', name)          # x * I * y
        self.Iyz = real('Iyz', name)          # y * I * z
        self.Ixz = real('Ixz', name)          # x * I * z
        self.J = nonnegative('J', name)       # Wheel spin inertia
        self.m = nonnegative('m', name)       # Gyrostat total mass
        self.R = nonnegative('R', name)       # Wheel radius
        self.r = nonnegative('r', name)       # Tire radius
        self.a = real('a', name)              # r^{mc/wc} * x
        self.b = real('b', name)              # r^{mc/wc} * y
        self.c = real('c', name)              # r^{mc/wc} * z
        self.d = real('d', name)              # r^{sa/wc} * y

    def __str__(self):
        return self.name + "_"

class WheelAssemblyGyrostatInputs(object):
    """A class to hold inputs (forces and/or torques) that are applied to the
    gyrostat."""

    def __init__(self, name):
        def nonnegative(s, name):
            return Symbol(s + '_' + name + '_input', real=True, nonnegative=True)

        def real(s, name):
            return Symbol(s + '_' + name + '_input', real=True)

        self.Tw = real('Tw', name)            # Wheel torque
        self.Tx = real('Tx', name)            # T * x (applied torque)
        self.Ty = real('Ty', name)            # T * y (applied torque)
        self.Tz = real('Tz', name)            # T * z (applied torque)
        self.Gx = real('Gx', name)            # G * x (ground contact force)
        self.Gy = real('Gy', name)            # G * y (ground contact force)
        self.Gz = real('Gz', name)            # G * z (ground contact force)
        self.Fx = real('Fx', name)            # F * x (force at COM)
        self.Fy = real('Fy', name)            # F * y (force at COM)
        self.Fz = real('Fz', name)            # F * z (force at COM)

    def __str__(self):
        return self.name + "_"
