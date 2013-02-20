from wheelassemblygyrostat import WheelAssemblyGyrostat
from sympy import symbols
from sympy.physics.mechanics import dynamicsymbols
from math import pi
rear = WheelAssemblyGyrostat('rear')
front = WheelAssemblyGyrostat('front')
ls, g, t = symbols('ls_ g_ t')
q = dynamicsymbols('q:8')

gyrostat_benchmark_parameters = {
    #Rear assembly:
    rear.Ixx :         7.684799791449106,
    rear.Iyy :         11.99931034482758,
    rear.Izz :         5.315110553378478,
    rear.Ixz :         4.262158094617231,
    rear.J   :                      0.12,
    rear.m   :                        87,
    rear.R   :                       0.3,
    rear.r   :                         0,
    rear.a   :        0.4599058376856175,
    rear.b   :       -0.4669419422355363,
    rear.c   :        0.9534570696121847,
    #Front assembly:
    front.Ixx :        0.4335379755311007,
    front.Iyy :        0.5746857142857142,
    front.Izz :        0.1481477387546135,
    front.Ixz :      0.005332503757935522,
    front.J   :                      0.28,
    front.m   :                         7,
    front.R   :                      0.35,
    front.r   :                         0,
    front.a   :     -0.003411905099535387,
    front.b   :       -0.2114010400161699,
    front.c   :       -0.0320714267276193,
    ls  :        0.2676445084476887,
    g   :                      9.81}

gyrostat_benchmark_reference_configuration = {q[1]:0.0, q[2]:pi/10.0, q[3]:0.0}

gyrostat_benchmark_external_forces = {
        rear.Tw : 0.0,
        rear.Tx : 0.0,
        rear.Ty : 0.0,
        rear.Tz : 0.0,
        rear.Fx : 0.0,
        rear.Fy : 0.0,
        rear.Fz : 0.0,
        front.Tw : 0.0,
        front.Tx : 0.0,
        front.Ty : 0.0,
        front.Tz : 0.0,
        front.Fx : 0.0,
        front.Fy : 0.0,
        front.Fz : 0.0,
        g : 9.81}
