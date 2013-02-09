"""
Convert the 25 parameters used in Meijaard et al. 2007 to the 23 used in my
model which employs a front and rear gyrostat and allows for arbitrary frame
inertia and mass center location.  In addition, the Meijaard model is extended
slightly to allow for toroidal tires, so a total of 27 parameters need to be
specified for the Meijaard model.

Six parameters are identical:
    Rear wheel assembly:    J, R, r     ==  IRyy, rR, tR
    Front wheel assembly:   J, R, r     ==  IFyy, rF, tF

Two parameters are trivial:
    Rear wheel assembly:    m           ==  mR + mB
    Front wheel assembly:   m           ==  mF + mH

The remaining 15 (23 - 16) parameters that need to be calculated are:
    Rear wheel assembly:    Ixx, Iyy, Izz, Ixz, a, b, c
    Front wheel assembly:   Ixx, Iyy, Izz, Ixz, a, b, c
    Steer axis offset:      ls
"""

from sympy import symbols, ccode, simplify, sin, cos, numbered_symbols, cse
from sympy.physics.mechanics import (Vector, ReferenceFrame, Point,
        inertia_of_point_mass, inertia)
Vector.simp = False
from wheelassemblygyrostat import WheelAssemblyGyrostat
import re

### Parameters used in Gyrostat formulation
Rear = WheelAssemblyGyrostat('r')
Front = WheelAssemblyGyrostat('f')

### Parameters used in Meijaard formulation
w, c, l = symbols('w.w w.c w.lambda')
rR, rF, tR, tF = symbols('w.rR w.rF w.tR w.tF')
mR, mF, mH, mB = symbols('w.mR w.mF w.mH w.mB')
IRxx, IRyy = symbols('w.IRxx w.IRyy')
IBxx, IByy, IBzz, IBxz = symbols('w.IBxx w.IByy w.IBzz w.IBxz')
IHxx, IHyy, IHzz, IHxz = symbols('w.IHxx w.IHyy w.IHzz w.IHxz')
IFxx, IFyy = symbols('w.IFxx w.IFyy')
xB, zB, xH, zH = symbols('w.xB w.zB w.xH w.zH')

# Orient reference frames
N = ReferenceFrame('N')                 # Newtonian frame
R = N.orientnew('R', 'Axis', [l, N.y])  # Rear frame, z aligned with steer axis

# Location mass centers using Meijaard parameters
P = Point('P')                                  # Rear contact
RO = P.locatenew('RO', -(rR + tR)*N.z)          # Rear wheel center
BO = P.locatenew('BO', xB*N.x + zB*N.z)         # Rear frame center
HO = P.locatenew('HO', xH*N.x + zH*N.z)         # Front frame center
FO = P.locatenew('FO', w*N.x - (rF + tF)*N.z)   # Front wheel center
Q = P.locatenew('Q', w*N.x)

# Mass centers of rear and front assembly relative to respective wheel centers
RO_BO_mc = RO.locatenew('RO_BO_mc', (mR*RO.pos_from(RO) +
    mB*BO.pos_from(RO)).express(R)/(mR + mB))
FO_HO_mc = FO.locatenew('FO_HO_mc', (mF*FO.pos_from(FO) +
    mH*HO.pos_from(FO)).express(R)/(mF + mH))

# Inertia dyads of bodies in Meijaard et al. 2007
IR_RO = inertia(R, IRxx, IRyy, IRxx)
IB_BO = inertia(N, IBxx, IByy, IBzz, 0, 0, IBxz)
IH_HO = inertia(N, IHxx, IHyy, IHzz, 0, 0, IHxz)
IF_FO = inertia(R, IFxx, IFyy, IFxx)

IB_BO = (R.x & IB_BO & R.x).expand()*(R.x|R.x) +\
        (R.y & IB_BO & R.y).expand()*(R.y|R.y) +\
        (R.z & IB_BO & R.z).expand()*(R.z|R.z) +\
        (R.x & IB_BO & R.z).expand()*(R.x|R.z) +\
        (R.z & IB_BO & R.x).expand()*(R.z|R.x)
IH_HO = (R.x & IH_HO & R.x).expand()*(R.x|R.x) +\
        (R.y & IH_HO & R.y).expand()*(R.y|R.y) +\
        (R.z & IH_HO & R.z).expand()*(R.z|R.z) +\
        (R.x & IH_HO & R.z).expand()*(R.x|R.z) +\
        (R.z & IH_HO & R.x).expand()*(R.z|R.x)

# Inertia of rear and front assemblies, expressed in R frame (same as F frame when
# steer is zero, as is the case in reference configuration)
IRear = IR_RO + IB_BO +\
        inertia_of_point_mass(mR, RO.pos_from(RO_BO_mc), R) +\
        inertia_of_point_mass(mB, BO.pos_from(RO_BO_mc), R)

IFront = IF_FO + IH_HO +\
         inertia_of_point_mass(mF, FO.pos_from(FO_HO_mc), R) +\
         inertia_of_point_mass(mH, HO.pos_from(FO_HO_mc), R)


expressions = [
    ("rear_.Ixx", R.x & IRear & R.x),
    ("rear_.Iyy", R.y & IRear & R.y),
    ("rear_.Izz", R.z & IRear & R.z),
    ("rear_.Ixz", R.x & IRear & R.z),
    ("rear_.J", IRyy), 
    ("rear_.m", mR + mB),
    ("rear_.R", rR),
    ("rear_.r", tR),
    ("rear_.a", RO_BO_mc.pos_from(RO) & R.x),
    ("rear_.b", RO_BO_mc.pos_from(RO) & R.z),
    ("rear_.c", (w + c)*cos(l) - (rR + tR)*sin(l)), # manual calculation
    ("front_.Ixx", R.x & IFront & R.x),
    ("front_.Iyy", R.y & IFront & R.y),
    ("front_.Izz", R.z & IFront & R.z),
    ("front_.Ixz", R.x & IFront & R.z),
    ("front_.J", IFyy),
    ("front_.m", mH + mF),
    ("front_.R", rF),
    ("front_.r", tF),
    ("front_.a", FO_HO_mc.pos_from(FO) & R.x),
    ("front_.b", FO_HO_mc.pos_from(FO) & R.z),
    ("front_.c", c*cos(l) - (rF + tF)*sin(l)), # manual calculation
    ("ls_", w*sin(l) + (rR + tR - rF - tF)*cos(l))]

eqns_to_cse = [rhs for lhs, rhs in expressions]
repl, redu = cse(eqns_to_cse, symbols=numbered_symbols("z"))

s  = '#include "bicycle.h"\n'
s += '#include "whipple.h"\n\n'
s += 'void Bicycle::set_parameters_from_whipple(const Whipple & w) {\n'
s += '  double * z = new double[{0}];\n\n'.format(len(repl))
for i, r in enumerate(repl):
    s += "  " + re.sub(r'z(\d+)', r'z[\1]', str(r[0])) + " = "
    s += re.sub(r'z(\d+)', r'z[\1]', ccode(r[1])) + ";\n"

s += "\n"
for i, (ex_i, red_i) in enumerate(zip(expressions, redu)):
    s += "  " + ex_i[0] + " = "
    s += re.sub(r'z(\d+)', r'z[\1]', ccode(red_i))
    s += ";\n"
s += "\n  delete [] z;\n}\n\n"

f = open("bicycle_convert_whipple.cc", "w")
f.write(s)
f.close()
