# # Derivation of equations of motion for an extended bicycle model
#
# This notebook derives symbolic equations of motion, along with a number of
# related quantities and functions that are useful for the analysis of the
# bicycle.  I refer to it as "extended" because it extends the Whipple bicycle
# model in four ways:  wheels with elliptical cross sections instead of
# knife-edged wheels, arbitrary frame and fork inertia (lateral symmetry is not
# assumed), arbitrary frame and fork mass center location, wheel axle
# misalignment angles (axle not assumed perpendicular to the traditional frame
# plane), and steer axis offsets which allow for wheels to not be in the same
# plane even when the steer angle is zero.  The contact model is such that
# unspecified forces are applied in the lateral and longitudinal directions of
# each wheel contact point; this approach permits either no-slip rolling
# assumptions to be made and allows for constraint force determination, or a
# tire contact force model can be added.  In the former case, the model has 3
# degrees of freedom; in the latter case, it has 7 degrees of freedom.
from sympy import symbols, Matrix, zeros, sqrt, trigsimp
from sympy.physics.mechanics import *
Vector.simp = False
from wheelassemblygyrostat import WheelAssemblyGyrostat
from utility import EigenMatrixCodeOutput

# Output code generator
code = EigenMatrixCodeOutput()

# ## Model description
# ### Parameters
# The model takes a front/rear symmetric view of the bicycle.  Eighteen
# parameters are defined which each have a front and rear analog. Geometric and
# inertial axial symmetry is assumed for both the front and the rear wheels.
# This symmetry assumption implies implies that the rear frame and wheel, as a
# pair, can be viewed as a cylindrical gyrostat: the mass center and combined
# inertia of the two bodies does not change as the wheel spins.  The same
# implication is true of the front fork and front wheel.  For brevity, I refer
# to the rear frame and rear wheel as the rear assembly; similarly, the front
# fork and front wheel are the front assembly.  While it may be helpful to
# think of a bicycle with a standard frame and fork, this model can represent
# other types of two wheeled vehicles such as Segways, Razor scooters, or other
# non-traditional two wheeled vehicles.
#
# Parameters associated with the rear assembly are subscripted with a 'r',
# those associated with the front assembly with a 'f'.  For example $m_r$ for
# the rear assembly mass, $m_f$ for the front assembly mass.  The parameters
# are:
#
# * $A, B, C, D, E, F$:  Inertia scalars of combined wheel and frame.
# * $J$: spin inertia of wheel.
# * $m$: combined mass of wheel and frame.
# * $r, a, b$: Wheel major radii, tire cross section major and minor radii (wheel is modelled as a revolved ellipse).
# * $c, d, e$: distances from wheel center to frame and wheel mass center.
# * $f, g, h$: distances from wheel center to steer axis point
# * $i, j$: axle misalignment angles; first about frame fixed $y$, then about rotated $x$.
#
# I manage these parameters by declaring two instances of the
# `WheelAssemblyGyrostat` class:

rear = WheelAssemblyGyrostat('rear')
front = WheelAssemblyGyrostat('front')
code.add_regex(r'([_0-9a-zA-Z]+)_(rear|front)', r'\2.\1')
#code.add_regex(r'(.+)_front', r'front.\1')

# Gravitational constant, time
ls, g, t = symbols('ls g t')

# ### Generalized coordinates
#
# The choice of coordinates are identical as to those used by [Meeijaard et
# al.][Meijaard2007], with one exception.  The rear frame pitch ($\theta_B$ in
# [Meeijaard et al.][Meijaard2007]) differs by a constant: $q_2 = \theta_B +
# \lambda$, where $q_2$ is the coordinate I use to denote the pitch.  This
# difference is due to the fact that the rear assembly body fixed frame has
# axes aligned with the steer axis and wheel symmetry axes, rather than being
# aligned with the inertial axes when the bicycle is in the "reference"
# configuration (zero lean and steer).
#
# I follow the convention of Thomas Kane of using $q$ to denote generalized
# coordinates, except that I begin with $q_0$ rather than $q_1$; this maps more
# naturally to languages which begin indices at 0 instead of 1, such as Python,
# C, and C++.
#
# The generalized coordinates are:
#
# * $q_0$:  Rear frame yaw angle
# * $q_1$:  Rear frame lean angle
# * $q_2$:  Rear frame pitch angle
# * $q_3$:  Steer angle
# * $q_4$:  Rear wheel angle
# * $q_5$:  Front wheel angle
# * $q_6$:  Rear wheel ground contact in $\mathbf{\hat{n}}_x$ direction
# * $q_7$:  Rear wheel ground contact in $\mathbf{\hat{n}}_y$ direction
# [Meijaard2007]: http://dx.doi.org/10.1098/rspa.2007.1857

q = dynamicsymbols('q:8')
qd = [qi.diff(t) for qi in q]

# ### Generalized speeds
#
# A total of 12 generalized speeds are used to describe the angular velocities
# of each body and the velocities of each mass center.  Six of these
# generalized speeds are associated with angular velocity are defined as $u_i
# \triangleq \dot{q}_i$ for $i = 0,\dots,5$.  I experimented with using the
# choice of generalized speeds advocated by [Mitiguy and Kane][Mitiguy1996] and
# [Mitiguy and Reckdahl][Mitiguy2001] but found them to not simplify the
# resulting equations of motion significantly, likely because of the velocity
# constraints complicate the equations of motion.  Additionally, when a bicycle
# is undergoing a steady turn, the three body fixed components of the rear
# frame angular velocity are non-zero, which is in contrast to the three
# generalized coordinate time derivatives need to describe the same angular
# velocity vector, of which two are zero in a steady turn (roll and pitch rate
# are zero, only the yaw rate is non-zero).  Choosing the generalized speeds
# advocated by Mitiguy would imply a that all six speeds would be non-zero in a
# steady turn, in contrast to having only three non-zero speeds when the
# traditional set of speeds is selected.  In particular, this has implications
# for the linearized equations of motion, because the Jacobian matrices can be
# evaluated more simply.
#
# The remaining 6 generalized speeds are associated with the velocity of the
# tire contact points.  The symbolic derivation makes the assumption that 4 of
# the 6 are actual degrees of freedom, while 2 are auxilliary generalized
# speeds used to obtain the normal forces acting at each tire contact point.
# If no-slip rolling investigations are to be performed, it is possible to
# manipulate the equations of motion for the model with 7 degrees of freedom
# such that the model only has 3 degrees of freedom and the constraint forces
# can be determined.  This will be detailed later in the numerical section.
# \begin{aligned}
#  u_6 & \triangleq \mathbf{v}^{RN} \cdot \mathbf{\hat{Y}}_x \\\\
#  u_7 &\triangleq \mathbf{v}^{RN} \cdot \mathbf{\hat{Y}}_y \\\\
#  u_8 &\triangleq \mathbf{v}^{RN} \cdot \mathbf{\hat{Y}}_z  \\\\
#  u_9 & \triangleq \mathbf{v}^{RN} \cdot \mathbf{\hat{n}}_z \\\\
# u_10 & \triangleq \mathbf{v}^{FN} \cdot \mathbf{\hat{n}}_x \\\\
# u_11 &\triangleq \mathbf{v}^{FN} \cdot \mathbf{\hat{n}}_z
# \end{aligned}
#
# [Mitiguy1996]: http://dx.doi.org/10.1177/027836499601500507
# [Mitiguy2001]: https://www1.aiaa.org/content.cfm?pageid=318&volume=24&issue=6&pubid=23&paperid=4849

u = dynamicsymbols('u:12')
ud = [ui.diff(t) for ui in u]
ud_zero_dict = {udi:0 for udi in ud}
qd_u_dict = {qdi : ui for qdi, ui in zip(qd[:6], u[:6])}

code.set_states(q+u, 'x_')

# ### Applied forces and torques
#
# The model includes joint torques at each revolute joint: rear wheel torque,
# steer torque, front wheel torque.  External forces and torques are applied to
# the front and rear assembly; three forces acting at the respective mass
# centers and three torques acting on each frame.  The contact point of each
# tire also has three forces acting on it; these force measure numbers are with
# respect to a reference frame with with one unit vector normal to the ground
# plane, another perpendicular to the line defined by the intersection of the
# wheel plane and ground plane, and the third perpendicular to the other two.

# Internal joint torques (control inputs)
T_rw, T_s, T_fw = symbols('T_rw T_s T_fw')
# Externally applied torques on rear and front frames
T_Rx, T_Ry, T_Rz, T_Fx, T_Fy, T_Fz = symbols('T_Rx T_Ry T_Rz T_Fx T_Fy T_Fz')
# Externally applied forces to rear and front mass centers
F_Rx, F_Ry, F_Rz, F_Fx, F_Fy, F_Fz = symbols('F_Rx F_Ry F_Rz F_Fx F_Fy F_Fz')
# Tire contact forces
G_Rx, G_Ry, G_Rz, G_Fx, G_Fy, G_Fz = symbols('G_Rx G_Ry G_Rz G_Fx G_Fy G_Fz')
# Control/disturbance input vector, 15 x 1
r = Matrix([T_rw, T_s, T_fw,                     # Internal joint torques
            G_Rx, G_Ry, G_Rz,                    # Rear tire contact forces
            F_Rx, F_Ry, F_Rz,                    # Rear frame external forces
            T_Rx, T_Ry, T_Rz,                    # Rear frame external torques
            G_Fx, G_Fy, G_Fz,                    # Front tire contact forces
            F_Fx, F_Fy, F_Fz,                    # Front frame external forces
            T_Fx, T_Fy, T_Fz])                   # Front frame external torques


# ### Reference frames
#
# I declare 5 reference frames:
print("Defining orientations...")
N = ReferenceFrame('N')                    # Inertial frame
Y = N.orientnew('Y', 'Axis', [q[0], N.z])  # Rear yaw frame (heading)
L = Y.orientnew('L', 'Axis', [q[1], Y.x])  # Rear lean frame (roll)
R = L.orientnew('R', 'Axis', [q[2], L.y])  # Fixed to rear frame
RW = R.orientnew('RW', 'Axis', [q[4], R.y])# Fixed to rear wheel
F = R.orientnew('F', 'Axis', [q[3], R.z])  # Fixed to front frame (fork)
FW = F.orientnew('FW', 'Axis', [q[5], F.y])# Fixed to front wheel

# Angular velocity of frames
print("Defining angular velocities...")
Y.set_ang_vel(N, u[0]*Y.z)
L.set_ang_vel(Y, u[1]*Y.x)
R.set_ang_vel(L, u[2]*L.y)
RW.set_ang_vel(R, u[4]*R.y)
F.set_ang_vel(R, u[3]*R.z)
FW.set_ang_vel(F, u[5]*F.y)

print("Defining rear wheel yaw frame unit vectors...")
wc_tc_uv_r = (Y.z - (R.y & Y.z)*R.y).normalize() # Wheel center to tire center
wc_tc_uv_r.simplify()
wyf_x_r = (R.y ^ Y.z).normalize()   # Wheel yaw frame, x, rear
wyf_y_r = Y.z ^ wyf_x_r             # Wheel yaw frame, y, rear

print("Defining front wheel yaw frame unit vectors...")
wc_tc_uv_f = (Y.z.express(F) - (F.y & Y.z)*F.y).normalize() # Wheel center to tire center
wc_tc_uv_f.simplify()
wyf_x_f = (F.y ^ Y.z).normalize()   # Wheel yaw frame, x, front
wyf_y_f = Y.z ^ wyf_x_f             # Wheel yaw frame, y, front


print("Defining positions of rear assembly points...")
gc_r = Point('gc_r')                               # Ground contact, rear
wc_r = gc_r.locatenew('wc_r',                      # Wheel center, rear
                      -rear.r*Y.z - rear.R*wc_tc_uv_r)
mc_r = wc_r.locatenew('mc_r',                      # Mass center, rear
                      rear.a*R.x + rear.b*R.z)
sa_r = wc_r.locatenew('sa_r', rear.c*R.x)          # Steer axis, rear

print("Defining positions of front assembly points...")
gc_f = Point('gc_f')                               # Ground contact, front
wc_f = gc_f.locatenew('wc_f',                      # Wheel center, front
                      -front.r*Y.z - front.R*wc_tc_uv_r)
mc_f = wc_f.locatenew('mc_f',                      # Mass center, front
                      front.a*F.x + front.b*F.z)
sa_f = wc_f.locatenew('sa_f', front.c*F.x)         # Steer axis, front

print("Forming configuration constraint...")
f_c = dot(Y.z, sa_r.pos_from(gc_r) + ls*R.z + gc_f.pos_from(sa_f))
code.generate(Matrix([f_c]), "Bicycle::f_c")

print("Forming velocities of rear assembly points...")
gc_r.set_vel(N, u[6]*wyf_x_r + u[7]*wyf_y_r + u[8]*Y.z)
wc_r.v2pt_theory(gc_r, N, RW)
mc_r.v2pt_theory(wc_r, N, R)
sa_r.v2pt_theory(wc_r, N, R)

print("Forming velocities of front assembly points...")
gc_f.set_vel(N, u[9]*wyf_x_f + u[10]*wyf_y_f + u[11]*Y.z)
wc_f.v2pt_theory(gc_f, N, FW)
mc_f.v2pt_theory(wc_f, N, F)
sa_f.v2pt_theory(wc_f, N, F)

print("Forming velocity constraint and coefficient matrices...")
f_v_vec = sa_r.vel(N) - sa_f.vel(N)
f_v = zeros(3,1)
B = zeros(3, len(u))
B_dq = [zeros(3, len(u)) for i in range(3)]        # List of partials w.r.t. 3 coordinates
B_hess = [[zeros(3,3) for j in range(len(u))] for i in range(3)]  # Hessian w.r.t 3 coords of each entry

for i, uv in enumerate(F):
    f_v[i, 0] = f_v_vec & uv
    for j, uj in enumerate(u):
        B[i, j] = f_v[i, 0].diff(uj)
        for k, qk in enumerate(q[1:3]):
            B_dq[k][i, j] = B[i, j].diff(qk)
            for l, ql in enumerate(q[1:3]):
                B_hess[i][j][k, l] = B_dq[k][i, j].diff(ql)

print("Generating constrint coefficient matrix code...")
code.generate(B, "Bicycle::f_v_coefficient")
#print("Generating constrint coefficient Jacobian matrix code...")
#code.generate(B_dq, "Bicycle::f_v_coefficient_dq")
#print("Generating constrint coefficient Hessian matrix code...")
#code.generate(B_hess, "Bicycle::f_v_coefficient_hessian")
#stop

print("Forming rear assembly partial angular velocities...")
R_N_pav, RW_N_pav, RW_R_pav = partial_velocity([R.ang_vel_in(N),
                                                RW.ang_vel_in(N),
                                                RW.ang_vel_in(R)], u, N)

print("Forming front assembly partial angular velocities...")
F_N_pav, FW_N_pav, FW_F_pav = partial_velocity([F.ang_vel_in(N),
                                                FW.ang_vel_in(N),
                                                FW.ang_vel_in(F)], u, N)

print("Forming rear assembly partial velocities...")
gc_r_pv, mc_r_pv = partial_velocity([gc_r.vel(N), mc_r.vel(N)], u, N)

print("Forming front assembly partial velocities...")
gc_f_pv, mc_f_pv = partial_velocity([gc_f.vel(N), mc_f.vel(N)], u, N)

print("Forming accelerations of rear assembly points...")
wc_r.set_acc(N, wc_r.vel(N).dt(N).subs(qd_u_dict))
mc_r.a2pt_theory(wc_r, N, R)
sa_r.a2pt_theory(wc_r, N, R)

print("Forming accelerations of front assembly points...")
wc_f.set_acc(N, wc_f.vel(N).dt(N).subs(qd_u_dict))
mc_f.a2pt_theory(wc_f, N, F)
sa_f.a2pt_theory(wc_f, N, F)

print("Defining applied forces and torques...")
F_gc_r = G_Rx*wyf_x_r + G_Ry*wyf_y_r + G_Rz*Y.z
F_gc_f = G_Fx*wyf_x_f + G_Fy*wyf_y_f + G_Fz*Y.z
F_mc_r = F_Rx*R.x + F_Ry*R.y + F_Rz*R.z + rear.m*g*Y.z
F_mc_f = F_Fx*F.x + F_Fy*F.y + F_Fz*F.z + front.m*g*Y.z
T_R = T_Rx*R.x + T_Ry*R.y + T_Rz*R.z - T_rw*R.y - T_s*R.z
T_F = T_Fx*F.x + T_Fy*F.y + T_Fz*F.z - T_fw*F.y + T_s*F.z
T_RW = T_rw*R.y
T_FW = T_fw*F.y

print("Defining inertia dyadic of each rigid body...")
I_r = inertia(R, rear.Ixx, rear.Iyy, rear.Izz, 0, 0, rear.Ixz)
I_f = inertia(F, front.Ixx, front.Iyy, front.Izz, 0, 0, front.Ixz)

print("Computing generalized active forces and generalized inertia forces...")
gaf = zeros(len(u), 1)
gif = zeros(len(u), 1)
a_r = mc_r.acc(N)
a_f = mc_f.acc(N)
w_r_n = R.ang_vel_in(N)
w_f_n = F.ang_vel_in(N)
w_rw_r = RW.ang_vel_in(R)
w_fw_f = FW.ang_vel_in(F)
alpha_r_n = R.ang_acc_in(N)
alpha_f_n = F.ang_acc_in(N)
alpha_rw_r = RW.ang_acc_in(R)
alpha_fw_f = FW.ang_acc_in(F)
f_r = zeros(len(u), len(r))
M = zeros(len(u), len(u))
f_star_qu = zeros(len(u), 1)
f_star_qu_dq = zeros(len(u), len(q))
f_star_qu_du = zeros(len(u), len(u))
for i in range(len(u)):
    # Generalized active forces
    gaf[i, 0] = ((F_gc_r & gc_r_pv[i])
               + (F_gc_f & gc_f_pv[i])
               + (F_mc_r & mc_r_pv[i])
               + (F_mc_f & mc_f_pv[i])
               + (T_R & R_N_pav[i])
               + (T_F & F_N_pav[i])
               + (T_RW & RW_N_pav[i])
               + (T_FW & FW_F_pav[i]))

    # Input coefficient matrix
    for j, rj in enumerate(r):
        f_r[i, j] = gaf[i, 0].diff(rj)

    # Generalized inertia forces
    gif[i, 0] = - (rear.m*(a_r & mc_r_pv[i])
                 + front.m*(a_f & mc_f_pv[i])
                 + (((I_r & alpha_r_n) + (w_r_n ^ (I_r & w_r_n)) +
                     rear.J*(alpha_rw_r + (w_r_n ^ w_rw_r))) & R_N_pav[i])
                 + (((I_f & alpha_f_n) + (w_f_n ^ (I_f & w_f_n)) +
                     front.J*(alpha_fw_f + (w_f_n ^ w_fw_f))) & F_N_pav[i])
                 + rear.J*((alpha_r_n + alpha_rw_r) & RW_R_pav[i])
                 + front.J*((alpha_f_n + alpha_fw_f) & FW_F_pav[i]))

    # Coriolis and centripel terms of generalized inertia forces
    f_star_qu[i, 0] = gif[i, 0].subs(ud_zero_dict)

    # Mass matrix and partial derivatives w.r.t u
    for j, (uj, udj) in enumerate(zip(u, ud)):
        f_star_qu_du[i, j] = f_star_qu[i, 0].diff(uj)
        M[i, j] = gif[i, 0].diff(udj)

    # Coriolis/centripetal partial derivatives w.r.t q
    for j, qj in enumerate(q):
        f_star_qu_dq[i, j] = f_star_qu[i, 0].diff(qj)

print("Generating mass matrix code...")
code.generate(M, "Bicycle::mass_matrix")
print("Generating input coefficient matrix code...")
code.generate(f_r, "Biycycle::input_coefficient_matrix")
print("Generating f_star_qu code...")
code.generate(f_star_qu, "Biycycle::f_star_qu")
print("Generating f_star_qu_dq code...")
code.generate(f_star_qu_dq, "Biycycle::f_star_qu_dq")
print("Generating f_star_qu_du code...")
code.generate(f_star_qu_du, "Biycycle::f_star_qu_du")

code.output("bicycle_generated.cc")
