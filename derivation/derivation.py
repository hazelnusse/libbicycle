## \class Derivation
#
# This script derives symbolic equations of motion, along with a number of
# related quantities and functions that are useful for the numerical analysis
# of bicycles.  I refer to it as "extended" because it extends the Whipple
# bicycle model in four ways:  wheels with elliptical cross sections instead of
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

from sympy import symbols, zeros, pi, S, Matrix, solve
from sympy.physics.mechanics import *
import numpy as np
Vector.simp = False
from wheelassemblygyrostat import (WheelAssemblyGyrostatParams,
        WheelAssemblyGyrostatInputs)
from utility import NumpyArrayOutput
from benchmark_parameters import (gyrostat_benchmark_parameters,
    gyrostat_benchmark_reference_configuration,
    gyrostat_benchmark_external_forces)

## Generate a OpenGL modelview matrix as a column major length 16 numpy array
#
# @param[in] view_frame A ReferenceFrame which will be the reference frame of
#            the view space
# @param[in] view_origin A Point to be used as the origin of the view space
# @param[in] object_frame A ReferenceFrame or a list of three basis vectors
#            to be used as basis vectors for the object space
# @param[in] object_origin A Point that will act as the origin of the model
#            space
# @returns A NumPy array of shape (16,) which represents an OpenGL modelview
#          matrix (in column major storage layout) mapping object coordinates
#          to view coordinates
def opengl_transformation_matrix(view_frame, view_origin,
                                 object_frame, object_origin):
    if isinstance(object_frame, ReferenceFrame):
        dcm = view_frame.dcm(object_frame)
    else:       # User supplied a list of three basic vectors
        dcm = zeros((3,3))
        for i, uv_cf in enumerate(view_frame):
            for j in range(3):
                dcm[i, j] = uv_cf & object_frame[j]

    r = [object_origin.pos_from(view_origin) & uv for uv in view_frame]
    m = np.zeros((16,), dtype=object)
    m[0] = S(dcm[0, 0])
    m[1] = S(dcm[1, 0])
    m[2] = S(dcm[2, 0])
    m[4] = S(dcm[0, 1])
    m[5] = S(dcm[1, 1])
    m[6] = S(dcm[2, 1])
    m[8] = S(dcm[0, 2])
    m[9] = S(dcm[1, 2])
    m[10] = S(dcm[2, 2])
    m[12] = S(r[0])
    m[13] = S(r[1])
    m[14] = S(r[2])
    m[15] = S(1)
    return m

## Model derivation and C++ code generation
#
# This function forms kinematic and dynamic quantities in symbolic form.  It
# then arranges them in numpy arrays of various shapes, and uses the
# NumpyArrayOutput class to generate output code.
def derivation():
    ## Model description
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
    # * $Ixx, Iyy, Izz, Ixy, Ixz$:  Inertia scalars of combined wheel and frame.
    # * $J$: spin inertia of wheel.
    # * $m$: combined mass of wheel and frame.
    # * $R, r: Wheel major radii, tire cross section radii
    # * $a, b$: distances from wheel center to mass center in x and z directions
    # * $c$: distances from wheel center to steer axis point in x direction
    rear_params = WheelAssemblyGyrostatParams('rear')
    front_params = WheelAssemblyGyrostatParams('front')
    rear_inputs = WheelAssemblyGyrostatInputs('rear')
    front_inputs = WheelAssemblyGyrostatInputs('front')

    # Gravitational constant, time
    ls, g, gt, t = symbols('ls_ g_ gt_ t')

    ## Generalized coordinates
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
    q_min = [q[1], q[2], q[3]]  # lean, pitch, steer
    qd = [qi.diff(t) for qi in q]

    ## Generalized speeds
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
    steady = {u[1]:0, u[2]:0, u[3]:0}
    no_slip = {ui:0 for ui in u[6:]}
    steady_no_slip = steady.copy()
    steady_no_slip.update(no_slip)

    # Generalized inertia forces have centripetal and coriolis portions which
    # are useful to be able to identify.  Each term in general takes the form
    # of f(q) * u_i * u_j.
    cross_terms = []
    unique = set([])
    for ui in u:
        for uj in u:
            if (ui*uj not in unique):
                unique.add(ui * uj)
                cross_terms.append((ui, uj))

    ## Applied forces and torques
    #
    # The model includes joint torques at each revolute joint: rear wheel torque,
    # steer torque, front wheel torque.  External forces and torques are applied to
    # the front and rear assembly; three forces acting at the respective mass
    # centers and three torques acting on each frame.  The contact point of each
    # tire also has three forces acting on it; these force measure numbers are with
    # respect to a reference frame with with one unit vector normal to the ground
    # plane, another perpendicular to the line defined by the intersection of the
    # wheel plane and ground plane, and the third perpendicular to the other two.
    T_s = symbols('steer_torque_')
    # Control/disturbance input vector, 22 x 1
    r = np.array(
        [rear_inputs.Tw,                       # Rear wheel torque
         rear_inputs.Tx, rear_inputs.Ty, rear_inputs.Tz,     # Torques on rear frame
         rear_inputs.Gx, rear_inputs.Gy, rear_inputs.Gz,     # Rear wheel ground contact forces
         rear_inputs.Fx, rear_inputs.Fy, rear_inputs.Fz,     # Forces at rear mass center
         front_inputs.Tw,                      # Front wheel torque
         front_inputs.Tx, front_inputs.Ty, front_inputs.Tz,  # Torques on front frame
         front_inputs.Gx, front_inputs.Gy, front_inputs.Gz,  # Front wheel ground contact forces
         front_inputs.Fx, front_inputs.Fy, front_inputs.Fz,  # Forces at front mass center
         T_s, g])                       # Steer torque, gravity

    # Reference frames
    print("Defining orientations...")
    N = ReferenceFrame('N')                    # Inertial frame
    H = N.orientnew('H', 'Axis', [-gt, N.x])   # Horizontal frame
    Y = N.orientnew('Y', 'Axis', [q[0], N.z])  # Rear yaw frame
    L = Y.orientnew('L', 'Axis', [q[1], Y.x])  # Rear roll frame
    R = L.orientnew('R', 'Axis', [q[2], L.y])  # Rear assembly fixed frame
    RW = R.orientnew('RW', 'Axis', [q[4], R.y])# Rear wheel fixed frame
    F = R.orientnew('F', 'Axis', [q[3], R.z])  # Front assembly fixed frame
    FW = F.orientnew('FW', 'Axis', [q[5], F.y])# Front wheel fixed frame
    # Camera related references frames
    cam_angles = symbols('azimuth elevation twist')
    cam_position = symbols('cam_x cam_y cam_z')
    # Intermediate frames
    C__ = Y.orientnew('C__', 'Axis', [-pi/S(2), Y.y])
    C_ = C__.orientnew('C_', 'Axis', [pi/S(2), C__.z])
    # When camera angles are all zero, camera should be pointed down Y.x axis
    C_az = C_.orientnew('C_az', 'Axis', [cam_angles[0], C_.y])      # Azimuth
    C_el = C_az.orientnew('C_el', 'Axis', [cam_angles[1], -C_az.x]) # Elevation
    C = C_el.orientnew('C', 'Axis', [cam_angles[2], -C_el.z])  # Twist
    # OpenGL camera has x right, y up, and is pointed down -z axis, the
    # following rotatations make this happen

    # Angular velocity of frames
    print("Defining angular velocities and angular accelerations...")
    # Rear frame
    R.set_ang_vel(N, (u[0]*Y.z + u[1]*L.x + u[2]*L.y).express(R))
    R.set_ang_acc(N, R.ang_vel_in(N).diff(t, R).subs(qd_u_dict))

    # Front frame
    F.set_ang_vel(N, R.ang_vel_in(N).express(F) + u[3]*F.z)
    F.set_ang_acc(N, F.ang_vel_in(N).diff(t, F).subs(qd_u_dict))

    # Rear wheel
    RW.set_ang_vel(N, R.ang_vel_in(N) + u[4]*R.y)
    RW.set_ang_vel(R, u[4]*R.y)
    RW.set_ang_acc(N, RW.ang_vel_in(N).diff(t, R).subs(qd_u_dict) +
                      (R.ang_vel_in(N) ^ RW.ang_vel_in(N)))
    # Front wheel
    FW.set_ang_vel(N, F.ang_vel_in(N) + u[5]*F.y)
    FW.set_ang_vel(F, u[5]*F.y)
    FW.set_ang_acc(N, FW.ang_vel_in(N).diff(t, F).subs(qd_u_dict)
                    + cross(F.ang_vel_in(N), FW.ang_vel_in(N)))
    FW.set_ang_acc(F, ud[5]*F.y)

    print("Defining rear wheel yaw frame unit vectors...")
    # Ground normal projected onto rear wheel plane
    Yz_R = Y.z.express(R)
    wc_tc_uv_r = (Yz_R - (R.y & Yz_R)*R.y)
    wc_tc_uv_r_mag = wc_tc_uv_r.magnitude()
    wc_tc_uv_r /= wc_tc_uv_r_mag
    wyf_x_r = (R.y ^ Yz_R).normalize()   # Wheel yaw frame, x, rear
    wyf_y_r = Yz_R ^ wyf_x_r             # Wheel yaw frame, y, rear

    print("Defining front wheel yaw frame unit vectors...")
    # Ground normal projected onto front wheel plane
    Yz_F = Y.z.express(F)
    wc_tc_uv_f = (Yz_F - (F.y & Yz_F)*F.y)
    wc_tc_uv_f_mag = wc_tc_uv_f.magnitude()
    wc_tc_uv_f /= wc_tc_uv_f_mag
    wyf_x_f = (F.y ^ Yz_F).normalize()   # Wheel yaw frame, x, front
    wyf_y_f = Yz_F ^ wyf_x_f             # Wheel yaw frame, y, front

    print("Defining kinematics of rear assembly points...")
    # Rear ground contact
    gc_r = Point('gc_r')
    gc_r.set_vel(N, u[6]*wyf_x_r + u[7]*wyf_y_r + u[8]*Yz_R)

    # Rear wheel center
    wc_r = gc_r.locatenew('wc_r', -(rear_params.r*Yz_R + rear_params.R*wc_tc_uv_r))
    wc_r.v2pt_theory(gc_r, N, RW)
    wc_r.set_acc(N, wc_r.vel(N).diff(t, R).subs(qd_u_dict)
                  + (R.ang_vel_in(N) ^ wc_r.vel(N)).subs(qd_u_dict))

    # Rear mass center
    mc_r = wc_r.locatenew('mc_r', rear_params.a*R.x + rear_params.b*R.y + rear_params.c*R.z)
    mc_r.v2pt_theory(wc_r, N, R)
    mc_r.a2pt_theory(wc_r, N, R)

    # Rear steer axis point
    sa_r = wc_r.locatenew('sa_r', rear_params.d*R.x)
    sa_r.v2pt_theory(wc_r, N, R)

    print("Defining kinematics of front assembly points...")
    # Front ground contact
    gc_f = Point('gc_f')
    gc_f.set_vel(N, u[9]*wyf_x_f + u[10]*wyf_y_f + u[11]*Yz_F)

    # Front wheel center
    wc_f = gc_f.locatenew('wc_f', -(front_params.r*Yz_F + front_params.R*wc_tc_uv_f))
    wc_f.v2pt_theory(gc_f, N, FW)
    wc_f.set_acc(N, wc_f.vel(N).diff(t, F).subs(qd_u_dict)
                  + (F.ang_vel_in(N) ^ wc_f.vel(N)).subs(qd_u_dict))

    # Front mass center
    mc_f = wc_f.locatenew('mc_f', front_params.a*F.x + front_params.b*F.y + front_params.c*F.z)
    mc_f.v2pt_theory(wc_f, N, F)
    mc_f.a2pt_theory(wc_f, N, F)

    # Front steer axis point
    sa_f = wc_f.locatenew('sa_f', front_params.d*F.x)
    sa_f.v2pt_theory(wc_f, N, F)

    # Connect assemblies together so we can form holonomic constraint
    sa_f.set_pos(sa_r, ls*R.z)

    # Generating kinematic positions in Yaw frame for debugging purposes
    wc_r_rel_gc_r = np.array([wc_r.pos_from(gc_r) & uv for uv in Y],
                              dtype=object)
    mc_r_rel_gc_r = np.array([mc_r.pos_from(gc_r) & uv for uv in Y],
                              dtype=object)
    sa_r_rel_gc_r = np.array([sa_r.pos_from(gc_r) & uv for uv in Y],
                              dtype=object)
    mc_f_rel_gc_r = np.array([mc_f.pos_from(gc_r) & uv for uv in Y],
                              dtype=object)
    wc_f_rel_gc_r = np.array([wc_f.pos_from(gc_r) & uv for uv in Y],
                              dtype=object)
    sa_f_rel_gc_r = np.array([sa_f.pos_from(gc_r) & uv for uv in Y],
                              dtype=object)
    # Position from rear contact to front contact
    gc_f_rel_gc_r = np.array([gc_f.pos_from(gc_r) & uv for uv in Y],
                              dtype=object)

    print("Defining postion of camera relative to rear wheel contact...")
    cam_origin = gc_r.locatenew('cam_origin', cam_position[0]*Y.x
                                            + cam_position[1]*Y.y
                                            + cam_position[2]*Y.x)

    print("Forming OpenGL transformation matrices...")
    gc_r_ogl = opengl_transformation_matrix(C, cam_origin,
                                            [wyf_x_r, wyf_y_r, Y.z], gc_r)
    wc_r_ogl = opengl_transformation_matrix(C, cam_origin, RW, wc_r)
    mc_r_ogl = opengl_transformation_matrix(C, cam_origin, R, mc_r)

    gc_f_ogl = opengl_transformation_matrix(C, cam_origin,
                                            [wyf_x_f, wyf_y_f, Y.z], gc_f)
    wc_f_ogl = opengl_transformation_matrix(C, cam_origin, FW, wc_f)
    mc_f_ogl = opengl_transformation_matrix(C, cam_origin, F, mc_f)

    NO = gc_r.locatenew('NO', -q[6]*N.x - q[7]*N.y)    # Inertial origin
    N_ogl    = opengl_transformation_matrix(C, cam_origin, N, NO)

    # Coordinates of rear wheel center for comparison to Basu-Mandal steady
    # turning results
    x_bm, y_bm, z_bm = symbols('x_bm y_bm z_bm')
    x_dot_bm, y_dot_bm, z_dot_bm = symbols('x_dot_bm y_dot_bm z_dot_bm')
    xyz_eqns = [x_bm - (wc_r.pos_from(NO) & N.x), 
                y_bm - (wc_r.pos_from(NO) & N.y),
                z_bm - (wc_r.pos_from(NO) & N.z)]
    soln = solve(xyz_eqns[:2], [q[6], q[7]])
    q6q7_from_bm = np.array([soln[q[6]], soln[q[7]]], dtype=object)
    wc_r_v_N = wc_r.vel(N)
    xyz_dot_bm = np.array([wc_r_v_N & uv for uv in N], dtype=object)

    print("Forming kinematic ODE's...")
    f_1 = np.zeros((8,), dtype=object)
    f_1_dq = np.zeros((8, 8), dtype=object)
    f_1_du = np.zeros((8, 12), dtype=object)
    for i, qdi in enumerate(qd[:6]):
        f_1[i] = -u[i]
        f_1_du[i, i] = S(-1)

    v_gc_rw = (gc_r.pos_from(wc_r).diff(t, R).subs(qd_u_dict) +
            (R.ang_vel_in(RW) ^ gc_r.pos_from(wc_r)))
    f_1[6] = -(v_gc_rw & N.x)
    f_1[7] = -(v_gc_rw & N.y)
    for j, qj in enumerate(q):
        f_1_dq[6, j] = f_1[6].diff(qj)
        f_1_dq[7, j] = f_1[7].diff(qj)

    for j, uj in enumerate(u):
        f_1_du[6, j] = f_1[6].diff(uj)
        f_1_du[7, j] = f_1[7].diff(uj)

    print("Forming configuration constraint...")
    f_c = np.array([Y.z & gc_f.pos_from(gc_r)])
    f_c_dq = np.array([f_c[0].diff(qi) for qi in q])

    print("Forming velocity constraint and partial derivative matrices...")
    f_v_vec = sa_r.vel(N) + (R.ang_vel_in(N) ^ sa_f.pos_from(sa_r)) - sa_f.vel(N)
    f_v = np.zeros((3,), dtype=object)
    f_v_du = np.zeros((3, 12), dtype=object)
    f_v_dudq = np.zeros((3, 12, 3), dtype=object)
    f_v_dudt = np.zeros((3, 12), dtype=object)
    f_v_dudtdq = np.zeros((3, 12, 3), dtype=object)

    for i, uv in enumerate(F):
        f_v[i] = f_v_vec & uv
        for qdi in qd:
            assert(f_v[i].diff(qdi) == 0)
        for j, uj in enumerate(u):
            f_v_du[i, j] = f_v[i].diff(uj)
            f_v_dudt[i, j] = f_v_du[i, j].diff(t).subs(qd_u_dict)
            for k, qk in enumerate(q_min):
                f_v_dudq[i, j, k] = f_v_du[i, j].diff(qk)
                f_v_dudtdq[i, j, k] = f_v_dudt[i, j].diff(qk)

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
    gc_f_pv, mc_f_pv = partial_velocity([gc_f.vel(N), mc_f.vel(N)], u, N)

    print("Defining applied forces and torques...")
    F_gc_r = rear_inputs.Gx*wyf_x_r + rear_inputs.Gy*wyf_y_r + rear_inputs.Gz*Y.z
    F_gc_f = front_inputs.Gx*wyf_x_f + front_inputs.Gy*wyf_y_f + front_inputs.Gz*Y.z
    F_mc_r = rear_inputs.Fx*R.x + rear_inputs.Fy*R.y + rear_inputs.Fz*R.z + rear_inputs.m*g*H.z
    F_mc_f = front_inputs.Fx*F.x + front_inputs.Fy*F.y + front_inputs.Fz*F.z + front_inputs.m*g*H.z
    T_R = rear_inputs.Tx*R.x + rear_inputs.Ty*R.y + rear_inputs.Tz*R.z - rear_inputs.Tw*R.y - T_s*R.z
    T_F = front_inputs.Tx*F.x + front_inputs.Ty*F.y + front_inputs.Tz*F.z - front_inputs.Tw*F.y + T_s*F.z
    T_RW = rear_inputs.Tw*R.y
    T_FW = front_inputs.Tw*F.y

    print("Defining inertia dyadic of each rigid body...")
    I_r = inertia(R, rear_params.Ixx, rear_params.Iyy, rear_params.Izz,
                     rear_params.Ixy, rear_params.Iyz, rear_params.Ixz)
    I_f = inertia(F, front_params.Ixx, front_params.Iyy, front_params.Izz,
                     front_params.Ixy, rear_params.Iyz, front_params.Ixz)

    # Form kinetic energy of rear and front assemblies
    print("Forming kinetic energy expressions")
    ke_rear_t = rear_params.m / 2.0 * (mc_r.vel(N).magnitude() ** 2)
    ke_rear_a = (((R.ang_vel_in(N) & (I_r & R.ang_vel_in(N)))
                  + rear_params.J * (RW.ang_vel_in(R).magnitude() ** 2)) / 2.0
                 + rear_params.J * (R.ang_vel_in(N) & RW.ang_vel_in(R)))

    ke_front_t = front_params.m / 2.0 * (mc_f.vel(N).magnitude() ** 2)
    ke_front_a = (((F.ang_vel_in(N) & (I_f & F.ang_vel_in(N)))
                   + front_params.J * (FW.ang_vel_in(F).magnitude() ** 2)) / 2.0
                  + front_params.J * (F.ang_vel_in(N) & FW.ang_vel_in(F)))

    pe_rear = -rear_params.m * g * (mc_r.pos_from(gc_r) & Y.z)
    pe_front = -front_params.m * g * (mc_f.pos_from(gc_f) & Y.z)

    ke_pe = np.array([ke_rear_t, ke_rear_a, ke_front_t, ke_front_a, pe_rear,
        pe_front])

    print("Computing rear and front wheel turn radii...")
    # gc_r_to_gc_f = R_path_r * wyf_y_r - R_path_f * wyf_y_f
    A = zeros((2, 2))
    A[0, 0] = (wyf_y_r & Y.x)
    A[0, 1] = (-wyf_y_f & Y.x)
    A[1, 0] = (wyf_y_r & Y.y)
    A[1, 1] = (-wyf_y_f & Y.y)
    rhs = Matrix([gc_f.pos_from(gc_r) & Y.x, gc_f.pos_from(gc_r) & Y.y])
    path_radii = np.array((A.inverse_ADJ() * rhs).transpose().tolist()[0],
                          dtype=object)

    print("Computing generalized active forces and generalized inertia forces...")
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
    # Partial derivatives w.r.t. are only taken with respect to q[1], q[2],
    # q[3], the lean, pitch and steer angles, since the others do not appear in
    # any of the dynamics
    gaf = np.zeros((len(u),), dtype=object)
    gaf_dq = np.zeros((len(u), len(q_min)), dtype=object)
    gaf_dr = np.zeros((len(u), len(r)), dtype=object)
    gif = np.zeros((len(u),), dtype=object)
    gif_ud_zero = np.zeros((len(u),), dtype=object)
    gif_ud_zero_steady = np.zeros((len(u),), dtype=object)
    gif_ud_zero_steady_dudu = np.zeros((len(u), len(cross_terms)), dtype=object)
    gif_dud = np.zeros((len(u), len(u)), dtype=object)
    gif_dud_dq = np.zeros((len(u), len(u), len(q_min)), dtype=object)
    gif_ud_zero_dq = np.zeros((len(u), len(q_min)), dtype=object)
    gif_ud_zero_du = np.zeros((len(u), len(u)), dtype=object)

    for i in range(len(u)):
        # Generalized active forces
        gaf[i] = ((F_gc_r & gc_r_pv[i])
                + (F_gc_f & gc_f_pv[i])
                + (F_mc_r & mc_r_pv[i])
                + (F_mc_f & mc_f_pv[i])
                + (T_R & R_N_pav[i])
                + (T_F & F_N_pav[i])
                + (T_RW & RW_N_pav[i])
                + (T_FW & FW_F_pav[i]))

        # Generalized inertia forces
        gif[i] = - (rear_params.m*(a_r & mc_r_pv[i])
                  + front_params.m*(a_f & mc_f_pv[i])
                  + (((I_r & alpha_r_n) + (w_r_n ^ (I_r & w_r_n)) +
                      rear_params.J*(alpha_rw_r + (w_r_n ^ w_rw_r))) & R_N_pav[i])
                  + (((I_f & alpha_f_n) + (w_f_n ^ (I_f & w_f_n)) +
                      front_params.J*(alpha_fw_f + (w_f_n ^ w_fw_f))) & F_N_pav[i])
                  + rear_params.J*((alpha_r_n + alpha_rw_r) & RW_R_pav[i])
                  + front_params.J*((alpha_f_n + alpha_fw_f) & FW_F_pav[i]))

        # Coriolis and centripel terms of generalized inertia forces
        gif_ud_zero[i] = gif[i].subs(ud_zero_dict)
        gif_ud_zero_steady[i] = gif_ud_zero[i].subs(steady)

        for j, (ui, uj) in enumerate(cross_terms):
            gif_ud_zero_steady_dudu[i, j] = gif_ud_zero_steady[i].diff(ui).diff(uj)
            if ui == uj:
                gif_ud_zero_steady_dudu[i, j] /= 2

        # Partial derivatives w.r.t. q
        for j, qj in enumerate(q_min):
            gaf_dq[i, j] = gaf[i].diff(qj)
            gif_ud_zero_dq[i, j] = gif_ud_zero[i].diff(qj)

        # Input coefficient matrix
        for j, rj in enumerate(r):
            gaf_dr[i, j] = gaf[i].diff(rj)

        # Mass matrix and partial derivatives w.r.t u
        for j, (uj, udj) in enumerate(zip(u, ud)):
            gif_ud_zero_du[i, j] = gif_ud_zero[i].diff(uj)
            gif_dud[i, j] = gif[i].diff(udj)
            # Partial derivatives of mass matrix w.r.t. q, needed for
            # linearization when du/dt != 0
            for k, qk in enumerate(q_min):
                gif_dud_dq[i, j, k] = gif_dud[i, j].diff(qk)

    # Identify non-zero columns
    non_zero_columns = []
    for j in range(len(cross_terms)):
        if gif_ud_zero_steady_dudu[:, j].any() != 0:
            non_zero_columns.append(j)

    gif_ud_zero_steady_dudu = gif_ud_zero_steady_dudu[:, non_zero_columns]
    cross_terms = np.array([ui*uj for ui, uj in cross_terms], dtype=object)[non_zero_columns]

    # Output code generation
    code = NumpyArrayOutput(class_name='Bicycle',
                            includes=['<cmath>', '"bicycle.h"'],
                            enclosing_namespace='bicycle',
                            using_declarations=['::std::sin', '::std::cos',
                                                '::std::pow', '::std::sqrt'])
    code.add_regex(r'([_0-9a-zA-Z]+)_(rear_params|front_params)', r'\2_params.\1')
    code.add_regex(r'([_0-9a-zA-Z]+)_(rear_inputs|front_inputs)', r'\2_inputs.\1')
    code.set_states(q+u, 'state_')

    print("Generating OpenGL modelview matrices...")
    code.generate(gc_r_ogl, "gc_r_ogl")
    code.generate(wc_r_ogl, "wc_r_ogl")
    code.generate(mc_r_ogl, "mc_r_ogl")
    code.generate(gc_f_ogl, "gc_f_ogl")
    code.generate(wc_f_ogl, "wc_f_ogl")
    code.generate(mc_f_ogl, "mc_f_ogl")
    code.generate(N_ogl, "N_ogl")

    print("Generating kinematics related outputs...")
    code.generate(wc_r_rel_gc_r, "rear_wheel_center_point")
    code.generate(mc_r_rel_gc_r, "rear_mass_center_point")
    code.generate(sa_r_rel_gc_r, "rear_steer_axis_point")
    code.generate(wc_f_rel_gc_r, "front_wheel_center_point")
    code.generate(mc_f_rel_gc_r, "front_mass_center_point")
    code.generate(sa_f_rel_gc_r, "front_steer_axis_point")
    code.generate(gc_f_rel_gc_r, "front_ground_contact_point")

    print("Generating wheel path radii outputs...")
    code.generate(path_radii, "path_radii")

    print("Generating kinetic and potential energy outputs...")
    code.generate(ke_pe, "ke_pe")

    print("Generating Basu-Mandal rear wheel center coordinates...")
    code.generate(q6q7_from_bm, "q6q7_from_bm")

    print("Generating Basu-Mandal rear wheel center velocity...")
    code.generate(xyz_dot_bm, "xyz_dot_bm")

    print("Generating configuration constraint (f_c) code...")
    code.generate(f_c, "f_c")
    print("Generating configuration constraint partial derivatives code" +
          "(f_c_dq) code...")
    code.generate(f_c_dq, "f_c_dq")

    print("Generating constraint coefficient matrix (f_v_du) code...")
    code.generate(f_v_du, "f_v_du")
    print("Generating constraint coefficient Jacobian matrix " +
          "(f_v_dudq) code...")
    code.generate(f_v_dudq, "f_v_dudq")
    print("Generating constraint coefficient matrix " +
          "time derivative (f_v_dudt) code...")
    code.generate(f_v_dudt, "f_v_dudt")
    print("Generating constraint coefficient Jacobian matrix " +
          "time derivative (f_v_dudtdq) code...")
    code.generate(f_v_dudtdq, "f_v_dudtdq")

    print("Generating kinematic differential equations (f_1, f_1_dq, f_1_du)")
    code.generate(f_1, "f_1")
    code.generate(f_1_dq, "f_1_dq")
    code.generate(f_1_du, "f_1_du")

    print("Generating mass matrix (gif_dud) code...")
    code.generate(gif_dud, "gif_dud")
    print("Generating mass matrix partial derivative (gif_dud_dq) code...")
    code.generate(gif_dud_dq, "gif_dud_dq")
    print("Generating coriolis/centripetal (gif_ud_zero) code...")
    code.generate(gif_ud_zero, "gif_ud_zero")
    print("Generating steady coriolis/centripetal (gif_ud_zero_steady) code...")
    code.generate(gif_ud_zero_steady, "gif_ud_zero_steady")

    print("Generating steady coriolis/centripetal (gif_ud_zero_steady_dudu) code...")
    code.generate(gif_ud_zero_steady_dudu, "gif_ud_zero_steady_dudu")

    print("Generating non-zero steady coriolis/centripetal (gif_ud_zero_steady_cross_terms) code...")
    code.generate(cross_terms, "gif_ud_zero_steady_cross_terms")

    print("Generating partial derivatives of coriolis/centripetal " +
          "(gif_ud_zero_dq, gif_ud_zero_du) code...")
    code.generate(gif_ud_zero_dq, "gif_ud_zero_dq")
    code.generate(gif_ud_zero_du, "gif_ud_zero_du")

    print("Generating generalized active forces (gaf) code...")
    code.generate(gaf, "gaf")
    print("Generating generalized active forces partial derivative matrices " +
          "(gaf_dq, gaf_dr) code...")
    code.generate(gaf_dq, "gaf_dq")
    code.generate(gaf_dr, "gaf_dr")

    code.output("bicycle_generated.cc")

if __name__ == "__main__":
    derivation()
