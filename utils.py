from firedrake import *
import numpy as np
from petsc4py import PETSc
print = PETSc.Sys.Print
from irksome import GaussLegendre, Dt, MeshConstant, TimeStepper

def i():
    return as_vector([1., 0., 0.])

def j():
    return as_vector([0., 1., 0.])

def k(twoD=False):
    if twoD:
        return as_vector([0., 1.])
    else:
        return as_vector([0., 0., 1.])

def buo_freq():
    return Constant(1.0e-2) ** 2

def Coriolis_param(no_rotation=False):
    '''
    f-plane Coriolis vector f * k_hat with f = 2 * Omega * sin(latitude).
    The Coriolis term in the momentum equation is cross(Coriolis_param(), u)
    (no extra factor of 2).
    '''
    if no_rotation:
        f = Constant(0.0)
    else:
        Omega = Constant(7.292e-5)
        latitude = pi / 3
        f = 2 * Omega * sin(latitude)
    return as_vector([0., 0., f])

def high_dim_mesh_hierarchy(mh, dim=3):
    '''
    Create a embedded mesh with the same structure in a higher dimension given a mesh hierarchy.
    '''
    meshes = []
    for m in mh:
        x, z = SpatialCoordinate(m)
        coord_fs = VectorFunctionSpace(m, "DG", 1, dim=dim)
        new_coord = assemble(interpolate(as_vector([x, 0, z]), coord_fs))
        new_mesh = Mesh(new_coord)
        new_mesh.init_cell_orientations(j())
        meshes.append(new_mesh)
    
    return HierarchyBase(meshes, mh.coarse_to_fine_cells,
                            mh.fine_to_coarse_cells,
                            mh.refinements_per_level, mh.nested)

def extrude_RT(mesh, k=1):
    CG = FiniteElement("CG", interval, k)
    DG = FiniteElement("DG", interval, k-1)
    CG_DG = TensorProductElement(CG, DG)
    RT_horiz = HDivElement(CG_DG)
    DG_CG = TensorProductElement(DG, CG)
    RT_vert = HDivElement(DG_CG)
    RT_e = RT_horiz + RT_vert
    return FunctionSpace(mesh, RT_e)

def W_theta(mesh, k=0):
    horiz = FiniteElement('DG', interval, k-1)
    vertical = FiniteElement('CG', interval, k)
    V_elt = TensorProductElement(horiz, vertical)
    return FunctionSpace(mesh, V_elt)

def SLB_velocity(u, p, b, w, dt, twoD=False):
    return (
        inner(w, u) * dx
        # + dt ** 2 * buo_freq() / 4 * inner(outer(k(twoD=twoD), k(twoD=twoD)) * u, w) * dx
        - dt * div(w) * p * dx
        - dt / 2 * inner(w, k(twoD=twoD)) * b * dx
        )

def SLB_buoyancy(u, b, q, dt, twoD=False):
    return (
        q * b * dx
        + dt * buo_freq() / 2 * q * inner(u, k(twoD=twoD)) * dx
        )

def SLB_pressure(u, phi):
    return (
        phi * div(u) * dx
        )


def LB_velocity_Irk(u, w, b, p, twoD=False):
    return (
            inner(w, Dt(u)) * dx
            + inner(w, cross(Coriolis_param(), u)) * dx
            - div(w) * p * dx
            - inner(w, k(twoD=twoD)) * b * dx
        )

def LB_buoyancy_Irk(b, q, u, twoD=False):
    return (
            q * Dt(b) * dx
            + buo_freq() * q * inner(k(twoD=twoD), u) * dx
        )

def LB_pressure_Irk(u, phi):
    return (
            phi * div(u) * dx
        )

def LB_velocity(unp1, un, unph, w, bnph, pnp1, n, dt, use_rotation=False, twoD=False, U_mean=0.0):
    eqn = inner(w, (unp1 - un)) * dx 
    eqn -= dt * div(w) * pnp1 * dx
    eqn -= dt * inner(w, k(twoD=twoD)) * bnph * dx
    if use_rotation:
        eqn += dt * inner(w, cross(Coriolis_param(), unph)) * dx

    if U_mean != 0:
        U_vec = as_vector([U_mean, 0]) if twoD else as_vector([U_mean, 0, 0])
        Un = 0.5 * (dot(U_vec, n) + abs(dot(U_vec, n)))
        eqn -= dt * inner(div(outer(U_vec, w)), unph) * dx
        eqn += dt * dot(jump(w), Un('+') * unph('+') - Un('-') * unph('-')) * (dS_v + dS_h)
    return eqn

def LB_buoyancy(bnp1, bn, q, unph, bnph, n, dt, twoD=False, U_mean=0.0):
    eqn = q * (bnp1 - bn) * dx
    eqn += dt * buo_freq() * q * inner(k(twoD=twoD), unph) * dx

    if U_mean != 0:
        U_vec = as_vector([U_mean, 0]) if twoD else as_vector([U_mean, 0, 0])
        Un = 0.5 * (dot(U_vec, n) + abs(dot(U_vec, n)))
        eqn -= dt * div(q * U_vec) * bnph * dx
        eqn += dt * jump(q) * (Un('+') * bnph('+') - Un('-') * bnph('-')) * (dS_v + dS_h)
    return eqn

def LB_pressure(unp1, phi):
    return (
            phi * div(unp1) * dx
        )


def unn_tool(unph, n):
    '''
    Given the Facet normal and a vector, return the form that needed for DG forms.
    '''
    return (
        0.5 * (dot(unph, n) + abs(dot(unph, n)))
    )

def Nonlinear_velocity(unp1, un, unph, w, bnph, pnp1, dt, n, use_rotation=False, twoD=False):
    unn = unn_tool(unph, n)
    eqn = inner(w, (unp1 - un)) * dx
    if use_rotation:
        eqn += dt * inner(w, cross(Coriolis_param(), unph)) * dx
    eqn -= dt * div(w) * pnp1 * dx
    eqn -= dt * inner(w, k(twoD=twoD)) * bnph * dx
    # Advective terms:
    eqn -= dt * inner(div(outer(unph, w)), unph) * dx
    eqn += dt * dot(jump(w), unn('+') * unph('+') - unn('-') * unph('-')) * (dS_v + dS_h)
    return eqn

def Nonlinear_buoyancy(bnp1, bn, bnph, q, unph, dt, n, twoD=False):
    unn = unn_tool(unph, n)
    eqn = q * (bnp1 - bn) * dx
    eqn -= dt * div(q * unph) * bnph * dx
    eqn += dt * jump(q) * (unn('+') * bnph('+') - unn('-') * bnph('-')) * (dS_v + dS_h)
    return eqn

def Nonlinear_pressure(unp1, phi):
    return (
            phi * div(unp1) * dx
        )



def Nonlinear_velocity_Irk(u, w, b, p, n, use_rotation=False, twoD=False):
    unn = unn_tool(u, n)
    eqn = inner(w, Dt(u)) * dx
    if use_rotation:
        eqn += inner(w, cross(Coriolis_param(), u)) * dx
    eqn -= div(w) * p * dx
    eqn -= inner(w, k(twoD=twoD)) * b * dx
    # Advective terms:
    eqn -= inner(div(outer(u, w)), u) * dx
    eqn += dot(jump(w), unn('+') * u('+') - unn('-') * u('-')) * (dS_v + dS_h)
    return eqn

def Nonlinear_buoyancy_Irk(b, q, u, n, twoD=False):
    unn = unn_tool(u, n)
    eqn = q * Dt(b) * dx
    eqn -= div(q * u) * b * dx
    eqn += jump(q) * (unn('+') * b('+') - unn('-') * b('-')) * (dS_v + dS_h)
    # Add linear stratification:
    eqn += buo_freq() * q * inner(k(twoD=twoD), u) * dx
    return eqn

def Nonlinear_pressure_Irk(u, phi):
    return (
            phi * div(u) * dx
        )

