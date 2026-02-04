from firedrake import *
import numpy as np
from petsc4py import PETSc
print = PETSc.Sys.Print

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
    if no_rotation:
        Omega = Constant(0.0)
    else:
        Omega = Constant(7.292e-5)
    theta = pi /3 # Latitude
    return as_vector([0, Omega * sin(theta), Omega * cos(theta)])

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

def extrude_RT(mesh, k=0):
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


def LB_velocity(unp1, un, unph, w, bnph, pnp1, dt, use_rotation=False, twoD=False):
    if use_rotation:
        return (
                inner(w, (unp1 - un)) * dx 
                + dt * inner(w, 2 * cross(Coriolis_param(), unph)) * dx
                - dt * div(w) * pnp1 * dx
                - dt * inner(w, k(twoD=twoD)) * bnph * dx
            )
    else:
        return (
                inner(w, (unp1 - un)) * dx 
                - dt * div(w) * pnp1 * dx
                - dt * inner(w, k(twoD=twoD)) * bnph * dx
            )

def LB_buoyancy(bnp1, bn, q, unph, dt, twoD=False):
    return (
            q * (bnp1 - bn) * dx
            + dt * buo_freq() * q * inner(k(twoD=twoD), unph) * dx
        )

def LB_pressure(unp1, phi):
    return (
            phi * div(unp1) * dx
        )

