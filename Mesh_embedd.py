from firedrake import *
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print

distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 2)}
old_m = PeriodicIntervalMesh(10, 1.0,distribution_parameters=distribution_parameters)
x, = SpatialCoordinate(old_m)
coord_fs = VectorFunctionSpace(old_m, "DG", 1, dim=2)
new_coord = assemble(interpolate(as_vector([x, 0.]), coord_fs))
m = Mesh(new_coord)
m.init_cell_orientations(as_vector([0.,1.]))
mh = MeshHierarchy(m, refinement_levels=2)
hierarchy = ExtrudedMeshHierarchy(mh, 1.0, layers=[20] * 3, extrusion_type='uniform')
mesh = hierarchy[-1]
finest_mesh_name = "finest"
mesh.name = finest_mesh_name

with CheckpointFile('sol_mesh.h5', 'w') as chk:
    chk.save_mesh(mesh)

x,y, = SpatialCoordinate(mesh)
print(mesh)