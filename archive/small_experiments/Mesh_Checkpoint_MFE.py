from firedrake import *

m = PeriodicIntervalMesh(10, 1.0)
mh = MeshHierarchy(m, refinement_levels=2)
hierarchy = ExtrudedMeshHierarchy(mh, 1.0, layers=[20] * 3, extrusion_type='uniform')
meshes = []
for m in hierarchy:
    x,z = SpatialCoordinate(m)
    coord_fs = VectorFunctionSpace(m, "DQ", 1, dim=3)
    new_coord = assemble(interpolate(as_vector([x, 0, z]), coord_fs))
    new_mesh = Mesh(new_coord)
    new_mesh.init_cell_orientations(as_vector([0.,1.,0.]))
    meshes.append(new_mesh)
new_mh = HierarchyBase(meshes, hierarchy.coarse_to_fine_cells, hierarchy.fine_to_coarse_cells, hierarchy.refinements_per_level, hierarchy.nested)
mesh = new_mh[-1]

with CheckpointFile('sol_mesh.h5', 'w') as chk:
    chk.save_mesh(mesh)

# todo check the mesh is having the correct DOF.