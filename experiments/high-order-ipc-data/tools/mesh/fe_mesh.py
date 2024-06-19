import numpy as np

import meshio
import igl

from weights.bases import gmsh_to_basis_order, basis_order_to_gmsh
from .attach_higher_order_nodes import attach_higher_order_nodes
from .utils import faces
from .write_obj import write_obj

tetra_type_to_order = {
    "tetra": 1,
    "tetra10": 2,
    "tetra20": 3,
    "tetra35": 4,
}

order_to_tetra_type = {
    1: "tetra",
    2: "tetra10",
    3: "tetra20",
    4: "tetra35",
}


class FEMesh:
    def __init__(self, filename, cell_i=0) -> None:
        mesh = meshio.read(filename)

        self.path = filename
        self.name = filename.stem

        self.V = mesh.points

        #################################
        # self.order = 1  # Added to set Default order
        # self.T = mesh.cells[cell_i].data if len(mesh.cells) > cell_i else mesh.cells[0].data
        # self.initialize_mesh_attributes()
        #################################

        assert(len(mesh.cells) > cell_i)

        for cell in mesh.cells:
            if cell.type == "triangle":
                print('Check cell data')
                self.T = cell.data
                break
            elif "tetra" in cell.type:
                print('check is tetra')
                self.order = tetra_type_to_order[cell.type]
                self.T = cell.data[:, gmsh_to_basis_order[self.order]]
                break
            # else:
            #     raise NotImplementedError(
            #         f"FEMesh not implemented for {cell.type}")

        *_, J = igl.remove_unreferenced(self.V, self.P1())
        self.in_vertex_ids = J
        self.in_faces = self.faces()
        self.in_edges = self.edges()

        # if self.order != 1:
        #     # TODO: Reorder the nodes to be [V; E; F; C] and set to input positions
        #     self.attach_higher_order_nodes(self.order)

    ################################################################
    def initialize_mesh_attributes(self):
        # Calculate boundary edges
        self.BE = igl.boundary_facets(self.T)
        # Calculate full edges
        self.E = igl.edges(self.T)
        # Map boundary edges to full edges
        self.BE2E = self.map_boundary_edges_to_edges()

    def map_boundary_edges_to_edges(self):
        # Create a dictionary to map full edges to their indices
        edge_to_index = {tuple(sorted(e)): i for i, e in enumerate(self.E)}
        # Map each boundary edge to its corresponding index in the full edge list
        BE2E = np.array([edge_to_index[tuple(sorted(e))] for e in self.BE])
        return BE2E
    ################################################################
    
    def dim(self):
        return self.V.shape[1]

    def n_nodes(self):
        return self.V.shape[0]

    def n_vertices(self):
        return self.in_vertex_ids.size

    def n_edges(self):
        return self.E.shape[0]

    def n_faces(self):
        return self.F.shape[0]

    def P1(self):
        return self.T[:, :self.dim()+1]

    def faces(self):
        return faces(self.P1())

    def edges(self):
        return igl.edges(self.P1())

    def boundary_faces(self):
        return igl.boundary_facets(self.P1())

    def attach_higher_order_nodes(self, order):
        """Insert higher order indices at end."""
        # Maintaint the same vertex positions throughout node shuffle
        if order == self.order:
            V_old = self.V.copy()
            T_old = self.T.copy()

        # Remove any currently attached HO nodes
        V, self.T, *_ = igl.remove_unreferenced(self.V, self.P1())

        # Replace the nodes and tets with the HO ones
        self.V, self.T = attach_higher_order_nodes(
            V, self.edges(), self.faces(), self.T, order)
        assert((self.V[:V.shape[0]] == V).all())

        if order == self.order:
            V_new = self.V.copy()
            for i, vi_new in np.ndenumerate(self.T):
                assert(i[1] > 3 or (self.V[vi_new] == V_old[T_old[i]]).all())
                V_new[vi_new] = V_old[T_old[i]]
            self.V = V_new
            # self.save("test.msh")
            # exit()

        self.order = order

    def save(self, filename):
        print(f"saving FEM mesh to {filename}")
        meshio.write(filename, meshio.Mesh(
            self.V, [(order_to_tetra_type[self.order],
                      self.T[:, basis_order_to_gmsh[self.order]])]),
                     file_format="gmsh")
