import numpy as np
import h5py
import scipy.sparse as sp

#%% Load mesh

#%% Save mesh
### Save mesh as .obj format
def save_mesh_to_obj(vertices, faces, filename="mesh.obj"):
    """ Save the mesh vertices and faces to a .obj file, adjusting for 1-based indexing in the .obj format.
    
    Parameters:
    vertices (np.array): The mesh vertices.
    faces (np.array): The mesh faces.
    filename (str): The name of the .obj file to save.
    """
    with open(filename, 'w') as file:
        file.write("# OBJ file\n")
        for vert in vertices:
            file.write(f"v {vert[0]} {vert[1]} 0.0\n")  # Assuming the vertices are 2D, z-coordinate is set to 0.0
        for face in faces:
            # Increment each index in the face by 1 for OBJ 1-based indexing
            face_str = ' '.join(str(index + 1) for index in face)
            file.write(f"f {face_str}\n")

#%% Save mapping as .hdf5 format
def save_mapping_to_hdf5(W, ordered_edges, ordered_faces, filename="mapping.hdf5"):
    """ Save the mapping to a .hdf5 file.
    
    Parameters:
    W: The mapping to save.
    filename (str): The name of the .hdf5 file to save.
    """
    # Convert W to a sparse matrix in COO format
    W_sparse = sp.coo_matrix(W)

    # Save to HDF5
    with h5py.File(filename, 'w') as f:
        # Save edges and faces with explicit data types
        f.create_dataset('ordered_edges', data=ordered_edges, dtype=np.int64)
        f.create_dataset('ordered_faces', data=ordered_faces, dtype=np.int64)
        
        # Save the weight triplets
        weight_group = f.create_group('weight_triplets')
        weight_group.attrs['shape'] = np.array(W_sparse.shape)  # Store the shape of the original matrix W
        weight_group.create_dataset('cols', data=W_sparse.col, dtype=np.int32)  # Column indices from COO
        weight_group.create_dataset('rows', data=W_sparse.row, dtype=np.int32)  # Row indices from COO
        weight_group.create_dataset('values', data=W_sparse.data, dtype=np.float64)  # Non-zero values from COO

    print(f"Data saved successfully in {filename}")