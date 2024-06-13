import numpy as np

#%% Generate Quad mesh
### Method 1: simply generate a grid mesh from a rectangle
def generate_quad_mesh(vertices, subdivisions):
    """ Generate a finer grid mesh for a quadrilateral defined by 'vertices'.
    
    Parameters:
    vertices (np.array): An array of shape (4,2), representing the four corners of the quad.
    subdivisions (int): The number of divisions along each edge.
    
    Returns:
    np.array: The generated mesh vertices.
    np.array: The faces of the mesh.
    """
    # Generate points along the edges
    edge_01 = np.linspace(vertices[0], vertices[1], subdivisions + 1)
    edge_32 = np.linspace(vertices[3], vertices[2], subdivisions + 1)
    
    vertices_list = []
    faces = []
    
    # Generate vertices
    for i in range(subdivisions + 1):
        start = edge_01[i]
        end = edge_32[i]
        verts_between = np.linspace(start, end, subdivisions + 1)
        vertices_list.append(verts_between)
    
    vertices_list = np.vstack(vertices_list)
    
    # Generate faces
    for i in range(subdivisions):
        for j in range(subdivisions):
            a = i * (subdivisions + 1) + j
            b = a + 1
            c = (i + 1) * (subdivisions + 1) + j + 1
            d = c - 1
            faces.append([a + 1, b + 1, c + 1, d + 1])  # OBJ indexing starts at 1
    
    return vertices_list, np.array(faces)

#%% Quad to Tri
### Method1: simply connect diagonal
# def quad_to_triangle_mesh(vertices, quad_faces):
#     """
#     Convert a mesh with quadrilateral faces to a mesh with triangular faces.

#     Parameters:
#     vertices : np.array
#         Array of vertices, where each row represents a vertex.
#     quad_faces : np.array
#         Array of quad faces, where each row contains four indices into the vertices array.

#     Returns:
#     triangle_faces : np.array
#         Array of triangular faces, where each row contains three indices into the vertices array.
#     """
#     # Initialize an empty list to store the triangle faces
#     triangle_faces = []

#     # Iterate over each quadrilateral face
#     for face in quad_faces:
#         # Ensure the face has exactly four vertices
#         if len(face) != 4:
#             raise ValueError("All faces in the input must be quadrilaterals with exactly four vertices.")

#         # Split the quadrilateral into two triangles
#         # Triangle 1: vertices 0, 1, 2
#         triangle_faces.append([face[0], face[1], face[2]])
#         # Triangle 2: vertices 0, 2, 3
#         triangle_faces.append([face[0], face[2], face[3]])

#     # Convert the list of triangle faces back to a numpy array for consistency
#     return np.array(triangle_faces)

def quad_to_triangle_mesh(vertices, quad_faces):
    """
    Convert a mesh with quadrilateral faces to a mesh with triangular faces,
    and return both the triangular faces and the unique edges.

    Parameters:
    vertices : np.array
        Array of vertices, where each row represents a vertex.
    quad_faces : np.array
        Array of quad faces, where each row contains four indices into the vertices array.

    Returns:
    triangle_faces : np.array
        Array of triangular faces, where each row contains three indices into the vertices array.
    edges : np.array
        Array of unique edges, where each row contains two indices.
    """
    triangle_faces = []
    edges = set()  # Use a set to prevent duplicate edges

    for face in quad_faces:
        # Ensure the face has exactly four vertices
        if len(face) != 4:
            raise ValueError("All faces in the input must be quadrilaterals with exactly four vertices.")

        # Define two triangles for each quad
        triangles = [
            [face[0], face[1], face[2]],
            [face[0], face[2], face[3]]
        ]
        triangle_faces.extend(triangles)

        # Add edges for these triangles, ensure each edge is added once
        # Use tuple(sorted()) to ensure that each edge is stored in a consistent order
        for tri in triangles:
            edges.update([
                tuple(sorted([tri[0], tri[1]])),
                tuple(sorted([tri[1], tri[2]])),
                tuple(sorted([tri[2], tri[0]]))
            ])

    # Convert the list of triangle faces back to a numpy array for consistency
    ordered_triangle_faces = np.array(triangle_faces, dtype=np.int64)
    ordered_edges = np.array(list(edges), dtype=np.int64)  # Convert set to a sorted numpy array

    return ordered_edges, ordered_triangle_faces