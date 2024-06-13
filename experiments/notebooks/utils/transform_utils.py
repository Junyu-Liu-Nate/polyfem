import numpy as np

### Compute the transformation matrix from the source points to the target points:
def compute_transformation_matrix(src, dst):
    """Compute the transformation matrix from src points to dst points."""
    A = np.array([
        [src[0][0], src[1][0], src[2][0], src[3][0]],
        [src[0][1], src[1][1], src[2][1], src[3][1]],
        [1, 1, 1, 1],
        [0, 0, 0, 0]
    ])
    B = np.array([
        [dst[0][0], dst[1][0], dst[2][0], dst[3][0]],
        [dst[0][1], dst[1][1], dst[2][1], dst[3][1]],
        [1, 1, 1, 1],
        [0, 0, 0, 0]
    ])

    M = np.dot(B, np.linalg.inv(A))
    return M[:3]

### This function maps a quadrilateral defined by four vertices to a unit square:
def transform_to_square(vertices):
    """Transform a given quadrilateral to a unit square with lower left corner at (0, 0) and length 1."""
    # Define the unit square corners
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    
    # Compute the transformation matrix to the unit square
    M = compute_transformation_matrix(vertices, square)
    
    # Transform vertices
    transformed_vertices = np.dot(M, np.vstack((vertices.T, np.ones(4))))
    return transformed_vertices[:2].T

### This function reverses the transformation, mapping the unit square back to the original quadrilateral:
def transform_back_to_quad(square_vertices, original_vertices):
    """Transform the unit square vertices back to the original quadrilateral vertices."""
    # Define the unit square corners
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    
    # Compute the transformation matrix back to the original quad
    M = compute_transformation_matrix(square, original_vertices)
    
    # Transform vertices
    transformed_vertices = np.dot(M, np.vstack((square_vertices.T, np.ones(4))))
    return transformed_vertices[:2].T