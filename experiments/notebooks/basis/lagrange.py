import numpy as np

def linear_triangle_basis(x, y):
    """Compute the linear basis functions for a right-angled triangle."""
    phi1 = 1 - x - y
    phi2 = x
    phi3 = y
    return np.array([phi1, phi2, phi3])

# def linear_triangle_basis(x, y, vertices):
#     """
#     Compute the linear basis functions for a triangle element.
    
#     Parameters:
#     x, y : float
#         The point at which to evaluate the basis functions
#     vertices : np.array
#         A 3x2 array containing the coordinates of the triangle vertices
        
#     Returns:
#     np.array : The values of the three basis functions at (x, y)
#     """
#     # Compute the area of the triangle
#     area = 0.5 * np.abs(np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0]))
    
#     # Initialize basis functions
#     phi = np.zeros(3)
    
#     for i in range(3):
#         j = (i + 1) % 3
#         k = (i + 2) % 3
        
#         # Compute the coefficients a, b, c for the linear function
#         a = (vertices[j, 1] - vertices[k, 1]) / (2 * area)
#         b = (vertices[k, 0] - vertices[j, 0]) / (2 * area)
#         c = (vertices[j, 0] * vertices[k, 1] - vertices[k, 0] * vertices[j, 1]) / (2 * area)
        
#         # Evaluate the basis function
#         phi[i] = a * x + b * y + c
    
#     return phi

# def gradient_linear_triangle_basis(vertices):
#     """
#     Compute the gradients of the linear basis functions for a triangle element.
    
#     Parameters:
#     vertices : np.array
#         A 3x2 array containing the coordinates of the triangle vertices
        
#     Returns:
#     np.array : A 3x2 array containing the gradients of the three basis functions
#     """
#     # Compute the area of the triangle
#     area = 0.5 * np.abs(np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0]))
    
#     # Initialize gradients
#     grad_phi = np.zeros((3, 2))
    
#     for i in range(3):
#         j = (i + 1) % 3
#         k = (i + 2) % 3
        
#         # Compute the coefficients a, b for the gradient
#         grad_phi[i, 0] = (vertices[j, 1] - vertices[k, 1]) / (2 * area)  # d/dx
#         grad_phi[i, 1] = (vertices[k, 0] - vertices[j, 0]) / (2 * area)  # d/dy
    
#     return grad_phi

# # Example usage:
# if __name__ == "__main__":
#     # Define a triangle
#     vertices = np.array([[0, 0], [1, 0], [0, 1]])
    
#     # Evaluate basis functions at a point
#     x, y = 0.3, 0.3
#     basis_values = linear_triangle_basis(x, y, vertices)
#     print("Basis function values at ({}, {}):".format(x, y))
#     print(basis_values)
#     print("Sum of basis functions:", np.sum(basis_values))
    
#     # Compute gradients
#     gradients = gradient_linear_triangle_basis(vertices)
#     print("\nBasis function gradients:")
#     print(gradients)