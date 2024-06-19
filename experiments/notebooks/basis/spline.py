import numpy as np

### Spline basis 2D
#%% Linear
def L_0(u):
    return 1 - u

def L_1(u):
    return u

def tensor_product_linear(u, v):
    """Compute the tensor product of 1D linear B-spline basis functions for 2D space."""
    basis_functions_1d_u = [L_0(u), L_1(u)]
    basis_functions_1d_v = [L_0(v), L_1(v)]
    basis_functions_2d = {}
    for I, fu in enumerate(basis_functions_1d_u):
        for J, fv in enumerate(basis_functions_1d_v):
            basis_functions_2d[(I, J)] = fu * fv
    return basis_functions_2d

#%% Quadratic
def Q_0(u):
    return (1 - u) ** 2 / 2

def Q_1(u):
    return (-2 * u ** 2 + 2 * u + 1) / 2

def Q_2(u):
    return u ** 2 / 2

def tensor_product_quadratic(u, v):
    """Compute the tensor product of 1D quadratic B-spline basis functions for 2D space."""
    basis_functions_1d_u = [Q_0(u), Q_1(u), Q_2(u)]
    basis_functions_1d_v = [Q_0(v), Q_1(v), Q_2(v)]
    basis_functions_2d = {}
    for I, fu in enumerate(basis_functions_1d_u):
        for J, fv in enumerate(basis_functions_1d_v):
            basis_functions_2d[(I, J)] = fu * fv
    return basis_functions_2d

#%% Cubic
def C_0(u):
    return (1 - u) ** 3 / 6

def C_1(u):
    return (3 * u ** 3 - 6 * u ** 2 + 4) / 6

def C_2(u):
    return (-3 * u ** 3 + 3 * u ** 2 + 3 * u + 1) / 6

def C_3(u):
    return u ** 3 / 6

def tensor_product_cubic(u, v):
    """Compute the tensor product of 1D cubic B-spline basis functions for 2D space."""
    basis_functions_1d_u = [C_0(u), C_1(u), C_2(u), C_3(u)]
    basis_functions_1d_v = [C_0(v), C_1(v), C_2(v), C_3(v)]
    basis_functions_2d = {}
    for I, fu in enumerate(basis_functions_1d_u):
        for J, fv in enumerate(basis_functions_1d_v):
            basis_functions_2d[(I, J)] = fu * fv
    return basis_functions_2d

### Flatten
def flatten_tensor_product(basis_dict):
    """
    Flatten the output of a tensor product of B-spline basis functions into a numpy array.

    Parameters:
    basis_dict : dict
        Dictionary with tuple keys (I, J) and values as the product of basis functions.

    Returns:
    flat_array : np.array
        Flat array containing the tensor product values.
    """
    # The size can be inferred from the keys as they will form a complete square/cube/etc. grid
    size = int(len(basis_dict) ** 0.5)  # Since it's 2D, we take the square root of the number of elements
    flat_array = np.zeros((size, size))

    for (i, j), value in basis_dict.items():
        flat_array[i, j] = value

    return flat_array.flatten()