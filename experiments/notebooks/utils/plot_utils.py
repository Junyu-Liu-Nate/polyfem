import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# def plot_mesh(vertices, faces):
#     """
#     Plot a mesh given vertices and faces using matplotlib.

#     Parameters:
#     vertices : np.array
#         Array of vertices, where each row represents a vertex.
#     faces : np.array
#         Array of faces, where each row contains indices into the vertices array that form a face.
#     """
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')

#     # Plotting each face
#     for face in faces:
#         verts = vertices[face]  # Get the actual vertices for the current face
#         tri = Poly3DCollection([verts], edgecolor='k', alpha=0.5, linewidths=1.5)
#         tri.set_facecolor([0.5, 0.5, 0.65])  # Set the face color (optional)
#         ax.add_collection3d(tri)

#         # Additionally, plot vertices
#         ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], color='r')

#     # Setting the aspect ratio
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     max_range = np.array([vertices[:, 0].max()-vertices[:, 0].min(), 
#                           vertices[:, 1].max()-vertices[:, 1].min(), 
#                           vertices[:, 2].max()-vertices[:, 2].min()]).max() / 2.0
#     mid_x = (vertices[:, 0].max()+vertices[:, 0].min()) * 0.5
#     mid_y = (vertices[:, 1].max()+vertices[:, 1].min()) * 0.5
#     mid_z = (vertices[:, 2].max()+vertices[:, 2].min()) * 0.5
#     ax.set_xlim(mid_x - max_range, mid_x + max_range)
#     ax.set_ylim(mid_y - max_range, mid_y + max_range)
#     ax.set_zlim(mid_z - max_range, mid_z + max_range)

#     plt.show()

def plot_mesh(vertices, faces):
    """
    Plot a mesh given vertices and faces using matplotlib. Handles both 2D and 3D vertices.

    Parameters:
    vertices : np.array
        Array of vertices, where each row represents a vertex.
    faces : np.array
        Array of faces, where each row contains indices into the vertices array that form a face.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Ensure vertices are in 3D by padding with zeros if necessary
    if vertices.shape[1] == 2:
        # If vertices are 2D, add a zero z-coordinate
        vertices = np.hstack([vertices, np.zeros((vertices.shape[0], 1))])

    # Plotting each face
    for face in faces:
        verts = vertices[face]  # Get the actual vertices for the current face
        tri = Poly3DCollection([verts], edgecolor='k', alpha=0.5, linewidths=1.5)
        tri.set_facecolor([0.5, 0.5, 0.65])  # Set the face color (optional)
        ax.add_collection3d(tri)

        # Additionally, plot vertices
        ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], color='r')

    # Setting the aspect ratio
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    max_range = np.array([vertices[:, 0].max()-vertices[:, 0].min(), 
                          vertices[:, 1].max()-vertices[:, 1].min(), 
                          vertices[:, 2].max()-vertices[:, 2].min()]).max() / 2.0
    mid_x = (vertices[:, 0].max()+vertices[:, 0].min()) * 0.5
    mid_y = (vertices[:, 1].max()+vertices[:, 1].min()) * 0.5
    mid_z = (vertices[:, 2].max()+vertices[:, 2].min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    plt.show()