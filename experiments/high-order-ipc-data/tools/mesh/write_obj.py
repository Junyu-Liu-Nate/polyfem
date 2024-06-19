# def write_obj(filename, V, E=None, F=None):
#     with open(filename, 'w') as f:
#         for v in V:
#             f.write("v {:.16f} {:.16f} {:.16f}\n".format(*v))
#         if E is not None:
#             for e in E:
#                 f.write("l {:d} {:d}\n".format(*(e + 1)))
#         if F is not None:
#             for face in F:
#                 f.write("f {:d} {:d} {:d}\n".format(*(face + 1)))

import numpy as np

def write_obj(filename, V, E=None, F=None):
    with open(filename, 'w') as f:
        for v in V:
            # Ensure the vertex has three components, appending 0 for the z-coordinate if necessary
            if len(v) == 2:  # Only x, y provided
                v = np.append(v, 0.0)  # Append 0 for the z-coordinate
            elif len(v) < 2:
                raise ValueError("Vertex must have at least two dimensions.")
            f.write("v {:.16f} {:.16f} {:.16f}\n".format(*v))
        
        if E is not None:
            for e in E:
                f.write("l {:d} {:d}\n".format(*(e + 1)))
        if F is not None:
            for face in F:
                f.write("f {:d} {:d} {:d}\n".format(*(face + 1)))