import numpy as np

from .units import Units

__all__ = ['rotation_matrix', 'axis_rotation_matrix']

def rotation_matrix(vec1, vec2):
    """
    Calculate the rotation matrix rotating *vec1* to *vec2*. Vectors can be any containers with 3 numerical values. They don't need to be normalized. Returns 3x3 numpy array.
    """
    vec1, vec2 = np.array(vec1), np.array(vec2)
    a = vec1/np.linalg.norm(vec1)
    b = vec2/np.linalg.norm(vec2)
    v1,v2,v3 = np.cross(a,b)
    M = np.array([[0, -v3, v2], [v3, 0, -v1], [-v2, v1, 0]])
    return (np.identity(3) + M + np.dot(M,M)/(1+np.dot(a,b)))


def axis_rotation_matrix(vector, angle, unit='radian'):
    """
    Calculate the rotation matrix rotating along the *vector* by *angle* expressed in *unit*.

    *vector* can be any container with 3 numerical values. They don't need to be normalized. A positive angle denotes counterclockwise rotation, when looking along *vector*. Returns 3x3 numpy array.
    """

    vector /= np.linalg.norm(vector)
    v0, v1, v2 = vector

    W = np.array([[0, -v2, v1],
                  [v2, 0, -v0],
                  [-v1, v0, 0]])

    angle = Units.convert(angle, unit, 'radian')
    a1 = np.sin(angle)
    a2 = 1.0 - np.cos(angle)

    return np.identity(3) + a1 * W + a2 * W@W
