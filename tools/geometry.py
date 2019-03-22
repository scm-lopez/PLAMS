import numpy as np
try:
    from scipy.spatial.distance import cdist
    scipy_present = True
except ImportError:
    scipy_present = False

from .units import Units

__all__ = ['rotation_matrix', 'axis_rotation_matrix', 'distance_array', 'dihedral']

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


def distance_array(array1, array2):
    """Calculates distance between each pair of points in *array1* and *array2*. Returns 2D ``numpy`` array.

    Uses fast ``cdist`` function if ``scipy`` is present, otherwise falls back to slightly slower ``numpy`` loop. Arguments should be 2-dimensional ``numpy`` arrays with the same second dimension. If *array1* is A x N and *array2* is B x N, the returned array is A x B.
    """
    return cdist(array1, array2) if scipy_present else np.array([np.linalg.norm(i - array2, axis=1) for i in array1])



def dihedral(p1, p2, p3, p4, unit='radian'):
    """Calculate the value of diherdal angle formed by points *p1*, *p2*, *p3* and *p4* in a 3D space. Arguments can be any containers with 3 numerical values, also instances of |Atom|. Returned value is always non-negative, measures the angle clockwise (looking along *p2-p3* vector) and is expressed in *unit*."""
    p1 = np.array([*p1], dtype=float)
    p2 = np.array([*p2], dtype=float)
    p3 = np.array([*p3], dtype=float)
    p4 = np.array([*p4], dtype=float)

    b0 = p1 - p2
    b1 = p3 - p2
    b2 = p4 - p3

    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    ret = np.arctan2(y, x)
    ret = 2*np.pi+ret if ret < 0 else ret
    return Units.convert(ret, 'radian', unit)
