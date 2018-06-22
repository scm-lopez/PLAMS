import numpy as np

__all__ = ['rotation_matrix']

def rotation_matrix(vec1, vec2):
    """
    Calculates the rotation matrix rotating *vec1* to *vec2*. Vectors can be any containers with 3 numerical values. They don't need to be normalized. Returns 3x3 numpy array.
    """
    a = np.array(vec1)/np.linalg.norm(vec1)
    b = np.array(vec2)/np.linalg.norm(vec2)
    v1,v2,v3 = np.cross(a,b)
    M = np.array([[0, -v3, v2], [v3, 0, -v1], [-v2, v1, 0]])
    return (np.identity(3) + M + np.dot(M,M)/(1+np.dot(a,b)))