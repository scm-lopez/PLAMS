import numpy as np

from .atom import Atom
from .molecule import Molecule
from ..core.private import sha256
from ..core.functions import add_to_class


__all__ = ['label']


possible_flags = ['BO', 'RS', 'EZ', 'DH', 'CO', 'H2']


def twist(v1, v2, v3, tolerance=1e-2):
    """
    Given 3 vectors in 3D space measure their "chirality" with *tolerance*.

    Returns a pair. The first element is an integer number measuring the orientation (clockwise vs counterclockwise) of *v1* and *v3* while looking along *v2*. Values 1 and -1 indicate this case and the second element of returned pair is ``None``. Value 0 indicates that *v1*, *v2*, and *v3* are coplanar, and the second element of the returned pair is indicating if two turns made by going *v1*->*v2*->*v3* are the same (left-left, right-right) or the opposite (left-right, right-left).
    """
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    v3 /= np.linalg.norm(v3)
    x = np.dot(v3, np.cross(v1,v2))
    if abs(x) <= tolerance:  #v1, v2, v3 are coplanar
        return 0, int(np.sign(np.dot(np.cross(v1,v2), np.cross(v2,v3))))
    return int(np.sign(x)), None


def unique_atoms(atomlist):
    """Filter *atomlist* (list or |Molecule|) for atoms with unique ``IDname``."""
    d = {}
    for atom in atomlist:
        if atom.IDname not in d:
            d[atom.IDname] = 0
        d[atom.IDname] += 1
    return [atom for atom in atomlist if d[atom.IDname] == 1]


def initialize(molecule):
    """Initialize atom labeling algorithm by setting ``IDname`` and ``IDdone`` attributes for all atoms in *molecule*."""
    for at in molecule:
        at.IDname = at.symbol
        at.IDdone = False


def clear(molecule):
    """Remove ``IDname`` and ``IDdone`` attributes from all atoms in *molecule*."""
    for at in molecule:
        del at.IDname
        del at.IDdone


def iterate(molecule, flags):
    """Perform one iteration of atom labeling alogrithm.

    First, mark all atoms that are unique and have only unique neighbors as "done". Then calculate new label for each atom that is not done. Return True if the number of different atom labels increased during this iteration.
    """
    names = len(set(at.IDname for at in molecule))
    unique = set(unique_atoms(molecule))

    for atom in molecule:
        if atom in unique and all(N in unique for N in atom.neighbors()):
            atom.IDdone = True
        if not atom.IDdone:
            atom.IDnew = new_name(atom, flags)

    for atom in molecule:
        if not atom.IDdone:
            atom.IDname = atom.IDnew

    new_names = len(set(atom.IDname for atom in molecule))
    return new_names > names #True means this iteration increased the number of distinct names


def new_name(atom, flags):
    """Compute new label for *atom*.

    The new label is based on the existing label of *atom*, labels of all its neighbors and (possibly) some additional conformational information. The labels of neighbors are not obtained directly by reading neighbor's ``IDname`` but rather by a process called "knocking". The *atom* knocks all its bonds. Each knocked bond returns an identifier describing the atom on the other end of the bond. The identifier is composed of knocked atom's ``IDname`` together with some additional information desribing the character of the bond and knocked atom's spatial environment. The exact behavior of this mechanism is adjusted by the contents of *flags* dictionary (see :func:`label_atoms` for details).

    TODO CO
    """

    knocks = [knock(atom, bond, flags) for bond in atom.bonds]
    knocks.sort(key=lambda x: x[0])

    more = []
    if flags['RS'] and len(knocks) == 4 and len(set(i[0] for i in knocks)) == 4:
        v1 = atom.vector_to(knocks[0][1])
        v2 = atom.vector_to(knocks[1][1])
        v3 = atom.vector_to(knocks[2][1])
        more.append(str(twist(v1,v2,v3)))
    return sha256('|'.join([atom.IDname] + [i[0] for i in knocks] + more))


def knock(A, bond, flags):
    S = bond.other_end(A)
    ret = S.IDname

    if flags['BO'] and bond.order != 1:
        ret += 'BO' + str(bond.order)

    if flags['EZ'] or flags['DH']:
        S_unique = unique_atoms(S.neighbors())
        if A in S_unique: #*A* is a unique neighbor of *S*
            S_unique.remove(A)
            S_unique.sort(key=lambda x: x.IDname)
            for b in S.bonds:
                N = b.other_end(S)
                if N in S_unique:
                    N_unique = unique_atoms(N.neighbors())
                    if S in N_unique:
                        N_unique.remove(S)
                    N_unique.sort(key=lambda x: x.IDname)
                    if N_unique:
                        F = N_unique[0]
                        v1 = A.vector_to(S)
                        v2 = S.vector_to(N)
                        v3 = N.vector_to(F)
                        t = twist(v1,v2,v3)
                        if flags['DH']:
                            ret += 'DH' + str(t)
                        elif flags['EZ'] and b.order == 2 and t[0] == 0:
                            #A-S=N-F are coplanar
                            ret += 'EZ' + str(t[1])
                        break

    return (ret, S)


def label_atoms(molecule, **kwargs):
    """Label atoms in *molecule*.

    Boolean keyword arguments:
    *   *BO* -- include bond orders
    *   *RS* -- include R/S stereoisomerism
    *   *EZ* -- include E/Z stereoisomerism
    *   *DH* -- include some dihedrals to detect alkane rotamers and syn-anti conformers in cycloalkanes
    *   *CO* -- include more spatial info to detect different conformation of coordination complexes (flat square, octaedr etc.)

    Diherdals considered with *DH* are all the dihedrals A-B-C-D such that A is a unique neighbor or B and D is a unique neighbor of C.

    TODO CO
    """
    initialize(molecule)
    while iterate(molecule, kwargs):
        pass
    return molecule


def molecule_name(molecule):
    """Compute the label of the whole *molecule* based on ``IDname`` attributes of all the atoms."""
    names = [atom.IDname for atom in molecule]
    names.sort()
    return sha256(' '.join(names))


def label(molecule, level=0, specify=None, keep_labels=False):
    """Compute the label of *molecule* using chosen *level* of detail.

    Possible levels are:

        *   **0**: only direct connectivity is considered, without bond orders (in other words, treat all bonds as single bonds)
        *   **1**: use connectivity and bond orders
        *   **2**: use connectivity, bond orders and some spatial information to distinguish R/S and E/Z isomers
        *   **3**: use all above, plus more spatial information to distinguish different rotamers and different types of coordination complexes

    For more precise control of what is taken into account while computing the label you can use *specify* argument. It should be a list of two-letter boolean flags recognized by :func:`label_atoms`. If you use *specify*, the *level* argument is ignored.

    This function, by default, erases ``IDname`` attributes of all atoms at the end. You can change this behavior with *keep_labels* argument.

    .. note::

        This method is a new PLAMS feature and its still somewhat experimental. The exact details of the algorithm can, and probably will, change in future. You are more than welcome to provide any feedback or feature requests.

    """
    flags = {i:False for i in possible_flags}
    if specify:
        for i in specify:
            flags[i] = True
    else:
        if level >= 1:
            flags['BO'] = True
        if level >= 2:
            flags['RS'] = True
            flags['EZ'] = True
        if level >= 3:
            flags['DH'] = True
            flags['CO'] = True

    if len(molecule.bonds) == 0:
        molecule.guess_bonds()

    label_atoms(molecule, **flags)
    ret = molecule_name(molecule)
    if not keep_labels:
        clear(molecule)
    return ret


#@add_to_class(Molecule)
#def identify(self, levels=(0,1,2,3)):
#    """TODO
#
#
#    """
#    return label(self, level=levels) if isinstance(levels, int) else tuple(label(self, level=i) for i in levels)




