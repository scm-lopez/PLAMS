from ..core.errors import MoleculeError
from ..core.settings import Settings
from ..tools.units import Units


__all__ = ['Bond']


class Bond:
    """A class representing a bond between two atoms.

    An instance of this class has the following attributes:

    *   ``atom1`` and ``atom2`` -- two instances of |Atom| that form this bond
    *   ``order`` -- order of the bond. It is either an integer number or the floating point value stored in ``Bond.AR``, indicating an aromatic bond
    *   ``mol`` -- |Molecule| this bond belongs to
    *   ``properties`` -- |Settings| instance storing all other  information about this bond (initially it is populated with *\*\*other*)

    .. note::

        Newly created bond is **not** added to ``atom1.bonds`` or ``atom2.bonds``. Storing information about |Bond| in |Atom| is relevant only in the context of the whole |Molecule|, so this information is updated by :meth:`~Molecule.add_bond`.

    """
    AR = 1.5
    def __init__(self, atom1=None, atom2=None, order=1, mol=None, **other):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        self.mol = mol
        self.properties = Settings(other)


    def __str__(self):
        """Return a string representation of this bond."""
        return '({})--{:1.1f}--({})'.format(str(self.atom1).strip(), self.order, str(self.atom2).strip())


    def __iter__(self):
        """Iterate over bonded atoms (``atom1`` first, then ``atom2``)."""
        yield self.atom1
        yield self.atom2


    def is_aromatic(self):
        """Check if this bond is aromatic."""
        return self.order == Bond.AR


    def length(self, unit='angstrom'):
        """Return bond length, expressed in *unit*."""
        return self.atom1.distance_to(self.atom2, result_unit=unit)


    def other_end(self, atom):
        """Return the atom on the other end of this bond with respect to *atom*. *atom* has to be one of the atoms forming this bond, otherwise an exception is raised.
        """
        if atom is self.atom1:
            return self.atom2
        elif atom is self.atom2:
            return self.atom1
        else:
            raise MoleculeError('Bond.other_end: invalid atom passed')


    def resize(self, atom, length, unit='angstrom'):
        """Change the length of the bond to *length*.

        One of two atoms forming this bond is moved along the bond in such a way that the new length is *length*, in *unit* (direction of the bond in space does not change). Atom indicated by *atom* has to be one of bond's atoms and it is the atom that is **not** moved.
        """
        ratio = 1.0 - Units.convert(length, unit, 'angstrom')/self.length()
        moving = self.other_end(atom)
        moving.translate(tuple(i*ratio for i in moving.vector_to(atom)))


