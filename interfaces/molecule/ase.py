from ...core.basemol import Molecule,Atom
from ...core.functions import add_to_class
from numpy import zeros as npz
from numpy import array as npa

__all__ = ['toASE', 'fromASE', 'readase', 'writease']
ase_present = False

try:
    from ase import Atom as aseAtom
    from ase import Atoms as aseAtoms
    ase_present = True
except ImportError:
    __all__ = []


if ase_present:
    @add_to_class(Molecule)
    def readase(self, f, geometry, **other):
        """Read Molecule using ASE engine

        The ``read`` function of the |Molecule| class passes a file descriptor into here, so in this case you must specify the *format* to be read by ASE.
        Otherwise use ``|Molecule|.readase()`` to avoid specifying the *format* manually.
        The ASE Atoms object then gets converted to a PLAMS Molecule and returned.
        All *other* options are passed to ``ASE.io.read()``.
        See https://wiki.fysik.dtu.dk/ase/ase/io/io.html on how to use it.

        NOTE: The *geometry* option is neglected. Use the corresponding ASE keyword to pick a certain geometry/image.
        The nomenclature of PLAMS and ASE is incompatible for reading multiple geometries, make sure that you only read single geometries with ASE! Reading multiple geometries is not supported, each frame needs to be read individually.
        """
        try:
            from ase import io as aseIO
        except ImportError:
            raise MoleculeError('Asked for ASE IO engine but could not load ASE.io module')

        aseMol = aseIO.read(f, **other)
        mol = fromASE(aseMol)
        #update self with the molecule read without overwriting e.g. settings
        self += mol
        #lattice does not survive soft update
        self.lattice = mol.lattice
        return


    @add_to_class(Molecule)
    def writease(self, f, **other):
        """Write molecular coordinates using ASE engine.

        The ``write`` function of the |Molecule| class passes a file descriptor into here, so in this case you must specify the *format* to be written by ASE.
        Otherwise use ``|Molecule|.writease()`` to avoid specifying the *format* manually.
        All *other* options are passed to ``ASE.io.write()``.
        See https://wiki.fysik.dtu.dk/ase/ase/io/io.html on how to use it.


        These two write the same content to the respective files:

        >>> molecule.write('filename.anyextension', outputformat='ase', format='gen')
        >>> molecule.writease('filename.gen')

        But this one will write a JSON dump because ASE does not know the filetype from a file descriptor passed by PLAMS:
        >>> molecule.writease('filename.gen', outputformat='ase')
        """

        aseMol = toASE(self)
        aseMol.write(f, **other)
        return

    Molecule._readformat['ase'] = Molecule.readase
    Molecule._writeformat['ase'] = Molecule.writease


def toASE(molecule):
    """Convert a PLAMS |Molecule| to an ASE molecule (``ase.Atoms`` instance). Translate coordinates, atomic numbers, and lattice vectors (if present). The order of atoms is preserved."""
    aseMol = aseAtoms()

    #iterate over PLAMS atoms
    for atom in molecule:

        #check if coords only consists of floats or ints
        if not all(isinstance(x, (int,float)) for x in atom.coords):
            raise ValueError("Non-Number in Atomic Coordinates, not compatible with ASE")

        #append atom to aseMol
        aseMol.append(aseAtom(atom.atnum, atom.coords))

    #get lattice info if any
    lattice = npz((3,3))
    pbc = [False,False,False]
    for i,vec in enumerate(molecule.lattice):

        #check if lattice only consists of floats or ints
        if not all(isinstance(x, (int,float)) for x in vec):
            raise ValueError("Non-Number in Lattice Vectors, not compatible with ASE")

        pbc[i] = True
        lattice[i] = npa(vec)

    #save lattice info to aseMol
    if any(pbc):
        aseMol.set_pbc(pbc)
        aseMol.set_cell(lattice)

    return aseMol




def fromASE(molecule):
    """Convert an ASE molecule to a PLAMS |Molecule|. Translate coordinates, atomic numbers, and lattice vectors (if present). The order of atoms is preserved."""
    plamsMol = Molecule()

    #iterate over ASE atoms
    for atom in molecule:
        #add atom to plamsMol
        plamsMol.add_atom(Atom(atnum=atom.number, coords=tuple(atom.position)))

    #add Lattice if any
    if any(molecule.get_pbc()):
        lattice = []
        #loop over three booleans
        for i,boolean in enumerate(molecule.get_pbc().tolist()):
            if boolean:
                lattice.append(tuple(molecule.get_cell()[i]))

        #write lattice to plamsMol
        plamsMol.lattice = lattice.copy()


    return plamsMol
