from ..core.basemol import Molecule,Atom
from numpy import zeros as npz
from numpy import array as npa

__all__ = ['toASE','fromASE']

try:
    from ase import Atom as aseAtom
    from ase import Atoms as aseAtoms
except:
    __all__ = []



def toASE(molecule):
    """
    Converts a PLAMS molecule to an ASE molecule. The following attributes are converted, conserving the order of atoms:
    -Coordinates
    -Atomic Number (Symbol is derived automaticaly)
    -Periodicity and Cell Vectors
    """
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
    """
    Converts an ASE molecule to a PLAMS molecule. The following attributes are converted, conserving the order of atoms:
    -Coordinates
    -Atomic Number (Symbol is derived automaticaly)
    -Periodicity and Cell Vectors
    """
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
