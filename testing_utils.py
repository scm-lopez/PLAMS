"""A module with a number of testing-related uilities."""

# NOTE: the absolute import paths ensure that the following import
# syntax can be used (this would otherwise result in an ImportError):
#
# >>> from scm.plams import testing_utils
# >>> from testing_utils import blablabla

__all__ = ['assert_mol_eq', 'SMILES_TUPLE', 'MOL_TUPLE', 'MOL', 'HDF5_TMP', 'HDF5_READ']

from typing import Tuple, Optional
from pathlib import Path

import numpy as np

from scm.plams.mol.molecule import Molecule
from scm.plams._data import PDB_TUPLE

try:
    from scm.plams.interfaces.hdf5.h5py_mol import HDF5Mol
except ImportError as ex:
    H5PY_EX: Optional[ImportError] = ex
else:
    H5PY_EX = None
    __all__.append('HDF5_MOL')

try:
    from scm.plams.interfaces.molecule.rdkit import readpdb
except ImportError:
    def readpdb(file):  # type: ignore
        ret = Molecule(file)
        ret.guess_bonds()
        return ret

#: A tuple of SMILES strings.
#: Contains SMILES representations of :data:`MOL_TUPLE`
SMILES_TUPLE: Tuple[str, ...] = (
    'CCCCCCC(=O)O',
    'COC(=O)[O-]',
    'CCCCCCC(=O)[O-]',
    'CCCCCCC(=O)[O-]',
    'CCCOCCC(=O)O',
    'C[O-]',
    'CCCC[O-]',
    'CCCCCCCCP(CCCCCCCC)CCCCCCCC',
    'CNC(=O)[O-]',
    'CCCNCCC(=O)O',
    'CN',
    'CC(=O)[O-]',
    'CCC(=O)[O-]',
    'CCCOCCC(=O)[O-]',
    'C[NH-]',
    'O=C[O-]',
    'CO',
    'CC[O-]',
    'CCCOCCC(=O)[O-]',
    'CCO',
    'CCCNCCC(=O)[O-]',
    'CCCCCCCCN(CCCCCCCC)CCCCCCCC',
    'CCCNCCC(=O)[O-]'
)

#: A tuple of PLAMS Molecules.
MOL_TUPLE: Tuple[Molecule, ...] = tuple(readpdb(f) for f in PDB_TUPLE)

#: A PLAMS Molecule.
MOL = MOL_TUPLE[5]

#: A path to a temporary (to-be created) hdf5 file.
HDF5_TMP = Path('unit_tests') / 'test_files' / '.structures.hdf5'

#: A path to a read-only hdf5 file.
HDF5_READ = Path('unit_tests') / 'test_files' / 'structures.hdf5'


def assert_mol_eq(mol1: Molecule, mol2: Molecule) -> None:
    """Assert that two molecules are equivalent by comparing four properties.

    * Cartesian coordinates.
    * Bonds and bond orders.
    * Atomic numbers.
    * Lattice vectors.

    """
    np.testing.assert_allclose(mol1, mol2, err_msg="coords")
    np.testing.assert_allclose(mol1.bond_matrix(), mol2.bond_matrix(),
                               atol=0.5, err_msg="bonds")

    atnum1 = np.fromiter((at.atnum for at in mol1), dtype=int, count=len(mol1))
    atnum2 = np.fromiter((at.atnum for at in mol2), dtype=int, count=len(mol2))
    np.testing.assert_array_equal(atnum1, atnum2, err_msg="atnum")

    if not (np.any(mol1.lattice) or np.any(mol2.lattice)):
        return
    else:
        np.testing.assert_allclose(mol1.lattice, mol2.lattice, err_msg="lattice")


if H5PY_EX is None:
    #: A HDF5Mol instance.
    HDF5_MOL = HDF5Mol.from_molecules(MOL_TUPLE)
