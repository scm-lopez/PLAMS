"""A module for constructing array-representations of .pdb files."""

import sys
from types import MappingProxyType
from typing import ClassVar, Mapping, Iterable, Optional, Iterator

import numpy as np

from .h5py_mol import HDF5Mol, InfoTuple
from ._annotations import DtypeLike, ArrayLike, _AttrNames, _GetTuple, _SetTuple
from ...mol.molecule import Molecule
from ...mol.atom import Atom

if sys.version_info[1] >= 7:
    from builtins import dict as OrderedDict
else:
    from collections import OrderedDict

try:
    import h5py
except ImportError:
    __all__ = []
else:
    __all__ = ['HDF5Pdb', 'ATOMS_DTYPE_PDB']


def _get_atom_info(at: Atom, i: int) -> _GetTuple:
    """Helper function for :meth:`HDF5Pdb.get_atom_info`."""
    symbol = at.symbol

    pdb = at.properties.get('pdb_info', {})
    pdb_get = pdb.get
    return (
        pdb_get('IsHeteroAtom', False),  # type: ignore
        pdb_get('SerialNumber', i),
        pdb_get('Name', symbol),
        pdb_get('AltLoc', ''),
        pdb_get('ResidueName', ''),
        pdb_get('ChainId', ''),
        pdb_get('ResidueNumber', 1),
        pdb_get('InsertionCode', ''),
        *at.coords,
        pdb_get('Occupancy', 1),
        pdb_get('TempFactor', 0),
        symbol,
        at.properties.get('charge', 0)
    )


def _set_atom_info(atom_array: Iterable[np.void]) -> Iterator[_SetTuple]:
    """Helper function for :meth:`HDF5Pdb.set_atom_info`."""
    for generic in atom_array:
        IsHeteroAtom, SerialNumber, Name, AltLoc, ResidueName, ChainId, ResidueNumber, InsertionCode, x, y, z, Occupancy, TempFactor, symbol, charge = generic.item()  # noqa: E501
        _pdb_info = {
            'IsHeteroAtom': IsHeteroAtom,
            'SerialNumber': SerialNumber,
            'Name': Name.decode(),
            'AltLoc': AltLoc.decode(),
            'ResidueName': ResidueName.decode(),
            'InsertionCode': InsertionCode.decode(),
            'ChainId': ChainId.decode(),
            'ResidueNumber': ResidueNumber,
            'Occupancy': Occupancy,
            'TempFactor': TempFactor
        }

        properties = {
            'charge': charge,
            'pdb_info': _pdb_info
        }
        yield properties, (x, y, z), symbol.decode()  # type: ignore


_ATOMS_MAPPING = OrderedDict(
    IsHeteroAtom='bool',
    SerialNumber='int16',
    Name='S4',
    AltLoc='S1',
    ResidueName='S3',
    ChainId='S1',
    ResidueNumber='int16',
    InsertionCode='S1',
    x='float64',
    y='float64',
    z='float64',
    Occupancy='float64',
    TempFactor='float64',
    symbol='S4',
    charge='int8'
)
#: The dtype of :attr:`HDF5Pdb.atoms`.
ATOMS_DTYPE_PDB = np.dtype(list(_ATOMS_MAPPING.items()))

_INFO = HDF5Mol.INFO.copy()
_INFO['atoms'] = _INFO['atoms']._replace(dtype=ATOMS_DTYPE_PDB)


class HDF5Pdb(HDF5Mol):
    """A :class:`~scm.plams.interfaces.hdf5.h5py_mol.HDF5Mol` subclass for the storage of :func:`~scm.plams.interfaces.molecule.rdkit.readpdb`- / :func:`~scm.plams.interfaces.molecule.rdkit.writepdb`-compatible molecules.

    Identical in all respect to its baseclass with the exception of :attr:`HDF5Pdb.atoms`,
    whose data type has been expanded with additional fields (see :data:`ATOMS_DTYPE_PDB`).

    """

    INFO: ClassVar[Mapping[_AttrNames, InfoTuple]] = MappingProxyType(_INFO)

    # The `atoms` property is redefined here for to sole reason to update its docstring

    @property
    def atoms(self) -> np.ndarray:
        """:class:`numpy.ndarray`, dtype :data:`ATOMS_DTYPE_PDB`, shape :math:`(n,)` : A ragged ndarray for keeping track of all atom-related information.

        Supports both getting and setting.
        Defaults to :func:`numpy.zeros` if :data:`None` is passed.

        """
        return self._atoms

    @atoms.setter
    def atoms(self, value: Optional[ArrayLike]) -> None:
        self._atoms = self._parse_array('atoms', value)

    @staticmethod
    def get_atom_info(mol: Molecule, dtype: DtypeLike) -> np.ndarray:
        """Helper method for :meth:`HDF5Pdb.from_molecules<HDF5Mol.from_molecules>`.

        Get data from all atoms in **mol** at the specified **dtype**.

        Parameters
        ----------
        mol : :class:`~scm.plams.mol.molecule.Molecule`
            A PLAMS molecule.
        dtype : :class:`numpy.dtype`
            The data type of the to-be returned array.
            See :data:`HDF5Pdb.INFO["atoms"].dtype<ATOMS_DTYPE_PDB>`.

        Returns
        -------
        :class:`numpy.ndarray`
            A 1-D array with the specified **dtype**.

        """
        count = len(mol.atoms)
        iterator = (_get_atom_info(at, i) for i, at in enumerate(mol, 1))
        return np.fromiter(iterator, dtype, count=count)

    @staticmethod
    def set_atom_info(mol: Molecule, data: np.ndarray) -> None:
        """Helper method for :meth:`HDF5Pdb.to_molecules()<HDF5Mol.to_molecules>`.

        Update all atoms in **mol** with **data**.

        Parameters
        ----------
        mol : :class:`~scm.plams.mol.molecule.Molecule`
            A PLAMS molecule.
        data : :class:`numpy.ndarray`
            The to-be set data.
            See :data:`HDF5Pdb.INFO["atoms"].dtype<ATOM_DTYPE_PDB>` for the data type.

        """
        iterator = _set_atom_info(data)
        for atom, (properties, coords, symbol) in zip(mol.atoms, iterator):
            atom.properties.update(properties)
            atom.coords = coords
            atom.symbol = symbol
