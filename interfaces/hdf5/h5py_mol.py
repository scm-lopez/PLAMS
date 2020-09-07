"""A module for constructing array-representations of molecules."""

import sys
import textwrap
import warnings
from types import MappingProxyType
from typing import NamedTuple, Optional
from itertools import repeat, zip_longest

import numpy as np

from ...mol.molecule import Molecule
from ...mol.atom import Atom
from ...mol.bond import Bond
from ._private import _parse_index, _alen, _to_array

# Normal dictionaries are ordered starting from python 3.7
if sys.version_info[1] >= 7:
    from builtins import dict as OrderedDict
else:
    from collections import OrderedDict

try:
    import h5py
except ImportError:
    __all__ = []
else:
    __all__ = ['HDF5Mol', 'InfoTuple', 'ATOMS_DTYPE', 'BONDS_DTYPE', 'MOL_DTYPE', 'SCALE_DTYPE']


class InfoTuple(NamedTuple):
    """A :func:`namedtuple<collections.namedtuple>` with three fields.

    Used as values in the :attr:`HDF5Mol.INFO` dictionary.
    """

    #: :class:`str`: The name of the corresponding :class:`h5py.Dataset`.
    name: str

    #: :class:`numpy.dtype`: The datatype of the array and corresponding Dataset.
    dtype: np.dtype

    #: :class:`bool`: Whether or :attr:`dtype` should be wrapped in a HDF5 variable-length datatype (see :func:`h5py.vlen_dtype`).
    vlen: bool = False

    #: :class:`int`: The dimensionality of the array and corresponding Dataset.
    ndim: int = 1

    #: :class:`str`, optional: The memory layout of the array.
    order: Optional[str] = None


class _FalseType:
    def __eq__(self, b): return False


#: An object which allways returns False when comparing.
#: Used to avoid :exc:`FutureWarning`s when comparing structured arrays to scalars.
_False = _FalseType()


_ATOMS_MAPPING = OrderedDict(
    symbol='S2',
    x='float64',
    y='float64',
    z='float64'
)
#: The dtype of :attr:`HDF5Mol.atoms`.
ATOMS_DTYPE = np.dtype(list(_ATOMS_MAPPING.items()))


_BONDS_MAPPING = OrderedDict(
    atom1='int64',
    atom2='int64',
    order='float64'
)
#: The dtype of :attr:`HDF5Mol.bonds`.
BONDS_DTYPE = np.dtype(list(_BONDS_MAPPING.items()))


_MOL_MAPPING = OrderedDict(
    lattice=('float64', (3, 3))
)
#: The dtype of :attr:`HDF5Mol.mol`.
MOL_DTYPE = np.dtype(list(_MOL_MAPPING.items()))


#: The (default) dtype of :attr:`HDF5Mol.scale`.
SCALE_DTYPE = np.dtype('int64')


class HDF5Mol:
    """An struct for holding array-like representions of PLAMS molecules.

    The :class:`HDF5Mol` class serves as an (intermediate) container
    facilitating the storage and interconversion
    between PLAMS molecules and the :mod:`h5py` interface.

    The methods implemented in this class can roughly be divided into four categories:

    * Molecule-interconversion: :meth:`~HDF5Mol.to_molecules` &
      :meth:`~HDF5Mol.from_molecules`.

    * Hdf5-interconversion: :meth:`~HDF5Mol.create_group`,
      :meth:`~HDF5Mol.to_hdf5` & :meth:`~HDF5Mol.from_hdf5`.

    * Set operatorations: :meth:`~HDF5Mol.intersection`,
      :meth:`~HDF5Mol.difference`, :meth:`~HDF5Mol.symmetric_difference` &
      :meth:`~HDF5Mol.union`.

    * Miscellaneous: :meth:`~HDF5Mol.keys`, :meth:`~HDF5Mol.values`,
      :meth:`~HDF5Mol.items`, :meth:`~HDF5Mol.__getitem__`,
      :meth:`~HDF5Mol.__len__`, :meth:`~HDF5Mol.copy` &
      :meth:`~HDF5Mol.concatenate`.

    Attributes
    ----------
    INFO : :class:`Mapping[str, InfoTuple]<typing.Mapping>`, class variable
        An ordered read-only dictionary with attribute names as keys and namedtuples as values.
        Used for describing the dataset names, data type and array type of all
        arrays embedded within the :class:`HDF5Mol` class.
        See the :class:`InfoTuple` namedtuple for more details.

    """

    __slots__ = ('__weakref__', '_atoms', '_bonds', '_mol', '_scale')

    INFO = MappingProxyType(OrderedDict(
        atoms=InfoTuple('atoms', ATOMS_DTYPE, vlen=True),
        bonds=InfoTuple('bonds', BONDS_DTYPE, vlen=True),
        mol=InfoTuple('mol', MOL_DTYPE, vlen=False),
        scale=InfoTuple('scale', SCALE_DTYPE, vlen=False)
    ))

    @property
    def atoms(self):
        """:class:`numpy.ndarray`, dtype :data:`ATOMS_DTYPE`, shape :math:`(n,)` : A ragged ndarray for keeping track of all atom-related information.

        Supports both getting and setting.
        Defaults to :func:`numpy.zeros` if :data:`None` is passed.

        """  # noqa: E501
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        size = max(_alen(i) for i in self.values())
        self._atoms = self._parse_array('atoms', value, size)

    @property
    def bonds(self):
        """:class:`numpy.ndarray`, dtype :data:`BONDS_DTYPE`, shape :math:`(n,)` : A ragged ndarray for keeping track of all bond-related information.

        Supports both getting and setting.
        Defaults to :func:`numpy.zeros` if :data:`None` is passed.

        """  # noqa: E501
        return self._bonds

    @bonds.setter
    def bonds(self, value):
        size = max(_alen(i) for i in self.values())
        self._bonds = self._parse_array('bonds', value, size)

    @property
    def mol(self):
        """:class:`numpy.ndarray`, dtype :data:`MOL_DTYPE`, shape :math:`(n,)` : An ndarray for keeping track of all molecule-related information.

        Supports both getting and setting.
        Defaults to :func:`numpy.zeros` if :data:`None` is passed.

        """  # noqa: E501
        return self._mol

    @mol.setter
    def mol(self, value):
        size = max(_alen(i) for i in self.values())
        self._mol = self._parse_array('mol', value, size)

    @property
    def scale(self):
        """:class:`numpy.ndarray`, dtype :data:`SCALE_DTYPE`, shape :math:`(n,)` : An ndarray containing unique molecular identifiers.

        Supports both getting and setting.
        Defaults to :func:`numpy.arange` if :data:`None` is passed.

        Notes
        -----
        The data type of this array is variable, depending on the passed value.

        """  # noqa: E501
        return self._scale

    @scale.setter
    def scale(self, value):
        size = max(_alen(i) for i in self.values())
        self._scale = self._parse_scale('scale', value, size)

    def __init_subclass__(cls) -> None:
        """Customize the creation of :class:`HDF5Mol` subclasses.

        Checks if a number of arrays have their dimensionality set to ``1``.
        """
        name = 'scale'
        if getattr(cls.INFO.get(name), 'ndim', 1) != 1:
            raise ValueError(f"Invalid dimensionality for {cls.__name__}.{name}; "
                             f"{name!r} requires 'ndim=1'")

    def __init__(self, atoms=None, bonds=None, mol=None, scale=None, validate=True, copy=True):
        """Initialize an instance.

        Parameters
        ----------
        atoms : array-like, shape :math:`(n,)`, optional
            See :attr:`atoms`.
        bonds : array-like, shape :math:`(n,)`, optional
            See :attr:`bonds`.
        mol : array-like, shape :math:`(n,)`, optional
            See :attr:`mol`.
        scale : array-like, shape :math:`(n,)`, optional
            See :attr:`scale`.
        validate : :class:`bool`
            If :data:`True` perform more thorough validation of the input arrays.
            :code:`validate = True` is recommend for general usage.
        copy : :class:`bool`
            If :data:`True`, always assign the passed arrays as copies.
            Only relevant if :code:`validate = True`.

        """
        array_dict = {'atoms': atoms, 'bonds': bonds, 'mol': mol, 'scale': scale}
        self._init(array_dict, validate, copy)

    def _init(self, array_dict, validate=True, copy=True):
        """Helper method for :meth:`HDF5Mol.__init__`."""
        if not validate:
            for name, ar in array_dict.items():
                try:
                    setattr(self, f'_{name}', ar)
                except AttributeError:
                    setattr(self, name, ar)
            return

        size = max((0 if ar is None else _alen(ar)) for _, ar in array_dict.items())
        for name, _ar in array_dict.items():
            func = self._parse_scale if name == 'scale' else self._parse_array
            ar = func(name, _ar, size=size, copy=copy)
            try:
                setattr(self, f'_{name}', ar)
            except AttributeError:
                setattr(self, name, ar)

    @classmethod
    def _parse_array(cls, name, value, size, copy=False):
        """Helper method for setting the :attr:`HDF5Mol.scale` property."""
        _, dtype, vlen, ndim, order = cls.INFO[name]
        return _to_array(value, size, dtype, vlen=vlen, copy=copy, ndim=ndim, order=order)

    @classmethod
    def _parse_scale(cls, name, value, size, copy=False):
        """Helper method for setting properties."""
        _, dtype, *_, order = cls.INFO[name]

        if value is None:
            return np.arange(size, dtype=dtype)

        ret = np.array(value, copy=copy, order=order)
        if not ret.ndim:  # 0D array; treat it as a scalar
            return np.full(size, ret, dtype=ret.dtype, order=order)
        else:
            if ret.shape != (size,):
                raise ValueError(f'The passed array should be of length {size}')
            return ret

    def __repr__(self):
        """Implement :class:`str(self)<str>` and :func:`repr(self)<repr>`."""
        wdith = max(len(k) for k in self.keys())

        def _str(k, v):
            dtype = v.dtype
            string_dtype = h5py.check_string_dtype(dtype)
            vlen_dtype = h5py.check_vlen_dtype(dtype)

            if string_dtype is not None:
                dtype_str = 'h5py.string_dtype(...)'
            elif vlen_dtype is not None:
                dt = str(vlen_dtype) if vlen_dtype.fields is None else '...'
                dtype_str = f'h5py.vlen_dtype({dt})'
            else:
                dtype_str = str(dtype) if dtype.fields is None else '...'

            return (f'{k:{wdith}} = {v.__class__.__module__}.{v.__class__.__name__}'
                    f'(..., shape={v.shape}, dtype={dtype_str})')

        ret = ',\n'.join(_str(k, v) for k, v in self.items())
        indent = 4 * ' '
        return f'{self.__class__.__name__}(\n{textwrap.indent(ret, indent)}\n)'

    def __reduce__(self):
        """Helper for :mod:`pickle`."""
        cls = type(self)
        return cls, (False, False, *self.values())

    def copy(self):
        """Return a deep copy of this instance.

        Returns
        -------
        :class:`HDF5Mol`
            A new HDF5Mol instance.

        """
        kwargs = {k: ar.copy() for k, ar in self.items()}
        cls = type(self)
        return cls(**kwargs, validate=False)  # type: ignore

    def __copy__(self):
        """Implement :func:`copy.copy(self)<copy.copy>`."""
        return self.copy()

    def __deepcopy__(self, memo=None):
        """Implement :func:`copy.deepcopy(self, memo=memo)<copy.deepcopy>`."""
        return self.copy()

    def __len__(self):
        """Implement :func:`len(self)<len>`.

        Returns
        -------
        :class:`int`
            The length of the arrays embedded within this instance
            (which are all of the same length).

        """
        return len(self.atoms)

    @staticmethod
    def _eq(vlen, ar1, ar2):
        """Helper method for :meth:`HDF5Mol.__eq__`."""
        if vlen:
            return all(np.all(i == j) for i, j in zip_longest(ar1, ar2, fillvalue=_False))
        else:
            return np.all(ar1 == ar2)

    def __eq__(self, value):
        """Implement :meth:`self == value<object.__eq__>`."""
        if type(self) is not type(value):
            return False

        INFO = self.INFO
        iterator = ((INFO[k].vlen, v, getattr(value, k)) for k, v in self.items())
        return all(self._eq(vlen, ar1, ar2) for vlen, ar1, ar2 in iterator)

    def __getitem__(self, index):
        """Implement :meth:`self[index]<object.__getitem__>`.

        Constructs a new :class:`HDF5Mol` instance by slicing all arrays with **index**.
        Follows the standard NumPy broadcasting rules:
        if an integer or slice is passed then the arrays are returned as views;
        they'll be copied otherwhise.

        Examples
        --------
        .. testsetup:: python

            >>> from scm.plams.testing_utils import HDF5_MOL as hdf5_mol

        .. code-block:: python

            >>> from scm.plams import HDF5Mol

            >>> hdf5_mol = HDF5Mol(...)  # doctest: +SKIP
            >>> print(hdf5_mol)
            HDF5Mol(
                atoms = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
                bonds = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
                mol   = numpy.ndarray(..., shape=(23,), dtype=...),
                scale = numpy.ndarray(..., shape=(23,), dtype=int64)
            )

            >>> hdf5_mol[0]
            HDF5Mol(
                atoms = numpy.ndarray(..., shape=(1,), dtype=h5py.vlen_dtype(...)),
                bonds = numpy.ndarray(..., shape=(1,), dtype=h5py.vlen_dtype(...)),
                mol   = numpy.ndarray(..., shape=(1,), dtype=...),
                scale = numpy.ndarray(..., shape=(1,), dtype=int64)
            )

            >>> hdf5_mol[:10]
            HDF5Mol(
                atoms = numpy.ndarray(..., shape=(10,), dtype=h5py.vlen_dtype(...)),
                bonds = numpy.ndarray(..., shape=(10,), dtype=h5py.vlen_dtype(...)),
                mol   = numpy.ndarray(..., shape=(10,), dtype=...),
                scale = numpy.ndarray(..., shape=(10,), dtype=int64)
            )

            >>> hdf5_mol[[0, 5, 7, 9, 10]]
            HDF5Mol(
                atoms = numpy.ndarray(..., shape=(5,), dtype=h5py.vlen_dtype(...)),
                bonds = numpy.ndarray(..., shape=(5,), dtype=h5py.vlen_dtype(...)),
                mol   = numpy.ndarray(..., shape=(5,), dtype=...),
                scale = numpy.ndarray(..., shape=(5,), dtype=int64)
            )

        Parameters
        ----------
        index : :class:`slice` or array-like
            An object for slicing the arrays in this instance.

        Returns
        -------
        :class:`HDF5Mol`
            A new HDF5Mol instance.

        """
        cls = type(self)
        idx = _parse_index(index, len(self))

        kwargs = {k: ar[idx] for k, ar in self.items()}
        return cls(**kwargs, validate=False)  # type: ignore

    @classmethod
    def keys(cls):
        """Yield the attribute names in this class.

        Examples
        --------
        .. code-block:: python

            >>> from scm.plams import HDF5Mol

            >>> for name in HDF5Mol.keys():
            ...     print(name)
            atoms
            bonds
            mol
            scale

        Yields
        ------
        :class:`str`
            The names of all attributes in this class.

        """
        return iter(cls.INFO.keys())

    def values(self):
        """Yield the attributes in this instance.

        Examples
        --------
        .. testsetup:: python

            >>> from scm.plams.testing_utils import HDF5_MOL as hdf5_mol

        .. code-block:: python

            >>> from scm.plams import HDF5Mol

            >>> hdf5_mol = HDF5Mol(...)  # doctest: +SKIP
            >>> for value in hdf5_mol.values():
            ...     print(object.__repr__(value))  # doctest: +ELLIPSIS
            <numpy.ndarray object at ...>
            <numpy.ndarray object at ...>
            <numpy.ndarray object at ...>
            <numpy.ndarray object at ...>

        Yields
        ------
        :class:`str`
            The values of all attributes in this instance.

        """
        return (getattr(self, k) for k in self.keys())

    def items(self):
        """Yield the attribute name/value pairs in this instance.

        Examples
        --------
        .. testsetup:: python

            >>> from scm.plams.testing_utils import HDF5_MOL as hdf5_mol

        .. code-block:: python

            >>> from scm.plams import HDF5Mol

            >>> hdf5_mol = HDF5Mol(...)  # doctest: +SKIP
            >>> for name, value in hdf5_mol.items():
            ...     print(name, '=', object.__repr__(value))  # doctest: +ELLIPSIS
            atoms = <numpy.ndarray object at ...>
            bonds = <numpy.ndarray object at ...>
            mol = <numpy.ndarray object at ...>
            scale = <numpy.ndarray object at ...>

        Yields
        ------
        :class:`str` and :class:`numpy.ndarray`
            The names and values of all attributes in this instance as 2-tuples.

        """
        return ((k, getattr(self, k)) for k in self.keys())

    def concatenate(self, *args):
        r"""Concatenate :math:`n` HDF5Mols into a single new instance.

        Examples
        --------
        .. testsetup:: python

            >>> from scm.plams.testing_utils import HDF5_MOL as pdb1

            >>> pdb2 = pdb3 = pdb1

        .. code-block:: python

            >>> from scm.plams import HDF5Mol

            >>> pdb1 = HDF5Mol(...)  # doctest: +SKIP
            >>> pdb2 = HDF5Mol(...)  # doctest: +SKIP
            >>> pdb3 = HDF5Mol(...)  # doctest: +SKIP
            >>> print(len(pdb1), len(pdb2), len(pdb3))
            23 23 23

            >>> pdb_new = pdb1.concatenate(pdb2, pdb3)
            >>> print(len(pdb_new))
            69

        Parameters
        ----------
        \*args : :class:`HDF5Mol`
            One or more HDF5Mols.

        Returns
        -------
        :class:`HDF5Mol`
            A new HDF5Mol cosntructed by concatenating **self** and **args**.

        """
        if not args:
            return self

        try:
            attr_list = [(k, v, [getattr(a, k) for a in args]) for k, v in self.items()]
        except AttributeError as ex:
            raise TypeError("'*args' expected one or more HDF5Mol instances") from ex

        cls = type(self)
        kwargs = {}
        for name, ar_self, ar_list in attr_list:

            # Ensure that the dtype of all passed scales is the same
            if name == 'scale':
                dtype = ar_self.dtype
                ar_list = [ar.astype(dtype, copy=False) for ar in ar_list]

            ar_new = np.concatenate((ar_self, *ar_list))
            kwargs[name] = ar_new
        return cls(**kwargs, validate=False)  # type: ignore

    @classmethod
    def from_molecules(cls, mol_list, scale=None, keys=None):
        """Convert a sequence of molecules into a new :class:`HDF5Mol` instance.

        Examples
        --------
        .. testsetup:: python

            >>> from scm.plams.testing_utils import (
            ...     MOL_TUPLE as mol_list
            ... )

        .. code-block:: python

            >>> from scm.plams import HDF5Mol, Molecule

            >>> mol_list = [Molecule(...), ...]  # doctest: +SKIP
            >>> HDF5Mol.from_molecules(mol_list)
            HDF5Mol(
                atoms = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
                bonds = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
                mol   = numpy.ndarray(..., shape=(23,), dtype=...),
                scale = numpy.ndarray(..., shape=(23,), dtype=int64)
            )

        Parameters
        ----------
        mol_list : :class:`Sequence[Molecule]<typing.Sequence>`
            A Sequence of PLAMS molecules.
        scale : array-like, optional
            An array-like object representing an user-specified index.
            Defaults to a simple range index if :data:`None` (see :func:`numpy.arange`).
        keys : ``"atoms"``, ``"bonds"`` and/or ``"mol"``, optional
            One or more strings with the names of the to-be updated arrays.
            For example, :code:`keys=["atoms"]` will only read the atoms in each molecule,
            :code:`keys=["atoms", "mol"]` will read both the molecule itself and its atoms,
            :code:`keys=["bonds"]` will only read the bonds, *etc.*.

        Returns
        -------
        :class:`HDF5Mol`
            A new HDF5Mol instance.

        """
        array_dict = {'scale': scale}
        return cls._from_molecules(mol_list, array_dict, keys)

    @classmethod
    def _from_molecules(cls, mol_list, array_dict, keys=None):
        """Helper function for :meth:`HDF5Mol.from_molecules`."""
        if keys is None:
            keys = cls.INFO.keys()

        # Parse the index
        kwargs = {}
        size = len(mol_list)
        for name, _ar in array_dict.items():
            if name == 'scale':
                kwargs[name] = cls._parse_scale(name, _ar, size=size)
            else:
                kwargs[name] = cls._parse_array(name, _ar, size=size)

        items = [('atoms', cls.get_atom_info),
                 ('bonds', cls.get_bond_info),
                 ('mol', cls.get_mol_info)]

        kwargs.update(
            {k: (cls._get_info(mol_list, k, v) if k in keys else None) for k, v in items}
        )
        validate = any(i is None for i in kwargs.values())
        return cls(**kwargs, validate=validate, copy=False)  # type: ignore

    @classmethod
    def _get_info(cls, mol_list, name, func):
        """Helper method for :meth:`HDF5Mol.from_molecules`.

        Used for managing the :meth:`get_atom_info`, :meth:`get_bond_info` and
        :meth:`get_mol_info` methods.
        """
        _, dtype, vlen, _ = cls.INFO[name]

        # np.fromiter() cannot be used for creating object arrays
        if vlen:
            vlen_dtype = h5py.vlen_dtype(dtype)
            array = np.empty_like(mol_list, dtype=vlen_dtype)
            array[:] = [func(mol, dtype) for mol in mol_list]
            return array
        else:
            count = len(mol_list)
            iterator = (func(mol, dtype) for mol in mol_list)
            return np.fromiter(iterator, dtype, count=count)

    @staticmethod
    def get_atom_info(mol, dtype):
        """Helper method for :meth:`HDF5Mol.from_molecules`.

        Get data from all atoms in **mol** at the specified **dtype**.

        Parameters
        ----------
        mol : :class:`~scm.plams.mol.molecule.Molecule`
            A PLAMS molecule.
        dtype : :class:`numpy.dtype`
            The data type of the to-be returned array.
            See :data:`HDF5Mol.INFO["atoms"].dtype<ATOMS_DTYPE>`.

        Returns
        -------
        :class:`numpy.ndarray`
            A 1-D array with the specified **dtype**.

        """
        count = len(mol.atoms)
        iterator = ((at.symbol, *at.coords) for at in mol)
        return np.fromiter(iterator, dtype, count=count)

    @staticmethod
    def get_bond_info(mol, dtype):
        """Helper method for :meth:`HDF5Mol.from_molecules`.

        Get data from all bonds in **mol** at the specified **dtype**.

        Parameters
        ----------
        mol : :class:`~scm.plams.mol.molecule.Molecule`
            A PLAMS molecule.
        dtype : :class:`numpy.dtype`
            The data type of the to-be returned array.
            See :data:`HDF5Mol.INFO["bonds"].dtype<BONDS_DTYPE>`.

        Returns
        -------
        :class:`numpy.ndarray`
            A 1-D array with the specified **dtype**.

        """
        mol.set_atoms_id()

        count = len(mol.bonds)
        iterator = ((b.atom1.id, b.atom2.id, b.order) for b in mol.bonds)
        ret = np.fromiter(iterator, dtype, count=count)

        mol.unset_atoms_id()
        return ret

    @staticmethod
    def get_mol_info(mol, dtype):
        """Helper method for :meth:`HDF5Mol.from_molecules`.

        Get data from **mol** at the specified **dtype**.

        Parameters
        ----------
        mol : :class:`~scm.plams.mol.molecule.Molecule`
            A PLAMS molecule.
        dtype : :class:`numpy.dtype`
            The data type of the to-be returned array
            See :data:`HDF5Mol.INFO["mol"].dtype<MOL_DTYPE>`.

        Returns
        -------
        :class:`numpy.ndarray`
            A 0-D array with the specified **dtype**.

        """
        data = mol.lattice
        return np.array(data if np.any(data) else 0, dtype=dtype)

    def to_molecules(self, mol=None, keys=None):
        """Create a molecule or list of molecules from this instance.

        Examples
        --------
        .. testsetup:: python

            >>> from scm.plams.testing_utils import (
            ...     HDF5_MOL as hdf5_mol,
            ...     MOL_TUPLE as mol_list,
            ...     MOL as mol
            ... )

        An example where one or more new molecules are created.

        .. code-block:: python

            >>> import numpy as np
            >>> from scm.plams import HDF5Mol, Molecule

            >>> hdf5_mol = HDF5Mol(...)  # doctest: +SKIP

            # Create three new molecules from `hdf5_mol`
            >>> hdf5_mol[0:3].to_molecules()  # doctest: +ELLIPSIS,+NORMALIZE_WHITESPACE
            [<scm.plams.mol.molecule.Molecule object at ...>,
             <scm.plams.mol.molecule.Molecule object at ...>,
             <scm.plams.mol.molecule.Molecule object at ...>]

        An example where one or more existing molecules are updated in-place.

        .. code-block:: python

            # Update all molecules in `mol_list` with info from `hdf5_mol`
            >>> mol_list = [Molecule(...), ...]  # doctest: +SKIP
            >>> mol_list_new = hdf5_mol[0:3].to_molecules(mol=mol_list)

            # Note that the molecules in `mol_list_new` are the same as in `mol_list`
            >>> iterator = zip(mol_list, mol_list_new)
            >>> all(m is m_new for m, m_new in iterator)
            True

        Parameters
        ----------
        mol : :class:`~scm.plams.mol.molecule.Molecule` or :class:`Sequence[Molecule]<typing.Sequence>`, optional
            A molecule or list of molecules.
            If not :data:`None` then the passed molecule(s) will be updated in-place.
        keys : ``"atoms"``, ``"bonds"`` and/or ``"mol"``, optional
            One or more strings with the names of the to-be exported arrays.
            For example, :code:`keys=["atoms"]` will only update the atoms in each molecule,
            :code:`keys=["atoms", "mol"]` will update both the molecule itself and its atoms,
            :code:`keys=["bonds"]` will only update the bonds, *etc.*.

        Returns
        -------
        :class:`List[Molecule]<typing.List>`
            A list of molecules.

        """  # noqa: E501
        if keys is None:
            keys = self.INFO.keys()

        if mol is None:
            mol_list: Iterable[Optional[Molecule]] = repeat(None)
        elif isinstance(mol, Molecule):
            mol_list = [mol]
        else:
            mol_list = mol

        mol_keys = {'atoms', 'bonds', 'mol'}
        values = (v for k, v in self.items() if k in mol_keys)
        iterator = zip(*values, mol_list)
        return [self._to_molecules(*args, keys, mol) for *args, mol in iterator]

    @classmethod
    def _to_molecules(cls, atoms, bonds, mol, keys, molecule=None):
        """Helper method for :meth:`HDF5Mol.to_molecules`: update/create a single molecule from the passed atoms and bonds."""  # noqa: E501
        # Parse the input molecule and check its atoms
        if molecule is None:
            not_atoms = True
            ret = Molecule()
        else:
            not_atoms = False
            ret = molecule

        # Check if the molecule has any bonds
        not_bonds = True if not ret.bonds else False

        # Update the molecule
        if 'mol' in keys:
            cls.set_mol_info(ret, mol)

        # Update atoms
        if 'atoms' in keys:
            if not_atoms:
                ret.atoms = [Atom(mol=ret) for _ in atoms]
            cls.set_atom_info(ret, atoms)

        # Update bonds
        if 'bonds' in keys:
            if not ret.bonds:
                ret.bonds = [Bond(mol=ret) for _ in bonds]
            else:
                for at in ret:
                    at.bonds = []
            cls.set_bond_info(ret, bonds)

        return ret

    @staticmethod
    def set_atom_info(mol, data):
        """Helper method for :meth:`HDF5Mol.to_molecules`.

        Update all atoms in **mol** with **data**.

        Parameters
        ----------
        mol : :class:`~scm.plams.mol.molecule.Molecule`
            A PLAMS molecule.
        data : :class:`numpy.ndarray`
            The to-be set data.
            See :data:`HDF5Mol.INFO["atoms"].dtype<ATOMS_DTYPE>` for the data type.

        """
        for atom, record in zip(mol.atoms, data):
            symbol, *coords = record.item()
            atom.symbol = symbol.decode()
            atom.coords = coords

    @staticmethod
    def set_bond_info(mol, data):
        """Helper method for :meth:`HDF5Mol.to_molecules`.

        Update all bonds in **mol** with **data**.

        Parameters
        ----------
        mol : :class:`~scm.plams.mol.molecule.Molecule`
            A PLAMS molecule.
        data : :class:`numpy.ndarray`
            The to-be set data.
            See :data:`HDF5Mol.INFO["bonds"].dtype<BONDS_DTYPE>` for the data type.

        """
        for bond, record in zip(mol.bonds, data):
            idx1, idx2, order = record.item()
            atom1 = mol[idx1]
            atom2 = mol[idx2]

            bond.order = order
            bond.atom1 = atom1
            bond.atom2 = atom2

            atom1.bonds.append(bond)
            atom2.bonds.append(bond)

    @staticmethod
    def set_mol_info(mol, data):
        """Helper method for :meth:`HDF5Mol.to_molecules`.

        Update **mol** with **data**.

        Parameters
        ----------
        mol : :class:`~scm.plams.mol.molecule.Molecule`
            A PLAMS molecule.
        data : :class:`numpy.void`
            The to-be set data.
            See :data:`HDF5Mol.INFO["mol"].dtype<MOL_DTYPE>` for the data type.

        """
        lattice, *_ = data.item()
        mol.lattice = lattice.tolist()

    def create_group(self, file, name,
                     atoms=None, bonds=None, mol=None, scale=None, create_group=True, **kwargs):
        r"""Create a h5py Group for storing :class:`HDF5Mol` instances.

        Examples
        --------
        .. testsetup:: python

            >>> import os
            >>> from scm.plams.testing_utils import (
            ...     HDF5_TMP as file,
            ...     HDF5_MOL as hdf5_mol
            ... )

            >>> if os.path.isfile(file):
            ...     os.remove(file)

        .. code-block:: python

            >>> import h5py
            >>> from scm.plams import HDF5Mol

            >>> file = str(...)  # doctest: +SKIP
            >>> hdf5_mol = HDF5Mol(...)  # doctest: +SKIP

            >>> with h5py.File(file, 'a') as f:
            ...     hdf5_mol.create_group(f, 'ligand')
            <HDF5 group "/ligand" (4 members)>

        .. testcleanup:: python

            >>> if os.path.isfile(file):
            ...     os.remove(file)

        Parameters
        ----------
        file : :class:`h5py.File` or :class:`h5py.Group`
            The h5py File or Group where the new Group will be created.
        name : :class:`str`
            The name of the to-be created Group.
        atoms : :class:`h5py.Dataset`, optional
            An pre-existing dataset for containg :attr:`HDF5Mol.atoms`.
            A soft link will be created if not :data:`None`.
        bonds : :class:`h5py.Dataset`, optional
            An pre-existing dataset for containg :attr:`HDF5Mol.bonds`.
            A soft link will be created if not :data:`None`.
        mol : :class:`h5py.Dataset`, optional
            An pre-existing dataset for containg :attr:`HDF5Mol.mol`.
            A soft link will be created if not :data:`None`.
        scale : :class:`h5py.Dataset`, optional
            An pre-existing dataset for containg :attr:`HDF5Mol.scale`.
            A soft link will be created if not :data:`None`.
        create_group : :class:`bool`
            If :data:`False`, create the datasets directly in the passed **file**,
            rather than first creating a group.
        \**kwargs : :data:`~typing.Any`
            Further keyword arguments for the creation of each dataset.
            Arguments already specified by default are
            ``name`` and ``dtype``; these cannot be overridden.
            See :meth:`h5py.Group.create_dataset` for more details.

        Returns
        -------
        :class:`h5py.Group`
            The newly created Group.

        """
        dset_dict = {'atoms': atoms, 'bonds': bonds, 'mol': mol, 'scale': scale}
        return self._create_group(file, name, dset_dict, create_group, **kwargs)

    def _create_group(self, file, name, dset_dict, create_group=True, **kwargs):
        """Helper method for :meth:`HDF5Mol.create_group`."""
        cls = type(self)
        if create_group:
            grp = file.create_group(name, track_order=True)
            grp.attrs['type'] = f'{cls.__module__}.{cls.__name__}'
        else:
            grp = file

        # Create the datasets
        for key, value in dset_dict.items():
            dset_name, dtype, vlen, ndim = cls.INFO[key]
            if key == 'scale':
                dtype = self.scale.dtype

            if vlen:
                dtype = h5py.vlen_dtype(dtype)

            if value is not None:
                # Create a softlink or softlink+ (if they originate from different files)
                if value.file == grp.file:
                    grp[dset_name] = h5py.SoftLink(value.name)
                else:
                    grp[dset_name] = h5py.ExternalLink(value.file.filename, value.name)

                # Check if the datatypes are compatible
                dset = value
                if dset.dtype != dtype:
                    msg = f"dtype mistmatch between {dset} and {cls.__name__}.{key}"
                    warnings.warn(msg, RuntimeWarning, stacklevel=2)

            else:
                kwargs.setdefault('shape', (0,) + (1,) * (ndim - 1))
                kwargs.setdefault('maxshape', (None,) * ndim)
                dset = grp.create_dataset(dset_name, dtype=dtype, **kwargs)

        # Set the scale
        try:
            scale_name = cls.INFO['scale'].name
        except KeyError:
            return
        scale_dset = grp[scale_name]
        scale_dset.make_scale(scale_name)

        # Use the index as a scale
        dset_iter = (grp[cls.INFO[k].name] for k in cls.keys() if k != 'scale')
        for dset in dset_iter:
            dset.dims[0].label = scale_name
            dset.dims[0].attach_scale(scale_dset)
        return grp

    def to_hdf5(self, group, index=None, keys=None):
        """Update all datasets in **group** positioned at **index** with its counterpart from **pdb**.

        Follows the standard broadcasting rules as employed by h5py.

        Examples
        --------
        .. testsetup:: python

            >>> import os
            >>> from scm.plams.testing_utils import (
            ...     HDF5_TMP as file,
            ...     HDF5_MOL as hdf5_mol
            ... )

            >>> if os.path.isfile(file):
            ...     os.remove(file)

        .. code-block:: python

            >>> import h5py
            >>> import numpy as np
            >>> from scm.plams import HDF5Mol

            >>> file = str(...)  # doctest: +SKIP
            >>> hdf5_mol = HDF5Mol(...)  # doctest: +SKIP
            >>> print(hdf5_mol)
            HDF5Mol(
                atoms = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
                bonds = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
                mol   = numpy.ndarray(..., shape=(23,), dtype=...),
                scale = numpy.ndarray(..., shape=(23,), dtype=int64)
            )

            >>> with h5py.File(file, 'a') as f:
            ...     dtype = hdf5_mol.scale.dtype
            ...     group = hdf5_mol.create_group(f, 'ligand')
            ...     hdf5_mol.to_hdf5(group, index=np.s_[0:23])
            ...
            ...     for dataset in group.values():
            ...         print(dataset)
            <HDF5 dataset "atoms": shape (23,), type "|O">
            <HDF5 dataset "bonds": shape (23,), type "|O">
            <HDF5 dataset "mol": shape (23,), type "|V72">
            <HDF5 dataset "scale": shape (23,), type "<i8">

        .. testcleanup:: python

            >>> if os.path.isfile(file):
            ...     os.remove(file)

        Important
        ---------
        If **index** is passed as a sequence of integers then, contrary to NumPy,
        they *will* have to be sorted.

        Parameters
        ----------
        group : :class:`h5py.Group`
            The to-be updated h5py group.
        index : :class:`slice` or array-like, optional
            An object for slicing all datasets in **group**.
            Note that, contrary to numpy, if a sequence of integers is provided
            then they *will* have to ordered.
        keys : ``"atoms"``, ``"bonds"``, ``"mol"`` and/or ``"scale"``, optional
            One or more strings with the names of the to-be exported arrays.
            For example, :code:`keys=["atoms"]` will only export the :attr:`~HDF5Mol.atoms` array,
            :code:`keys=["atoms", "mol"]` will export both the :attr:`~HDF5Mol.atoms` and
            :attr:`~HDF5Mol.mol` arrays, *etc.*.

        """
        if keys is None:
            keys = self.INFO.keys()

        try:
            idx = idx_max = index.__index__()  # type: ignore
        except (TypeError, AttributeError):
            idx = _parse_index(index, None)
            if isinstance(idx, slice):
                idx_max = idx.stop or len(self)
            else:
                generic = idx.dtype.type
                if issubclass(generic, np.integer):
                    idx_max = idx.max() + 1
                elif issubclass(generic, np.bool_):
                    idx_max = len(self)
                else:
                    raise TypeError("'index' expected a boolean or integer array; "
                                    f"observed type: {idx_max.__class__.__name__}[{generic.__name__}]")

        # Update the datasets
        self._update_hdf5_shape(group, idx_max)
        items = ((group[self.INFO[k].name], ar) for k, ar in self.items() if k in keys)
        for dataset, ar in items:
            dataset[idx] = ar

    def _update_hdf5_shape(self, group, min_len):
        """Update the shape of all datasets in **group** such their length is at least as long as **min_len**.

        The length of all datasets will be set equal to the ``index`` dimensional cale.

        This method is automatically called by :meth:`HDF5Mol.update_hdf5`.

        Parameters
        ----------
        group : :class:`h5py.Group`
            The h5py Group with the to-be reshaped datasets.
        min_len : :class:`int`
            The minimum desired length of each dataset.

        """  # noqa: E501
        dset_iter = ((k, group[self.INFO[k].name]) for k in self.keys())
        for name, dataset in dset_iter:
            shape = (min_len, *getattr(self, name).shape[1:])
            for i, j in enumerate(shape):
                if dataset.shape[i] < j:
                    dataset.resize(j, axis=0)

    @classmethod
    def from_hdf5(cls, group, index = None, keys=None):
        """Construct a new HDF5Mol from the passed hdf5 **group**.

        Examples
        --------
        .. testsetup:: python

            >>> from scm.plams.testing_utils import HDF5_READ as file

        .. code-block:: python

            >>> import h5py
            >>> import numpy as np
            >>> from scm.plams import HDF5Mol

            >>> file = str(...)  # doctest: +SKIP

            >>> with h5py.File(file, 'r') as f:
            ...     group = f['ligand']
            ...     HDF5Mol.from_hdf5(group)
            HDF5Mol(
                atoms = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
                bonds = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
                mol   = numpy.ndarray(..., shape=(23,), dtype=...),
                scale = numpy.ndarray(..., shape=(23,), dtype=int64)
            )

        Parameters
        ----------
        group : :class:`h5py.Group`
            The to-be read h5py group.
        index : :class:`slice` or array-like, optional
            An object for slicing all datasets in **group**.
        keys : ``"atoms"``, ``"bonds"``, ``"mol"`` and/or ``"scale"``, optional
            One or more strings with the names of the to-be updated arrays.
            For example, :code:`keys=["atoms"]` will only update the :attr:`~HDF5Mol.atoms` array,
            :code:`keys=["atoms", "mol"]` will update both the :attr:`~HDF5Mol.atoms` and
            :attr:`~HDF5Mol.mol` arrays, *etc.*.

        Returns
        -------
        :class:`HDF5Mol`
            A new HDF5Mol instance constructed from **group**.

        """
        if keys is None:
            keys = cls.INFO.keys()

        try:
            idx = index.__index__()  # type: ignore
        except (AttributeError, TypeError):
            idx = _parse_index(index, None)

        def _get_array(name):
            return group[cls.INFO[name].name][idx]

        kwargs = {k: (_get_array(k) if k in keys else None) for k in cls.keys()}
        validate = any(i is None for i in kwargs.values())
        return cls(**kwargs, validate=validate, copy=False)  # type: ignore

    def intersection(self, value):
        """Construct a new HDF5Mol by the intersection of **self** and **value**.

        Examples
        --------
        .. testsetup:: python

            >>> from scm.plams.testing_utils import HDF5_MOL as hdf5_mol

            >>> hdf5_mol2 = range(4)

        An example where one or more new molecules are created.

        .. code-block:: python

            >>> from scm.plams import HDF5Mol

            >>> hdf5_mol = HDF5Mol(..., range(23))  # doctest: +SKIP
            >>> hdf5_mol2 = HDF5Mol(..., scale=range(4))  # doctest: +SKIP

            >>> print(hdf5_mol.scale)
            [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22]

            >>> hdf5_mol_new = hdf5_mol.intersection(hdf5_mol2)
            >>> print(hdf5_mol_new.scale)
            [0 1 2 3]

        Parameters
        ----------
        value : :class:`HDF5Mol` or array-like
            Another HDF5Mol or an array-like object representing :attr:`HDF5Mol.scale`.
            Note that both **value** and **self.scale** should consist of unique elements.

        Returns
        -------
        :class:`HDF5Mol`
            A new instance by intersecting :attr:`self.scale<HDF5Mol.scale>` and **value**.

        See Also
        --------
        :meth:`set.intersection<frozenset.intersection>`
            Return the intersection of two sets as a new set.

        """
        idx1, idx2 = self._get_index(value)
        _, i, _ = np.intersect1d(idx1, idx2, assume_unique=True, return_indices=True)
        return self[i]

    def difference(self, value):
        """Construct a new HDF5Mol by the difference of **self** and **value**.

        Examples
        --------
        .. testsetup:: python

            >>> from scm.plams.testing_utils import HDF5_MOL as hdf5_mol

            >>> hdf5_mol2 = range(10, 30)

        .. code-block:: python

            >>> from scm.plams import HDF5Mol

            >>> hdf5_mol = HDF5Mol(..., scale=range(23))  # doctest: +SKIP
            >>> hdf5_mol2 = HDF5Mol(..., scale=range(10, 30))  # doctest: +SKIP

            >>> print(hdf5_mol.scale)
            [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22]

            >>> hdf5_mol_new = hdf5_mol.difference(hdf5_mol2)
            >>> print(hdf5_mol_new.scale)
            [0 1 2 3 4 5 6 7 8 9]

        Parameters
        ----------
        value : :class:`HDF5Mol` or array-like
            Another HDF5Mol or an array-like object representing :attr:`HDF5Mol.scale`.
            Note that both **value** and **self.scale** should consist of unique elements.

        Returns
        -------
        :class:`HDF5Mol`
            A new instance as the difference of :attr:`self.scale<HDF5Mol.scale>`
            and **value**.

        See Also
        --------
        :meth:`set.difference<frozenset.difference>`
            Return the difference of two or more sets as a new set.

        """
        idx1, idx2 = self._get_index(value)
        i = np.in1d(idx1, idx2, assume_unique=True, invert=True)
        return self[i]

    def symmetric_difference(self, value):
        """Construct a new HDF5Mol by the symmetric difference of **self** and **value**.

        Examples
        --------
        .. testsetup:: python

            >>> import numpy as np
            >>> from scm.plams import HDF5Mol
            >>> from scm.plams.testing_utils import HDF5_MOL as hdf5_mol

            >>> a = np.random.rand(20)
            >>> hdf5_mol2 = HDF5Mol(atoms=a, bonds=a, mol=a, scale=np.arange(10, 30))

        .. code-block:: python

            >>> from scm.plams import HDF5Mol

            >>> hdf5_mol = HDF5Mol(..., range(23))  # doctest: +SKIP
            >>> hdf5_mol2 = HDF5Mol(..., scale=range(10, 30))  # doctest: +SKIP

            >>> print(hdf5_mol.scale)
            [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22]

            >>> hdf5_mol_new = hdf5_mol.symmetric_difference(hdf5_mol2)
            >>> print(hdf5_mol_new.scale)
            [ 0  1  2  3  4  5  6  7  8  9 23 24 25 26 27 28 29]

        Parameters
        ----------
        value : :class:`HDF5Mol`
            Another HDF5Mol.
            Note that both **value.scale** and **self.scale** should consist of unique elements.

        Returns
        -------
        :class:`HDF5Mol`
            A new instance as the symmetric difference of :attr:`self.scale<HDF5Mol.scale>`
            and **value**.

        See Also
        --------
        :meth:`set.symmetric_difference<frozenset.symmetric_difference>`
            Return the symmetric difference of two sets as a new set.

        """
        pdb1 = self.difference(value)
        pdb2 = value.difference(self)
        return pdb1.concatenate(pdb2)

    def union(self, value):
        """Construct a new HDF5Mol by the union of **self** and **value**.

        Examples
        --------
        .. testsetup:: python

            >>> import numpy as np
            >>> from scm.plams import HDF5Mol
            >>> from scm.plams.testing_utils import HDF5_MOL as hdf5_mol

            >>> hdf5_mol2 = HDF5Mol(scale=np.arange(10, 30))

        .. code-block:: python

            >>> from scm.plams import HDF5Mol

            >>> hdf5_mol = HDF5Mol(..., range(23))  # doctest: +SKIP
            >>> hdf5_mol2 = HDF5Mol(..., scale=range(10, 30))  # doctest: +SKIP

            >>> print(hdf5_mol.scale)
            [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22]

            >>> hdf5_mol_new = hdf5_mol.union(hdf5_mol2)
            >>> print(hdf5_mol_new.scale)
            [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
             24 25 26 27 28 29]

        Parameters
        ----------
        value : :class:`HDF5Mol`
            Another HDF5Mol.
            Note that both **value** and **self.scale** should consist of unique elements.

        Returns
        -------
        :class:`HDF5Mol`
            A new instance as the union of :attr:`self.scale<HDF5Mol.scale>` and **value**.

        See Also
        --------
        :meth:`set.union<frozenset.union>`
            Return the union of sets as a new set.

        """
        cls = type(self)
        if not isinstance(value, cls):
            raise TypeError(f"'value' expected a {cls.__name__} instance; "
                            f"observed type: {value.__class__.__name__}")

        idx1 = self.scale
        idx2 = value.scale.astype(idx1.dtype, copy=False)
        i = np.in1d(idx2, idx1, assume_unique=True, invert=True)

        ret = value[i]
        return self.concatenate(ret)

    def _get_index(self, value):
        """Parse and return the :attr:`~HDF5Mol.scale` of **self** and **value**."""
        scale1 = self.scale
        scale2 = np.asarray(getattr(value, 'scale', value), dtype=scale1.dtype)
        return scale1, scale2
