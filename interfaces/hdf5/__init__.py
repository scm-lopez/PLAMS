"""Test.

.. currentmodule:: scm.plams

The :class:`HDF5Mol` class serves as an intermediate container between the PLAMS
:class:`~scm.plams.mol.molecule.Molecule` class and the h5py_ interface to the hdf5 format.

.. _h5py: https://www.h5py.org/

.. SeeAlso:: h5py_

    .. image:: https://badge.fury.io/py/h5py.svg
        :target: https://badge.fury.io/py/h5py

    The h5py package provides both a high- and low-level interface to the HDF5 library from Python.
    The low-level interface is intended to be a complete wrapping of the HDF5 API,
    while the high-level component supports access to HDF5 files,
    datasets and groups using established Python and NumPy concepts.

    A strong emphasis on automatic conversion between Python (Numpy) datatypes and
    data structures and their HDF5 equivalents vastly simplifies the process of reading and
    writing data from Python.


Basic Usage
-----------
An example which shows some basic hdf5/Molecule interconversion:

1.  Define the initial molecules of interest.
2.  Create a h5py group and export the molecules to aforemointed group.
3.  Import molecules from the group.

.. testsetup:: python

    >>> import os
    >>> from scm.plams.testing_utils import (
    ...     MOL_TUPLE as mol_list,
    ...     HDF5_TMP as filename,
    ...     assert_mol_eq
    ... )

    >>> if os.path.isfile(filename):
    ...     os.remove(filename)

.. code-block:: python
    :emphasize-lines: 4, 8, 14

    >>> import h5py
    >>> from scm.plams import Molecule, HDF5Mol

    >>> # 1: Define the input molecules
    >>> mol_list = [Molecule(...), ...]  # doctest: +SKIP
    >>> filename = str(...)  # doctest: +SKIP

    >>> # 2: Export molecules to the hdf5 format
    >>> hdf5_mol = HDF5Mol.from_molecules(mol_list)
    >>> with h5py.File(filename, 'a') as f:
    ...     group = hdf5_mol.create_group(f, "ligand")
    ...     hdf5_mol.to_hdf5(group)

    >>> # 3: Import molecules to the hdf5 format
    >>> with h5py.File(filename, 'r') as f:
    ...     group = f["ligand"]
    ...     hdf5_mol2 = HDF5Mol.from_hdf5(group)
    ...     mol_list2 = hdf5_mol2.to_molecules()

.. testcleanup::

    >>> if os.path.isfile(filename):
    ...     os.remove(filename)

    >>> for m_new, m in zip(mol_list2, mol_list):
    ...     assert_mol_eq(m_new, m)

Class Structure
---------------
At its basics the :class:`HDF5Mol` class consists of four numpy arrays which are
responsible for containg all atom, bond and molecule-related information.
The fourth array, ``"scale"``, is used for holding user-specified identifiers for
each molecule and will be explained in more detail in `Dimensional Scales`_.

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

In order to hold the varied data contained within a molecule these arrays
are make use of so called structured data types.
The matching structured arrays are divided into one or more named fields,
each field containing its own data type and, in a sense,
they can be thought of as a dictionary of arrays.
This is illustrated in the example below, where a structured array is created
for holding atomic symbols and the Cartesian coordinates.

See `Structured arrays <https://numpy.org/doc/stable/user/basics.rec.html>`_
for more details.

.. code-block:: python

    >>> import numpy as np

    >>> dtype = [("symbol", "S2"),
    ...          ("x", float),
    ...          ("y", float),
    ...          ("z", float)]

    >>> array = np.array(("H", 0, 1, 2), dtype=dtype, ndmin=1)
    >>> array
    array([(b'H', 0., 1., 2.)],
          dtype=[('symbol', 'S2'), ('x', '<f8'), ('y', '<f8'), ('z', '<f8')])

    >>> print(array["symbol"])
    [b'H']

    >>> print(array["x"])
    [0.]

See `Arbitrary vlen data <https://h5py.readthedocs.io/en/stable/special.html#arbitrary-vlen-data>`_
for more details.

.. code-block:: python

    >>> with np.printoptions(threshold=4, edgeitems=2):
    ...     print(hdf5_mol.atoms[:2])  # doctest: +NORMALIZE_WHITESPACE
    [array([(b'C', 3.861, -0.415,  0.046),
            (b'C', 2.506,  0.297, -0.033),
            ...,
            (b'H', 3.936, -0.956,  1.015),
            (b'H', 3.936, -1.156, -0.779)],
           dtype=[('symbol', 'S2'), ('x', '<f8'), ('y', '<f8'), ('z', '<f8')])
     array([(b'C', 3.606,  0.203,  0.576),
            (b'O', 2.237,  0.297,  0.822),
            ...,
            (b'H', 3.749, -0.162, -0.461),
            (b'H', 4.095,  1.194,  0.695)],
           dtype=[('symbol', 'S2'), ('x', '<f8'), ('y', '<f8'), ('z', '<f8')])]

Dimensional Scales
------------------
The example below serves as an illustration.
Herein a list of SMILES strings is used both as dimensional scale and
to construct the molecules themselves.

Dimensional scales in hdf5 can be thought of as arrays of meta-data,
providing additional context to the elements of one or more other arrays.
In this sence they share similarities with dictionary keys,
though the link between keys and values is implicit.

See `Dimension Scales <https://docs.h5py.org/en/latest/high/dims.html>`_ for more details.

.. testsetup:: python

    >>> import os
    >>> from scm.plams.testing_utils import (
    ...     SMILES_TUPLE as smiles_list,
    ...     HDF5_TMP as filename
    ... )

    >>> if os.path.isfile(filename):
    ...     os.remove(filename)

.. code-block:: python

    >>> import h5py
    >>> import numpy as np
    >>> from scm.plams import Molecule, HDF5Mol, from_smiles

    >>> # Create a list of 23 SMILES strings and correpsonding molecules
    >>> smiles_list = ["CN", "CO", "CCO", ...]  # doctest: +SKIP
    >>> mol_list = [from_smiles(smiles) for smiles in smiles_list]

    >>> # Use (variable-length) smiles strings as scale
    >>> dtype = h5py.string_dtype('ascii')
    >>> scale = np.array(smiles_list, dtype=bytes).astype(dtype)

.. testcleanup:: python

    >>> if os.path.isfile(filename):
    ...     os.remove(filename)

The SMILES strings are stored using the h5py variable length data type
(:func:`h5py.string_dtype`), which grants additional (much needed) flexibility
over the alternative: fixed-length byte strings (:class:`numpy.string_`).

See `Strings in HDF5 <https://h5py.readthedocs.io/en/stable/strings.html>`_ for more details.

After the molecules and (optionally) the scale are defined
the new :class:`HDF5Mol` instance can be constructed using
the :meth:`HDF5Mol.from_molecules` method.
The resulting HDF5Mol instance is, in essence, a fancy container for holding four numpy arrays,
the latter being reflected in its attributes.

.. code-block:: python

    >>> hdf5_mol = HDF5Mol.from_molecules(mol_list, scale=scale)
    >>> print(hdf5_mol)
    HDF5Mol(
        atoms = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
        bonds = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
        mol   = numpy.ndarray(..., shape=(23,), dtype=...),
        scale = numpy.ndarray(..., shape=(23,), dtype=h5py.string_dtype(...))
    )

HDF5 Groups
-----------
With the new HDF5Mol instance present a h5py Group can be created,
which will hold the Datasets for storing the molecules.
Groups and Datasets are the two central objects behind the h5py interface,
respectively behaving very much like (nested) dictionaries and numpy arrays.

See :class:`h5py.Group` and :class:`h5py.Dataset` for more details.

.. code-block:: python

    >>> filename = str(...)  # doctest: +SKIP

    >>> # Create a group by the name of "ligand"
    >>> with h5py.File(filename, 'a') as f:
    ...     group = hdf5_mol.create_group(f, "ligand")
    ...
    ...     print(group)
    ...     print(list(group.values()))  # doctest: +NORMALIZE_WHITESPACE
    <HDF5 group "/ligand" (4 members)>
    [<HDF5 dataset "atoms": shape (0,), type "|O">,
     <HDF5 dataset "bonds": shape (0,), type "|O">,
     <HDF5 dataset "mol": shape (0,), type "|V72">,
     <HDF5 dataset "scale": shape (0,), type "|O">]

Additional options for the group construction are available
in :meth:`HDF5Mol.create_group`.
For example, one might want to use one or more pre-existing datasets
instead of reconstructing them from scratch.
This is possible via the ``"atoms"``, ``"bonds"``, *etc.* keywords,
which will create a soft link to whichever datasets are passed.

See `Dict interface and links <https://h5py.readthedocs.io/en/stable/high/group.html#dict-interface-and-links>`_
for more details.

By the default each HDF5Mol attribute is stored under an identically-named dataset.
This behavior is configurable via the named tuples stored :class:`HDF5Mol.INFO[...]<InfoTuple>`,
which contain, among others, the dataset names.
This is illustrated in the example below,
where the names of the :attr:`~HDF5Mol.atoms` and :attr:`~HDF5Mol.mol` datasets
are respectivelly changed into ``"atoms_dset"`` and ``"path/to/mol"``.
The second name change deserves special attention,
as the way h5py creates datasets and groups share many similarties with file paths.
Namelly, ``"path/to/mol"`` is automatically expanded to groups and datasets
(:code:`f["path"]["to"]["mol"] = ...` in dictionary notation).

See `Creating groups <https://h5py.readthedocs.io/en/stable/high/group.html#creating-groups>`_ for more details
(note that this applies to both groups and datasets.)

.. code-block:: python

    >>> class SubClass(HDF5Mol):
    ...     INFO = HDF5Mol.INFO.copy()
    ...     INFO["atoms"] = INFO["atoms"]._replace(name="atoms_dset")
    ...     INFO["mol"] = INFO["mol"]._replace(name="path/to/mol")

    >>> sub_cls = SubClass()

    >>> with h5py.File(filename, 'a') as f:
    ...     group = sub_cls.create_group(f, "sub_cls")
    ...
    ...     print(group["atoms_dset"])
    ...     print(group["path"]["to"]["mol"])  # Equivalent to group["path/to/mol"]
    <HDF5 dataset "atoms_dset": shape (0,), type "|O">
    <HDF5 dataset "mol": shape (0,), type "|V72">

Molecule Importing/Exporting
----------------------------
Finally the molecules can be exported to the newly created group.
The relevant datasets are automatically resized such that accomodate
the indices as specified in the ``index`` parameter of :meth:`HDF5Mol.to_hdf5`.
In the example below this is illustrated with the :code:`[0:23]` slice.

.. code-block:: python

    >>> with h5py.File(filename, 'a') as f:
    ...     group = f['ligand']
    ...
    ...     # Define which elements in the group's datasets should be updated
    ...     index = np.s_[0:23]
    ...     hdf5_mol.to_hdf5(group, index=index)
    ...
    ...     print(list(group.values()))  # doctest: +NORMALIZE_WHITESPACE
    [<HDF5 dataset "atoms": shape (23,), type "|O">,
     <HDF5 dataset "bonds": shape (23,), type "|O">,
     <HDF5 dataset "mol": shape (23,), type "|V72">,
     <HDF5 dataset "scale": shape (23,), type "|O">]

Converselly molecules can be imported, rather than exported, from the group.
The relevant :meth:`HDF5Mol.to_molecules` method supports both the creation of
new molecules as well as the (in-place) updating of pre-existing molecules via
the ``mol`` argument.

.. testsetup:: python

    >>> import builtins
    >>> import reprlib

    >>> repr_ = reprlib.Repr()
    >>> repr_.maxother = 1000
    >>> def print(obj):
    ...     builtins.print(repr_.repr(obj))

.. code-block:: python

    >>> with h5py.File(filename, 'r+') as f:
    ...     group = f['ligand']
    ...
    ...     index = np.s_[0:23]
    ...     hdf5_mol2 = HDF5Mol.from_hdf5(group, index=index)
    ...     mol_list2 = hdf5_mol2.to_molecules()
    ...
    ...     print(mol_list2)  # doctest: +ELLIPSIS,+NORMALIZE_WHITESPACE
    [<scm.plams.mol.molecule.Molecule object at ...>,
     <scm.plams.mol.molecule.Molecule object at ...>,
     <scm.plams.mol.molecule.Molecule object at ...>,
     ...]

.. testcleanup:: python

    >>> print = builtins.print

If one wants to, for example, only use the ``HDF5Mol`` for
constructing the molecule's atoms, thus ignoring any bonds,
then this can be accomplished with the ``keys`` parameter.
Accepted values are ``"atoms"``, ``"bonds"`` and ``"mol"``,
named after the respective :class:`HDF5Mol` attributes.
Note that the ``keys`` parameter is available to all functions
used for Molecule/hdf5 interconversion.

.. code-block:: python

    >>> # Export everything
    >>> mol_list1 = hdf5_mol.to_molecules()
    >>> print(mol_list1[5])  # doctest: +NORMALIZE_WHITESPACE
      Atoms:
        1         C     -0.063899      0.013017     -0.012822
        2         O      1.326996     -0.105806      0.098365
        3         H     -0.290404      0.952594     -0.555502
        4         H     -0.456427     -0.888371     -0.529138
        5         H     -0.516266      0.028566      0.999097
      Bonds:
       (1)--1.0--(2)
       (1)--1.0--(3)
       (1)--1.0--(4)
       (1)--1.0--(5)
      Lattice:
            0.0000000000     0.0000000000     0.0000000000
            0.0000000000     0.0000000000     0.0000000000
            0.0000000000     0.0000000000     0.0000000000
    <BLANKLINE>

    # Export just the atoms
    >>> mol_list2 = hdf5_mol.to_molecules(keys=["atoms"])
    >>> print(mol_list2[5])  # doctest: +NORMALIZE_WHITESPACE
      Atoms:
        1         C     -0.063899      0.013017     -0.012822
        2         O      1.326996     -0.105806      0.098365
        3         H     -0.290404      0.952594     -0.555502
        4         H     -0.456427     -0.888371     -0.529138
        5         H     -0.516266      0.028566      0.999097
    <BLANKLINE>

    # Export both the atoms and mol (i.e. the lattice)
    >>> mol_list3 = hdf5_mol.to_molecules(keys=["atoms", "mol"])
    >>> print(mol_list3[5])  # doctest: +NORMALIZE_WHITESPACE
      Atoms:
        1         C     -0.063899      0.013017     -0.012822
        2         O      1.326996     -0.105806      0.098365
        3         H     -0.290404      0.952594     -0.555502
        4         H     -0.456427     -0.888371     -0.529138
        5         H     -0.516266      0.028566      0.999097
      Lattice:
            0.0000000000     0.0000000000     0.0000000000
            0.0000000000     0.0000000000     0.0000000000
            0.0000000000     0.0000000000     0.0000000000
    <BLANKLINE>

Storing Additional Properties
-----------------------------
:class:`HDF5Mol` is designed with subclassing in mind for cases when, for example,
one wants to store additional molecular properties.
The subclassing consists of three general steps:

1.  Create a new datatype for the property of interest.
    In the example below the :attr:`HDF5Mol.mol` dtype is expanded with
    a new "charge" field for containing the molecular charge.
2.  Set the new dtype in the :attr:`HDF5Mol.INFO["mol"]<InfoTuple>` named tuple.
    Changes to :attr:`HDF5Mol.atoms` and :attr:`HDF5Mol.bonds` should, respectively,
    be set under the ``"atoms"`` and ``"bonds"`` keys.
3.  Update the respective ``get`` and ``set`` methods
    (*e.g.* :meth:`~HDF5Mol.get_mol_info` and :meth:`~HDF5Mol.set_mol_info`).
    These methods are, as the name implies, responsible for reading/writting the relevant
    data to/from any passed molecules.

.. Warning
    Altering the :attr:`~InfoTuple.vlen` and :attr:`~InfoTuple.ndim` fields of pre-defined
    attributes (_e.g._ :attr:`~HDF5Mol.atoms` or :attr:`~HDF5Mol.scale`)
    is ill-advised and can lead to undefined behavior.

.. code-block:: python
    :emphasize-lines: 4, 11, 15

    >>> import numpy as np
    >>> from scm.plams import HDF5Mol

    >>> # Step 1: Define the new data type
    >>> mol_dtype = np.dtype([
    ...     ('lattice', 'float64', (3, 3)),
    ...     ('charge', 'float64')
    ... ])

    >>> class SubClass(HDF5Mol):
    ...     # Step 2: Set the new data type
    ...     INFO = HDF5Mol.INFO.copy()
    ...     INFO["mol"] = INFO["mol"]._replace(dtype=mol_dtype)
    ...
    ...     # Step 3: Update the `get` and `set` methods
    ...     @staticmethod
    ...     def get_mol_info(mol, dtype):
    ...         lattice = mol.lattice if np.any(mol.lattice) else 0
    ...         charge = mol.properties.get('charge', 0)
    ...         return np.array((lattice, charge), dtype=dtype)
    ...
    ...     @staticmethod
    ...     def set_mol_info(mol, data):
    ...         lattice, charge = data.item()
    ...         mol.lattice = lattice.tolist()
    ...         mol.properties.charge = charge

.. testcleanup:: python

    >>> from scm.plams.testing_utils import MOL_TUPLE

    >>> hdf5_mol = SubClass.from_molecules(MOL_TUPLE)
    >>> assert hdf5_mol.mol.dtype == mol_dtype


Another example is the :class:`HDF5Pdb` subclass, which has its :attr:`~HDF5Pdb.atoms`
modified such that it can store all necessary information for constructing protein databank files
(*i.e.* the .pdb format; see :func:`~scm.plams.interfaces.molecule.rdkit.readpdb` and
:func:`~scm.plams.interfaces.molecule.rdkit.writepdb`).


Adding or Removing Attributes
-----------------------------
.. code-block:: python
    :emphasize-lines: 4, 16, 19, 24, 30, 37

    >>> import numpy as np
    >>> from scm.plams import HDF5Mol, InfoTuple

    >>> # 1: Define the data types of the new attributes
    >>> solv_dtype = np.dtype('float64')

    >>> class HDF5Solv(HDF5Mol):
    ...     # 2: Create a slot for the new attribute
    ...     __slots__ = ['solv']
    ...
    ...     # 3: Update the `INFO` dictionary
    ...     INFO = HDF5Mol.INFO.copy()
    ...     INFO['solv'] = InfoTuple('solv', solv_dtype, ndim=2)
    ...
    ...     # 4: Update the signature of `__init__()`
    ...     def __init__(self, atoms=None, bonds=None, mol=None, scale=None,  solv=None,
    ...                  validate=True, copy=True):
    ...         array_dict = {'atoms': atoms, 'bonds': bonds, 'mol': mol, 'scale': scale, 'solv': solv}
    ...         self._init(array_dict, validate=validate, copy=copy)
    ...
    ...     # 5: Update the signature of `create_group()`
    ...     def create_group(self, file, name, atoms=None, bonds=None, mol=None, scale=None, solv=None,
    ...                      create_group=True, **kwargs):
    ...         dset_dict = {'atoms': atoms, 'bonds': bonds, 'mol': mol, 'scale': scale, 'solv': solv}
    ...         self._create_group(file, name dset_dict, create_group, **kwargs)
    ...
    ...     # 6: Update the signature of `from_molecules()`
    ...     @classmethod
    ...     def from_molecules(cls, mol_list, scale=None, solv=None, keys=None):
    ...         array_dict = {'scale': scale, 'solv': solv}
    ...         return cls._from_molecules(mol_list, array_dict, keys=keys)

The new ``HDF5Solv`` class in action.

.. code-block:: python

    >>> solv_shape = (len(mol_list), 10)
    >>> solv = np.zeros(solv_shape)

    >>> hdf5_solv = HDF5Solv.from_molecules(mol_list, solv=solv)
    >>> print(hdf5_solv)
    HDF5Solv(
        atoms = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
        bonds = numpy.ndarray(..., shape=(23,), dtype=h5py.vlen_dtype(...)),
        mol   = numpy.ndarray(..., shape=(23,), dtype=...),
        scale = numpy.ndarray(..., shape=(23,), dtype=int64),
        solv  = numpy.ndarray(..., shape=(23, 10), dtype=float64)
    )

.. testcleanup:: python

    >>> with h5py.File(filename, 'r+') as f:
    ...     grp = hdf5_solv.create_group(f, HDF5Solv.__name__)
    ...     hdf5_solv.to_hdf5(grp, np.s_[0:23])
    ...     hdf5_solv2 = HDF5Solv.from_hdf5(grp)
    ...     assert hdf5_solv == hdf5_solv2

    >>> if os.path.isfile(filename):
    ...     os.remove(filename)


Index
-----
.. autosummary::
    HDF5Mol

Index: Miscellaneous Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    HDF5Mol.keys
    HDF5Mol.values
    HDF5Mol.items
    HDF5Mol.__getitem__
    HDF5Mol.__len__
    HDF5Mol.copy
    HDF5Mol.concatenate

Index: Molecule Interconversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    HDF5Mol.from_molecules
    HDF5Mol.to_molecules

Index: HDF5 Interconversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    HDF5Mol.create_group
    HDF5Mol.from_hdf5
    HDF5Mol.to_hdf5

Index: Set Operations
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    HDF5Mol.intersection
    HDF5Mol.difference
    HDF5Mol.symmetric_difference
    HDF5Mol.union

Index: Subclassing
~~~~~~~~~~~~~~~~~~
.. autosummary::
    HDF5Pdb
    HDF5Mol.set_atom_info
    HDF5Mol.set_bond_info
    HDF5Mol.set_mol_info
    HDF5Mol.get_atom_info
    HDF5Mol.get_bond_info
    HDF5Mol.get_mol_info

Index: InfoTuple
~~~~~~~~~~~~~~~~
.. autosummary::
    InfoTuple

Index: Data Type
~~~~~~~~~~~~~~~~
.. autosummary::
    ATOMS_DTYPE
    ATOMS_DTYPE_PDB
    BONDS_DTYPE
    MOL_DTYPE
    SCALE_DTYPE

API
---
.. autoclass:: HDF5Mol
    :members: atoms, bonds, mol, scale, INFO, __init__

API: Miscellaneous Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: HDF5Mol.keys
.. automethod:: HDF5Mol.values
.. automethod:: HDF5Mol.items
.. automethod:: HDF5Mol.__getitem__
.. automethod:: HDF5Mol.__len__
.. automethod:: HDF5Mol.copy
.. automethod:: HDF5Mol.concatenate

API: Molecule Interconversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: HDF5Mol.from_molecules
.. automethod:: HDF5Mol.to_molecules

API: HDF5 Interconversion
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: HDF5Mol.create_group
.. automethod:: HDF5Mol.from_hdf5
.. automethod:: HDF5Mol.to_hdf5

API: Set Operations
~~~~~~~~~~~~~~~~~~~
.. automethod:: HDF5Mol.intersection
.. automethod:: HDF5Mol.difference
.. automethod:: HDF5Mol.symmetric_difference
.. automethod:: HDF5Mol.union

API: Subclassing
~~~~~~~~~~~~~~~~
.. autoclass:: HDF5Pdb
    :members: atoms, bonds, mol, scale, __init__

Important
---------
The six methods described below are intended as hooks for the purpose of subclassing.
They should, generally speaking, not be called on their own.


.. automethod:: HDF5Mol.set_atom_info
.. automethod:: HDF5Mol.set_bond_info
.. automethod:: HDF5Mol.set_mol_info
.. automethod:: HDF5Mol.get_atom_info
.. automethod:: HDF5Mol.get_bond_info
.. automethod:: HDF5Mol.get_mol_info

API: InfoTuple
~~~~~~~~~~~~~~
.. autoclass:: InfoTuple
    :members:

API: Data Type
~~~~~~~~~~~~~~
.. autodata:: ATOMS_DTYPE
    :annotation: = numpy.dtype(...)

    The dtype of :attr:`HDF5Mol.atoms`.

    .. code-block:: python

        >>> import numpy as np

        >>> ATOMS_DTYPE = np.dtype([
        ...     ('symbol', 'S2'),
        ...     ('x', 'float64'),
        ...     ('y', 'float64'),
        ...     ('z', 'float64')
        ... ])

    .. testcleanup:: python

        >>> from scm.plams import ATOMS_DTYPE as ref
        >>> assert ATOMS_DTYPE == ref

.. autodata:: ATOMS_DTYPE_PDB
    :annotation: = numpy.dtype(...)

    The dtype of :attr:`HDF5Pdb.atoms`.

    .. code-block:: python

        >>> import numpy as np

        >>> ATOMS_DTYPE_PDB = np.dtype([
        ...     ('IsHeteroAtom', 'bool'),
        ...     ('SerialNumber', 'int16'),
        ...     ('Name', 'S4'),
        ...     ('AltLoc', 'S1'),
        ...     ('ResidueName', 'S3'),
        ...     ('ChainId', 'S1'),
        ...     ('ResidueNumber', 'int16'),
        ...     ('InsertionCode', 'S1'),
        ...     ('x', 'float64'),
        ...     ('y', 'float64'),
        ...     ('z', 'float64'),
        ...     ('Occupancy', 'float64'),
        ...     ('TempFactor', 'float64'),
        ...     ('symbol', 'S4'),
        ...     ('charge', 'int8')
        ... ])

    .. testcleanup:: python

        >>> from scm.plams import ATOMS_DTYPE_PDB as ref
        >>> assert ATOMS_DTYPE_PDB == ref

.. autodata:: BONDS_DTYPE
    :annotation: = numpy.dtype(...)

    The dtype of :attr:`HDF5Mol.bonds`.

    .. code-block:: python

        >>> import numpy as np

        >>> BONDS_DTYPE = np.dtype([
        ...     ('atom1', 'int64'),
        ...     ('atom2', 'int64'),
        ...     ('order', 'float64')
        ... ])

    .. testcleanup:: python

        >>> from scm.plams import BONDS_DTYPE as ref
        >>> assert BONDS_DTYPE == ref

.. autodata:: MOL_DTYPE
    :annotation: = numpy.dtype(...)

    The dtype of :attr:`HDF5Mol.mol`.

    .. code-block:: python

        >>> import numpy as np

        >>> MOL_DTYPE = np.dtype([
        ...     ('lattice', 'float64', (3, 3))
        ... ])

    .. testcleanup:: python

        >>> from scm.plams import MOL_DTYPE as ref
        >>> assert MOL_DTYPE == ref

.. autodata:: SCALE_DTYPE
    :annotation: = numpy.dtype(...)

    The dtype of :attr:`HDF5Mol.scale`.

    .. code-block:: python

        >>> import numpy as np

        >>> SCALE_DTYPE = np.dtype('int64')

    .. testcleanup:: python

        >>> from scm.plams import SCALE_DTYPE as ref
        >>> assert SCALE_DTYPE == ref

"""
