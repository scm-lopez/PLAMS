"""Tests for :mod:`scm.plams.h5py_tools`."""

import os
import copy
import dill as pickle
from shutil import copyfile
from tempfile import TemporaryDirectory

import pytest
import numpy as np

try:
    import h5py
    from scm.plams import Molecule, HDF5Mol, HDF5Pdb
    from scm.plams.testing_utils import assert_mol_eq, HDF5_MOL, HDF5_TMP, HDF5_READ, MOL_TUPLE
    H5PY_EX = None
except ImportError as ex1:
    H5PY_EX = ex1

try:
    import rdkit
    from scm.plams import readpdb, writepdb
    RDKIT_EX = None
except ImportError as ex2:
    RDKIT_EX = ex2

    def readpdb(file):
        ret = Molecule(file)
        ret.guess_bonds()
        return ret


@pytest.mark.skipif(H5PY_EX, reason=str(H5PY_EX))
def test_pickle() -> None:
    """Test :meth:`HDF5Mol.__reduce__`."""
    dumps = pickle.dumps(HDF5_MOL)
    loads = pickle.loads(dumps)
    assert loads == HDF5_MOL


@pytest.mark.skipif(H5PY_EX, reason=str(H5PY_EX))
def test_eq() -> None:
    """Test :meth:`HDF5Mol.__eq__`."""
    hdf5_mol1 = HDF5_MOL[:]
    hdf5_mol2 = HDF5_MOL[0]

    assert hdf5_mol1 == HDF5_MOL
    assert hdf5_mol2 != HDF5_MOL
    assert 1 != HDF5_MOL


@pytest.mark.skipif(H5PY_EX, reason=str(H5PY_EX))
def test_copy() -> None:
    """Test :meth:`HDF5Mol.__copy__`, :meth:`HDF5Mol.copy` and :meth:`HDF5Mol.__deepcopy__`."""
    hdf5_mol1 = copy.copy(HDF5_MOL)
    hdf5_mol2 = copy.deepcopy(HDF5_MOL)
    hdf5_mol3 = HDF5_MOL.copy()

    assert HDF5_MOL == hdf5_mol1
    assert HDF5_MOL == hdf5_mol2
    assert HDF5_MOL == hdf5_mol3


@pytest.mark.skipif(H5PY_EX, reason=str(H5PY_EX))
def test_to_molecules() -> None:
    """Test :meth:`HDF5Mol.to_molecules`."""
    with h5py.File(HDF5_READ, 'r') as f:
        grp = f['ligand']
        hdf5_mol = HDF5Mol.from_hdf5(grp)

    mol_list = hdf5_mol.to_molecules()

    for mol, ref in zip(mol_list, MOL_TUPLE):
        assert_mol_eq(mol, ref)


@pytest.mark.skipif(H5PY_EX, reason=str(H5PY_EX))
def test_index_other() -> None:
    copyfile(HDF5_READ, HDF5_TMP)

    try:
        with h5py.File(HDF5_TMP, 'r+') as f:
            scale = f['ligand/scale']
            group = HDF5_MOL.create_group(f, 'test', scale=scale)
            dset = group['atoms']
            assert scale == dset.dims[0]['scale']
    finally:
        if os.path.isfile(HDF5_TMP):
            os.remove(HDF5_TMP)

@pytest.mark.skipif(H5PY_EX, reason=str(H5PY_EX))
@pytest.mark.skipif(RDKIT_EX, reason=str(RDKIT_EX))
def test_pdb() -> None:
    try:
        pdb = HDF5Pdb.from_molecules(MOL_TUPLE)
        with h5py.File(HDF5_TMP, 'a') as f:
            grp = pdb.create_group(f, 'ligand')
            pdb.to_hdf5(grp, None)
            mol_list = HDF5Pdb.from_hdf5(grp).to_molecules()

        for mol, ref in zip(mol_list, MOL_TUPLE):
            assert_mol_eq(mol, ref)

        with TemporaryDirectory() as tmp_dir:
            for i, mol in enumerate(mol_list):
                filename = os.path.join(tmp_dir, f'mol{i}.pdb')
                writepdb(mol, filename)

    finally:
        if os.path.isfile(HDF5_TMP):
            os.remove(HDF5_TMP)
