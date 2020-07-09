import os
from pathlib import Path

import numpy as np
try:
    import dill as pickle
except ImportError:
    import pickle

from scm.plams import Molecule, Atom, MoleculeError


PATH = Path('unit_tests') / 'xyz'
BENZENE = Molecule(PATH / 'benzene.xyz')
BENZENE.guess_bonds()


def test_index():
    """Test :meth:`Molecule.index`."""
    atom = BENZENE[1]
    bond = BENZENE[1, 2]
    atom_test = Atom(coords=[0, 0, 0], symbol='H')

    assert BENZENE.index(atom) == 1
    assert BENZENE.index(bond) == (1, 2)

    try:
        BENZENE.index(None)  # None is of invalid type
    except MoleculeError:
        pass
    else:
        raise AssertionError("'BENZENE.index(None)' failed to raise a 'MoleculeError'")

    try:
        BENZENE.index(atom_test)  # atom_test is not in BENZENE
    except MoleculeError:
        pass
    else:
        raise AssertionError("'BENZENE.index(atom_test)' failed to raise a 'MoleculeError'")


def test_set_integer_bonds():
    """Test :meth:`Molecule.set_integer_bonds`."""
    ref1 = np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1, 1, 1], dtype=float)
    ref2 = np.array([1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1], dtype=float)

    benzene = BENZENE.copy()
    np.testing.assert_array_equal([b.order for b in benzene.bonds], ref1)

    benzene.set_integer_bonds()
    np.testing.assert_array_equal([b.order for b in benzene.bonds], ref2)


def test_round_coords():
    """Test :meth:`Molecule.round_coords`."""
    benzene = BENZENE.copy()
    ref1 = np.array([[ 1., -1.,  0.],
                     [ 1.,  1.,  0.],
                     [ 0.,  1.,  0.],
                     [-1.,  1.,  0.],
                     [-1., -1.,  0.],
                     [ 0., -1.,  0.],
                     [ 2., -1.,  0.],
                     [ 2.,  1.,  0.],
                     [ 0.,  2.,  0.],
                     [-2.,  1.,  0.],
                     [-2., -1.,  0.],
                     [ 0., -2.,  0.]])
    ref2 = np.array([[ 1.19, -0.69,  0.  ],
                     [ 1.19,  0.69,  0.  ],
                     [ 0.  ,  1.38,  0.  ],
                     [-1.19,  0.69,  0.  ],
                     [-1.19, -0.69,  0.  ],
                     [-0.  , -1.38,  0.  ],
                     [ 2.13, -1.23, -0.  ],
                     [ 2.13,  1.23, -0.  ],
                     [ 0.  ,  2.46, -0.  ],
                     [-2.13,  1.23, -0.  ],
                     [-2.13, -1.23, -0.  ],
                     [-0.  , -2.46, -0.  ]])

    benzene2 = round(benzene)
    np.testing.assert_array_equal(benzene2, ref1)

    benzene.round_coords(decimals=2)
    np.testing.assert_allclose(benzene, ref2)


IMMUTABLE_TYPE = (int, float, tuple)
ATTR_EXCLUDE = frozenset({'mol', 'bonds', 'atoms', 'atom1', 'atom2'})


def _compare_attrs(obj1, obj2, eval_eq=True):
    assert obj1 is not obj2

    for name, attr in vars(obj1).items():
        if name in ATTR_EXCLUDE:
            continue

        attr_ref = getattr(obj2, name)
        if eval_eq:
            assert attr == attr_ref
        if not isinstance(attr, IMMUTABLE_TYPE):
            assert attr is not attr_ref


def test_copy(mol_list=None):
    """Tests for :meth:`Molecule.copy`, :meth:`Molecule.__copy__` and :meth:`Molecule.__deepcopy__`."""
    if mol_list is None:
        mol_list = [
            BENZENE.copy(),
            BENZENE.__copy__(),
            BENZENE.__deepcopy__({})
        ]

    for mol in mol_list:
        _compare_attrs(mol, BENZENE, eval_eq=False)

        for at, at_ref in zip(mol.atoms, BENZENE.atoms):
            _compare_attrs(at, at_ref)

        for bond, bond_ref in zip(mol.bonds, BENZENE.bonds):
            _compare_attrs(bond, bond_ref)


def test_set_get_state():
    """Tests for :meth:`Molecule.__setstate__` and :meth:`Molecule.__getstate__`."""
    mol = BENZENE.copy()
    dill_new = PATH / 'benzene_new.dill'

    # An old pickle file created before the introduction of Molecule.__getstate__
    dill_old = PATH / 'benzene_old.dill'

    try:
        with open(dill_new, 'wb') as f:
            pickle.dump(mol, f)
        with open(dill_new, 'rb') as f:
            mol_new = pickle.load(f)

        with open(dill_old, 'rb') as f:
            mol_old = pickle.load(f)

        test_copy(mol_list=[mol_new, mol_old])
    finally:
        os.remove(dill_new) if os.path.isfile(dill_new) else None