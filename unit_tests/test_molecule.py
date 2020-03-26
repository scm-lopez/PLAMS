from pathlib import Path

import numpy as np
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
