from pathlib import Path

from scm.plams import Molecule, PT

PATH = Path('.') / 'xyz'

PT.set_connectors('Mg', 4)
m1 = Molecule(PATH / 'chlorophyl1.xyz')
m2 = Molecule(PATH / 'chlorophyl2.xyz')


def testYES():
    for i in range(2): assert m1.label(i) == m2.label(i)


def testNO():
    for i in range(2,5): assert m1.label(i) != m2.label(i)
