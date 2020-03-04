from scm.plams import Molecule
from scm.plams.mol.identify import label

m1 = Molecule('xyz/RS1.xyz')
m2 = Molecule('xyz/RS2.xyz')

def testYES():
    for i in range(3): assert m1.label(i) == m2.label(i)

def testNO():
    for i in range(3,5): assert m1.label(i) != m2.label(i)