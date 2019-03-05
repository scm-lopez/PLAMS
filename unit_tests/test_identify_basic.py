from scm.plams import Molecule, PT

PT.set_connectors('Mg', 4)

m1 = Molecule('xyz/chlorophyl1.xyz')
m2 = Molecule('xyz/chlorophyl2.xyz')

def testYES():
    for i in range(2): assert m1.label(i) == m2.label(i)

def testNO():
    for i in range(2,5): assert m1.label(i) != m2.label(i)