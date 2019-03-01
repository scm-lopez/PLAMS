from scm.plams import Molecule, label

m1 = Molecule('xyz/EZ1.xyz')
m2 = Molecule('xyz/EZ2.xyz')

def testEZ0(): assert label(m1,0) == label(m2,0)
def testEZ1(): assert label(m1,1) == label(m2,1)
def testEZ2(): assert label(m1,2) != label(m2,2)
def testEZ3(): assert label(m1,3) != label(m2,3)