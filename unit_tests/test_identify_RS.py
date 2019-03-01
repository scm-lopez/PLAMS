from scm.plams import Molecule, label

m1 = Molecule('xyz/RS1.xyz')
m2 = Molecule('xyz/RS2.xyz')

def testRS0(): assert label(m1,0) == label(m2,0)
def testRS1(): assert label(m1,1) == label(m2,1)
def testRS2(): assert label(m1,2) != label(m2,2)
def testRS3(): assert label(m1,3) != label(m2,3)