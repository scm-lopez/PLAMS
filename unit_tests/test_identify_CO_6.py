from scm.plams import Molecule, label

m1 = Molecule('xyz/CO_6_1.xyz')
m2 = Molecule('xyz/CO_6_2.xyz')

n1 = Molecule('xyz/CO_6_3.xyz')
n2 = Molecule('xyz/CO_6_4.xyz')
n3 = Molecule('xyz/CO_6_5.xyz')

def testYES():
    for i in range(4): assert m1.label(i) == m2.label(i)
    for i in range(4): assert n1.label(i) == n2.label(i) == n3.label(i)

def testNO():
    assert m1.label(4) != m2.label(4)
    assert n1.label(4) != n2.label(4)
    assert n1.label(4) != n3.label(4)
    assert n2.label(4) != n3.label(4)



