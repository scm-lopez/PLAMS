from scm.plams import Molecule, label

m1 = Molecule('xyz/CO_flat4_1.xyz')
m2 = Molecule('xyz/CO_flat4_2.xyz')
m3 = Molecule('xyz/CO_flat4_3.xyz')
m4 = Molecule('xyz/CO_flat4_4.xyz')


def testYES():
    for i in range(4): assert m1.label(i) == m2.label(i)
    for i in range(4): assert m3.label(i) == m4.label(i)

def testNO():
    assert m1.label(4) != m2.label(4)
    assert m3.label(4) != m4.label(4)



