import os

_HERE = os.path.abspath(os.path.dirname(__file__))

#: A tuple with the absolute paths to all .pdb files in this directory.
PDB_TUPLE = tuple(
    os.path.join(_HERE, f) for f in os.listdir(_HERE) if f.endswith('pdb')
)
