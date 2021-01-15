Trajectories
-------------------------

In this chapter the PLAMS module responsible for reading and writing molecular trajectories is presented.
The main output of a molecular simulation is often a series of conformations of a molecular system,
as successively produced by the software over time.
Example applications that produce such a trajectory are molecular dynamics (MD) simulations, 
and geometry optimizations.
The conformations consist of coordinates for all the atoms in the system, and - depending on the format -
information on its periodicity, atomic elements and the location of covalent bonds.
The most common file format storing these trajectories is the XYZ format (``xyz``), which is a human readable format
that stores elements and coordinates for every conformation.
Other examples are the compact binary DCD format (``dcd``), which stores only coordinates and periodic lattice vectors,
and the binary RKF format (``rkf``) which is the native format of the Amsterdam Modeling Suite, and is able to store
any desired information on the system.
The vast majority of molecular modeling software packages will be able to produce the trajectory in at least one of these formats.

Through the ``Molecule`` object PLAMS provides the possibility to analyze and alter single conformations of a molecular system.
The trajectories module provides the possibility to extract such a conformation from a trajectory file 
and load it into a ``Molecule`` object,
as well as write (altered) conformations to a new trajectory file of one in the three formats described above.


.. toctree::

    xyz
    rkf
    dcd
